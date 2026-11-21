// Renders a Scene on the GPU with wgpu (which uses Metal on macOS).
//
// The Scene's flat arrays are converted into GPU-layout structs, uploaded once
// as storage buffers, and the compute shader in raytrace.wgsl is dispatched
// once per "pass". Each pass adds a few samples per pixel into an accumulation
// buffer; keeping passes short stops the OS from killing a long-running GPU job.
use crate::camera::Camera;
use crate::material::{ Material, Texture };
use crate::scene::Scene;
use crate::utils::{ Colour, Vec3 };
use bytemuck::{ Pod, Zeroable };
use wgpu::util::DeviceExt;

// --- GPU data layouts --------------------------------------------------------
// WGSL aligns a vec3f to 16 bytes, so each vec3 is followed by a 4-byte scalar
// (or explicit padding) to keep these byte-for-byte identical to the structs in
// raytrace.wgsl. bytemuck's Pod lets us copy slices of them straight into buffers.

#[repr(C)]
#[derive(Copy, Clone, Pod, Zeroable)]
struct GpuParams {
    pixel00: [f32; 3],
    samples_per_pass: u32,
    delta_u: [f32; 3],
    width: u32,
    delta_v: [f32; 3],
    height: u32,
    center: [f32; 3],
    max_depth: u32,
    defocus_u: [f32; 3],
    defocus_angle: f32,
    defocus_v: [f32; 3],
    seed: u32,
}

#[repr(C)]
#[derive(Copy, Clone, Pod, Zeroable)]
struct GpuSphere {
    center: [f32; 3],
    radius: f32,
    material: u32,
    _pad: [u32; 3],
}

#[repr(C)]
#[derive(Copy, Clone, Pod, Zeroable)]
struct GpuMaterial {
    albedo: [f32; 3],
    fuzz: f32,
    kind: u32,
    texture: u32,
    ior: f32,
    _pad: u32,
}

#[repr(C)]
#[derive(Copy, Clone, Pod, Zeroable)]
struct GpuTexture {
    colour0: [f32; 3],
    inv_scale: f32,
    colour1: [f32; 3],
    kind: u32,
    image_offset: u32,
    image_width: u32,
    image_height: u32,
    _pad: u32,
}

#[repr(C)]
#[derive(Copy, Clone, Pod, Zeroable)]
struct GpuBvhNode {
    min: [f32; 3],
    first: u32,
    max: [f32; 3],
    count: u32,
}

// Must match the constants in raytrace.wgsl.
const MATERIAL_LAMBERTIAN: u32 = 0;
const MATERIAL_METAL: u32 = 1;
const MATERIAL_DIELECTRIC: u32 = 2;
const TEXTURE_SOLID: u32 = 0;
const TEXTURE_CHECKER: u32 = 1;
const TEXTURE_IMAGE: u32 = 2;

const WORKGROUP_SIZE: u32 = 8; // @workgroup_size(8, 8) in the shader

fn arr(v: &Vec3) -> [f32; 3] {
    [v.x(), v.y(), v.z()]
}

fn convert_material(m: &Material) -> GpuMaterial {
    let mut g = GpuMaterial::zeroed();
    match *m {
        Material::Lambertian { albedo } => {
            g.kind = MATERIAL_LAMBERTIAN;
            g.texture = albedo;
        }
        Material::Metal { albedo, fuzz } => {
            g.kind = MATERIAL_METAL;
            g.albedo = arr(&albedo);
            g.fuzz = fuzz;
        }
        Material::Dielectric { index_of_refraction } => {
            g.kind = MATERIAL_DIELECTRIC;
            g.ior = index_of_refraction;
        }
    }
    g
}

// Packs all images into one array of RGBA8 pixels, returning it along with
// each image's (offset, width, height) within it.
fn pack_images(images: &[image::RgbImage]) -> (Vec<u32>, Vec<(u32, u32, u32)>) {
    let mut pixels = Vec::new();
    let mut placements = Vec::new();
    for image in images {
        placements.push((pixels.len() as u32, image.width(), image.height()));
        pixels.extend(image.pixels().map(|p| u32::from_le_bytes([p[0], p[1], p[2], 255])));
    }
    if pixels.is_empty() {
        pixels.push(0); // GPU buffers can't be empty
    }
    (pixels, placements)
}

fn convert_texture(t: &Texture, image_placements: &[(u32, u32, u32)]) -> GpuTexture {
    let mut g = GpuTexture::zeroed();
    match *t {
        Texture::Solid(colour) => {
            g.kind = TEXTURE_SOLID;
            g.colour0 = arr(&colour);
        }
        Texture::Checker { inverse_scale, odd, even } => {
            g.kind = TEXTURE_CHECKER;
            g.inv_scale = inverse_scale;
            g.colour0 = arr(&odd);
            g.colour1 = arr(&even);
        }
        Texture::Image { image } => {
            let (offset, width, height) = image_placements[image as usize];
            g.kind = TEXTURE_IMAGE;
            g.image_offset = offset;
            g.image_width = width;
            g.image_height = height;
        }
    }
    g
}

impl Camera {
    fn gpu_params(&self, samples_per_pass: u32, seed: u32) -> GpuParams {
        GpuParams {
            pixel00: arr(&self.pixel00_location),
            samples_per_pass,
            delta_u: arr(&self.pixel_delta_u),
            width: self.image_width as u32,
            delta_v: arr(&self.pixel_delta_v),
            height: self.image_height as u32,
            center: arr(&self.camera_center),
            max_depth: self.max_depth as u32,
            defocus_u: arr(&self.defocus_disk_u),
            defocus_angle: self.defocus_angle,
            defocus_v: arr(&self.defocus_disk_v),
            seed,
        }
    }

    // Renders on the GPU and writes a PPM to stdout.
    // `samples_per_pass` trades GPU job length against per-pass overhead.
    pub fn render_gpu(&self, scene: &Scene, samples_per_pass: u32) {
        let pixels = pollster::block_on(self.render_pixels_gpu(scene, samples_per_pass));
        self.write_ppm(&pixels);
    }

    async fn render_pixels_gpu(&self, scene: &Scene, samples_per_pass: u32) -> Vec<Colour> {
        // --- Device setup ---
        let instance = wgpu::Instance::default();
        let adapter = instance
            .request_adapter(&wgpu::RequestAdapterOptions {
                power_preference: wgpu::PowerPreference::HighPerformance,
                ..Default::default()
            })
            .await
            .expect("no suitable GPU adapter found");
        eprintln!("Using GPU: {}", adapter.get_info().name);
        let (device, queue) = adapter
            .request_device(&wgpu::DeviceDescriptor::default())
            .await
            .expect("failed to create GPU device");

        // --- Scene data ---
        let spheres: Vec<GpuSphere> = scene.spheres.iter()
            .map(|s| GpuSphere { center: arr(&s.center), radius: s.radius, material: s.material, _pad: [0; 3] })
            .collect();
        let materials: Vec<GpuMaterial> = scene.materials.iter().map(convert_material).collect();
        let (image_pixels, image_placements) = pack_images(&scene.images);
        let textures: Vec<GpuTexture> = scene.textures.iter()
            .map(|t| convert_texture(t, &image_placements))
            .collect();
        let bvh: Vec<GpuBvhNode> = scene.bvh.iter()
            .map(|n| GpuBvhNode { min: arr(&n.bbox.min), first: n.first, max: arr(&n.bbox.max), count: n.count })
            .collect();
        assert!(!bvh.is_empty(), "call Scene::build_bvh before rendering");

        let storage = |label: &str, contents: &[u8]| {
            device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
                label: Some(label),
                contents,
                usage: wgpu::BufferUsages::STORAGE,
            })
        };
        let sphere_buffer = storage("spheres", bytemuck::cast_slice(&spheres));
        let material_buffer = storage("materials", bytemuck::cast_slice(&materials));
        let texture_buffer = storage("textures", bytemuck::cast_slice(&textures));
        let bvh_buffer = storage("bvh", bytemuck::cast_slice(&bvh));
        let image_buffer = storage("image_pixels", bytemuck::cast_slice(&image_pixels));

        let params_buffer = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("params"),
            size: std::mem::size_of::<GpuParams>() as u64,
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });

        // One vec4f per pixel; wgpu zero-initialises new buffers.
        let pixel_count = (self.image_width * self.image_height) as u64;
        let accum_size = pixel_count * 16;
        let accum_buffer = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("accum"),
            size: accum_size,
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC,
            mapped_at_creation: false,
        });
        // The GPU-only accum buffer can't be read by the CPU directly, so it's
        // copied into this mappable one at the end.
        let readback_buffer = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("readback"),
            size: accum_size,
            usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });

        // --- Pipeline ---
        let shader = device.create_shader_module(wgpu::include_wgsl!("raytrace.wgsl"));
        let pipeline = device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
            label: Some("raytrace"),
            layout: None, // derive the bind group layout from the shader
            module: &shader,
            entry_point: Some("main"),
            compilation_options: Default::default(),
            cache: None,
        });
        let buffers = [
            &params_buffer, &sphere_buffer, &material_buffer, &texture_buffer,
            &bvh_buffer, &image_buffer, &accum_buffer,
        ];
        let entries: Vec<wgpu::BindGroupEntry> = buffers.iter().enumerate()
            .map(|(binding, buffer)| wgpu::BindGroupEntry {
                binding: binding as u32,
                resource: buffer.as_entire_binding(),
            })
            .collect();
        let bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("raytrace"),
            layout: &pipeline.get_bind_group_layout(0),
            entries: &entries,
        });

        // --- Render, one pass at a time ---
        let workgroups_x = (self.image_width as u32).div_ceil(WORKGROUP_SIZE);
        let workgroups_y = (self.image_height as u32).div_ceil(WORKGROUP_SIZE);
        let total_samples = self.samples_per_pixel as u32;
        let samples_per_pass = samples_per_pass.max(1);
        let mut samples_done = 0;
        let mut pass = 0;
        while samples_done < total_samples {
            let samples = samples_per_pass.min(total_samples - samples_done);
            // A fresh seed per pass so passes take different samples.
            queue.write_buffer(&params_buffer, 0, bytemuck::bytes_of(&self.gpu_params(samples, pass)));

            let mut encoder = device.create_command_encoder(&Default::default());
            {
                let mut compute = encoder.begin_compute_pass(&Default::default());
                compute.set_pipeline(&pipeline);
                compute.set_bind_group(0, &bind_group, &[]);
                compute.dispatch_workgroups(workgroups_x, workgroups_y, 1);
            }
            queue.submit([encoder.finish()]);
            // Wait for this pass, so we report real progress and never queue a
            // long backlog of GPU work.
            device.poll(wgpu::PollType::wait_indefinitely()).expect("GPU poll failed");

            samples_done += samples;
            pass += 1;
            eprint!("\rSamples: {:4} / {}", samples_done, total_samples);
        }
        eprintln!();

        // --- Read the result back ---
        let mut encoder = device.create_command_encoder(&Default::default());
        encoder.copy_buffer_to_buffer(&accum_buffer, 0, &readback_buffer, 0, accum_size);
        queue.submit([encoder.finish()]);

        let slice = readback_buffer.slice(..);
        slice.map_async(wgpu::MapMode::Read, |result| result.expect("failed to map readback buffer"));
        device.poll(wgpu::PollType::wait_indefinitely()).expect("GPU poll failed");
        let data = slice.get_mapped_range().expect("readback buffer not mapped");
        let sums: &[[f32; 4]] = bytemuck::cast_slice(&data);
        sums.iter().map(|c| Colour::new(c[0], c[1], c[2])).collect()
    }
}
