use crate::bvh::{ self, BvhNode };
use crate::material::{ Material, Texture };
use crate::utils::{ Colour, HitRecord, Point3, Ray, Sphere, Vec3 };

// Everything the renderer needs, as flat arrays that refer to each other by
// u32 index instead of by pointer. This is the layout a GPU wants: each Vec
// becomes one buffer (or, for images, a texture).
pub struct Scene {
    pub spheres: Vec<Sphere>,
    pub materials: Vec<Material>,
    pub textures: Vec<Texture>,
    pub images: Vec<image::RgbImage>,
    pub bvh: Vec<BvhNode>
}

impl Scene {
    pub fn new() -> Self {
        Self { spheres: Vec::new(), materials: Vec::new(), textures: Vec::new(), images: Vec::new(), bvh: Vec::new() }
    }

    // --- Building -------------------------------------------------------
    // Each add_* returns the index of what was added, for use by later calls.

    pub fn add_texture(&mut self, texture: Texture) -> u32 {
        self.textures.push(texture);
        (self.textures.len() - 1) as u32
    }

    pub fn image_texture(&mut self, image: image::RgbImage) -> u32 {
        self.images.push(image);
        let image = (self.images.len() - 1) as u32;
        self.add_texture(Texture::Image { image })
    }

    pub fn add_material(&mut self, material: Material) -> u32 {
        self.materials.push(material);
        (self.materials.len() - 1) as u32
    }

    pub fn lambertian(&mut self, colour: Colour) -> u32 {
        let texture = self.add_texture(Texture::Solid(colour));
        self.lambertian_textured(texture)
    }

    pub fn lambertian_textured(&mut self, texture: u32) -> u32 {
        self.add_material(Material::Lambertian { albedo: texture })
    }

    pub fn metal(&mut self, albedo: Colour, fuzz: f32) -> u32 {
        self.add_material(Material::Metal { albedo, fuzz: fuzz.min(1.0) })
    }

    pub fn dielectric(&mut self, index_of_refraction: f32) -> u32 {
        self.add_material(Material::Dielectric { index_of_refraction })
    }

    pub fn add_sphere(&mut self, center: Point3, radius: f32, material: u32) {
        self.spheres.push(Sphere::new(center, radius, material));
    }

    // Call once all spheres are added. Reorders `spheres`.
    pub fn build_bvh(&mut self) {
        self.bvh = bvh::build(&mut self.spheres);
    }

    // --- Rendering ------------------------------------------------------

    // Finds the closest hit by walking the BVH with an explicit stack
    // (GPU code can't recurse).
    pub fn hit(&self, ray: &Ray, t_min: f32, t_max: f32) -> Option<HitRecord> {
        let mut stack = [0u32; 64];
        let mut stack_len = 1; // start with the root, node 0
        let mut closest_t = t_max;
        let mut closest_sphere = None;

        while stack_len > 0 {
            stack_len -= 1;
            let node = &self.bvh[stack[stack_len] as usize];
            // Using closest_t (not t_max) skips boxes entirely behind the best hit so far.
            if !node.bbox.hit(ray, t_min, closest_t) {
                continue;
            }
            if node.is_leaf() {
                let first = node.first as usize;
                for sphere in &self.spheres[first..first + node.count as usize] {
                    if let Some(t) = sphere.intersect(ray, t_min, closest_t) {
                        closest_t = t;
                        closest_sphere = Some(sphere);
                    }
                }
            } else {
                stack[stack_len] = node.first + 1;
                stack[stack_len + 1] = node.first;
                stack_len += 2;
            }
        }

        closest_sphere.map(|sphere| sphere.hit_record(ray, closest_t))
    }

    pub fn texture_value(&self, texture: u32, u: f32, v: f32, p: &Vec3) -> Colour {
        self.textures[texture as usize].value(u, v, p, &self.images)
    }
}
