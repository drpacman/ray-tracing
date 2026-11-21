use crate::utils::{ Vec3, Ray, Colour, write_colour };
use crate::scene::Scene;
use crate::rng::Rng;
use rayon::prelude::*;
use std::io::{self, BufWriter, Write};
use std::sync::atomic::{AtomicUsize, Ordering};

pub struct CameraOpts {
    pub aspect_ratio : f32,
    pub image_width : i32,
    pub samples_per_pixel : i32,
    pub max_depth: u8,
    pub vfov: f32,
    pub look_from: Vec3,
    pub look_at: Vec3,
    pub vup: Vec3,
    pub defocus_angle: f32,
    pub focus_distance: f32
}

pub struct Camera {
    pub(crate) samples_per_pixel: i32,
    pub(crate) max_depth: u8,
    pub(crate) image_width : i32,
    pub(crate) image_height: i32,
    pub(crate) camera_center: Vec3,
    pub(crate) pixel_delta_u: Vec3,
    pub(crate) pixel_delta_v: Vec3,
    pub(crate) pixel00_location: Vec3,
    pub(crate) defocus_disk_u: Vec3,
    pub(crate) defocus_disk_v: Vec3,
    pub(crate) defocus_angle: f32
}

fn degrees_to_radians(degrees: f32) -> f32 {
    degrees * std::f32::consts::PI / 180.0
}

impl Camera {

    // Renders on the CPU and writes a PPM to stdout.
    pub fn render(&self, scene: &Scene) {
        let pixels = self.render_pixels(scene);
        self.write_ppm(&pixels);
    }

    // Returns the sum of all samples for each pixel, row by row.
    fn render_pixels(&self, scene: &Scene) -> Vec<Colour> {
        // Render scanlines in parallel across all cores. Rayon's collect()
        // preserves order, so rows come back top-to-bottom.
        let rows_done = AtomicUsize::new(0);
        let rows: Vec<Vec<Colour>> = (0..self.image_height)
            .into_par_iter()
            .map(|i| {
                let row = (0..self.image_width)
                    .map(|j| {
                        // Each pixel gets its own generator, seeded from its index,
                        // exactly as a GPU thread would.
                        let mut rng = Rng::for_pixel((i * self.image_width + j) as u32, 0);
                        let mut pixel_colour = Colour::new(0.0, 0.0, 0.0);
                        for _ in 0..self.samples_per_pixel {
                            let r = self.get_ray(i, j, &mut rng);
                            pixel_colour = pixel_colour + &Camera::ray_colour(r, self.max_depth, scene, &mut rng);
                        }
                        pixel_colour
                    })
                    .collect();
                let done = rows_done.fetch_add(1, Ordering::Relaxed) + 1;
                eprint!("\rScanlines remaining: {:5}", self.image_height as usize - done);
                row
            })
            .collect();
        eprintln!("\rDone.                     ");
        rows.into_iter().flatten().collect()
    }

    // Writes summed pixel colours as a PPM to stdout, through a buffered, locked stdout.
    pub(crate) fn write_ppm(&self, pixels: &[Colour]) {
        let stdout = io::stdout();
        let mut out = BufWriter::new(stdout.lock());
        writeln!(out, "P3\n{} {}\n255", self.image_width, self.image_height).unwrap();
        for pixel_colour in pixels {
            write_colour(&mut out, pixel_colour, self.samples_per_pixel).unwrap();
        }
        out.flush().unwrap();
    }

    fn defocus_disk_sample(&self, rng: &mut Rng) -> Vec3 {
        let p = Vec3::random_in_unit_disk(rng);
        return self.camera_center + &(self.defocus_disk_u * p.x()) + &(self.defocus_disk_v * p.y());
    }

    fn get_ray(&self, i: i32, j: i32, rng: &mut Rng) -> Ray {
        let pixel_center = self.pixel00_location + &(self.pixel_delta_u * (j as f32)) + &(self.pixel_delta_v * (i as f32));
        let pixel_sample = pixel_center + &self.pixel_sample_square(rng);
        let ray_origin = if self.defocus_angle <= 0.0 { self.camera_center } else { self.defocus_disk_sample(rng) };
        let ray_direction = pixel_sample - &ray_origin;
        Ray::new( ray_origin, ray_direction )
    }

    fn pixel_sample_square(&self, rng: &mut Rng) -> Vec3 {
       let px = -0.5 + rng.next_f32();
       let py = -0.5 + rng.next_f32();
       self.pixel_delta_u * px + &(self.pixel_delta_v * py)
    }

    pub fn initialize(opts: CameraOpts) -> Camera {
        let mut image_height = (opts.image_width as f32 / opts.aspect_ratio) as i32;
        if image_height < 1 {
            image_height = 1;
        }

        let camera_center : Vec3 = opts.look_from;

        // Determine viewport dimensions.
        let theta = degrees_to_radians(opts.vfov);
        let h = (theta/2.0).tan();
        let viewport_height : f32 = 2.0 * h * opts.focus_distance;
        let viewport_width : f32 = (opts.image_width as f32/image_height as f32) * viewport_height;

        // Calculate the u,v,w unit basis vectors for the camera coordinate frame.
        let w = (opts.look_from - &opts.look_at).unit_vector();
        let u = opts.vup.cross(&w).unit_vector();
        let v = w.cross(&u);

        // Calculate the vectors across the horizontal and down the vertical viewport edges.
        let viewpoint_u : Vec3 = u * viewport_width;
        let viewpoint_v : Vec3 = (v * viewport_height) * -1.0;
        let pixel_delta_u : Vec3 = &viewpoint_u / (opts.image_width as f32);
        let pixel_delta_v : Vec3 = &viewpoint_v / (image_height as f32);

        // Calculate the location of the upper left pixel.
        let viewpoint_top_left : Vec3 = camera_center - &(w * opts.focus_distance) - &(&viewpoint_u / 2.0) - &(&viewpoint_v / 2.0);
        let pixel00_location : Vec3 = viewpoint_top_left + &(&pixel_delta_u / 2.0) + &(&pixel_delta_v / 2.0);

        // Calculate the camera defocus disk basis vectors.
        let defocus_radius = opts.focus_distance * (degrees_to_radians(opts.defocus_angle/2.0)).tan();

        Camera {
            samples_per_pixel: opts.samples_per_pixel,
            max_depth: opts.max_depth,
            image_width: opts.image_width,
            image_height: image_height,
            camera_center: camera_center,
            pixel_delta_u: pixel_delta_u,
            pixel_delta_v: pixel_delta_v,
            pixel00_location: pixel00_location,
            // defocus settings
            defocus_disk_u: u * defocus_radius,
            defocus_disk_v: v * defocus_radius,
            defocus_angle: opts.defocus_angle
        }
    }

    // Iterative version of the book's recursive ray_colour. Rather than
    // multiplying each bounce's attenuation into the result on the way back up
    // the call stack, we carry the running product (the "throughput") forward.
    // The result is the same, but there's no recursion, which GPU code can't do.
    fn ray_colour(ray: Ray, max_depth: u8, scene: &Scene, rng: &mut Rng) -> Colour {
        let mut ray = ray;
        let mut throughput = Colour::new(1.0, 1.0, 1.0);
        for _ in 0..max_depth {
            match scene.hit(&ray, 0.001, f32::INFINITY) {
                Some(record) => {
                    let material = &scene.materials[record.material as usize];
                    match material.scatter(&ray, &record, scene, rng) {
                        Some((attenuation, scattered_ray)) => {
                            throughput = throughput * &attenuation;
                            ray = scattered_ray;
                        }
                        None => return Colour::new(0.0, 0.0, 0.0)
                    }
                }
                None => {
                    // background colour
                    let unit_direction = ray.dir.unit_vector();
                    let t = 0.5 * (unit_direction.y() + 1.0);
                    let background = Vec3::new(1.0, 1.0, 1.0) * (1.0 - t) + &(Vec3::new(0.5, 0.7, 1.0) * t);
                    return throughput * &background;
                }
            }
        }
        // Exceeded the bounce limit: no more light is gathered.
        Colour::new(0.0, 0.0, 0.0)
    }
}
