use crate::utils::{ Vec3, Ray, Hittable, Colour, write_colour };
use ordered_float::OrderedFloat;
use rand::Rng;

pub struct CameraOpts {
    pub aspect_ratio : f64,
    pub image_width : i32,
    pub samples_per_pixel : i32
}

pub struct Camera {
    samples_per_pixel: i32,
    aspect_ratio : f64,
    image_width : i32,
    image_height: i32,
    camera_center: Vec3,
    pixel_delta_u: Vec3,
    pixel_delta_v: Vec3,
    pixel00_location: Vec3
}

impl Camera {

    pub fn render(&self, world: &dyn Hittable) {
        // write to provider output stream, not println
        println!("P3\n{} {}\n255", self.image_width, self.image_height);
        for i in 0..self.image_height{
            eprintln!("\rScanlines remaining: {}", self.image_height - i);
            for j in 0..self.image_width {
                let mut pixel_colour = Colour::new(0.0,0.0,0.0);
                for _ in 0..self.samples_per_pixel {
                    let r = self.get_ray(i, j);
                    pixel_colour = pixel_colour + &Camera::ray_colour(&r, 5, world);
                }
                write_colour(&pixel_colour, self.samples_per_pixel);
            }
        }
    }

    fn get_ray(&self, i: i32, j: i32) -> Ray {
        let pixel_center = self.pixel00_location + &(self.pixel_delta_u * (j as f64)) + &(self.pixel_delta_v * (i as f64));
        let pixel_sample = pixel_center + &self.pixel_sample_square();
        let ray_direction = pixel_sample - &self.camera_center;
        Ray::new( self.camera_center, ray_direction )           
    }

    fn pixel_sample_square(&self) -> Vec3 {
       let mut rng = rand::thread_rng();
       let px = -0.5 + rng.gen_range(0.0..1.0);
       let py = -0.5 + rng.gen_range(0.0..1.0);
       self.pixel_delta_u * px + &(self.pixel_delta_v * py)
    }

    pub fn initialize(opts: CameraOpts) -> Camera {
        let mut image_height = (opts.image_width as f64 / opts.aspect_ratio) as i32;
        if image_height < 1 {
            image_height = 1;
        }

        let viewport_height : f64 = 2.0;
        let viewport_width : f64 = (opts.image_width as f64/image_height as f64) * viewport_height;
        let camera_center : Vec3 = Vec3::new(0.0, 0.0, 0.0);
        let focal_length : f64 = 1.0;

        let viewpoint_u : Vec3 = Vec3::new(viewport_width, 0.0, 0.0);
        let viewpoint_v : Vec3 = Vec3::new(0.0, -viewport_height, 0.0);
        let pixel_delta_u : Vec3 = &viewpoint_u / (opts.image_width as f64);
        let pixel_delta_v : Vec3 = &viewpoint_v / (image_height as f64);
        let viewpoint_top_left : Vec3 = camera_center - &(&viewpoint_u / 2.0) - &(&viewpoint_v / 2.0) - &Vec3::new(0.0, 0.0, focal_length);
        let pixel00_location : Vec3 = viewpoint_top_left + &(&pixel_delta_u / 2.0) + &(&pixel_delta_v / 2.0);
        Camera {
            samples_per_pixel: opts.samples_per_pixel,
            aspect_ratio: opts.aspect_ratio,
            image_width: opts.image_width,
            image_height: image_height,
            camera_center: camera_center,
            pixel_delta_u: pixel_delta_u,
            pixel_delta_v: pixel_delta_v,
            pixel00_location: pixel00_location
        }
    }

    fn ray_colour(ray: &Ray, depth: u8, world: &dyn Hittable) -> Colour {
        if depth == 0 {
            return Colour::new(0.0, 0.0, 0.0);
        }

        let range = OrderedFloat(0.001)..OrderedFloat(f64::INFINITY);
        match world.hit(ray, &range) {
            Some(record) => {
                match record.material.scatter(ray, &record) {
                    Some((attenuation, scattered_ray)) => {
                        return attenuation * &Self::ray_colour(&scattered_ray, depth-1, world);
                    }
                    None => return Colour::new(0.0, 0.0, 0.0)
                }
            },
            None  => {           
                // return background colour
                let unit_direction = ray.dir.unit_vector();
                let t = 0.5 * (unit_direction.y() + 1.0);
                Vec3::new(1.0, 1.0, 1.0) * (1.0 - t) + &(Vec3::new(0.5, 0.7, 1.0) * t)
            }   
        }     
    }
}