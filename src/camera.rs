use crate::utils::{ Vec3, Ray, Hittable, Colour, write_colour };
use ordered_float::OrderedFloat;
use rand::Rng;

pub struct CameraOpts {
    pub aspect_ratio : f64,
    pub image_width : i32,
    pub samples_per_pixel : i32,
    pub max_depth: u8,
    pub vfov: f64,
    pub look_from: Vec3,
    pub look_at: Vec3,
    pub vup: Vec3,
    pub defocus_angle: f64,
    pub focus_distance: f64
}

pub struct Camera {
    samples_per_pixel: i32,
    max_depth: u8,    
    image_width : i32,
    image_height: i32,
    camera_center: Vec3,
    pixel_delta_u: Vec3,
    pixel_delta_v: Vec3,
    pixel00_location: Vec3,
    defocus_disk_u: Vec3,
    defocus_disk_v: Vec3,
    defocus_angle: f64
}

fn degrees_to_radians(degrees: f64) -> f64 {
    degrees * std::f64::consts::PI / 180.0
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
                    pixel_colour = pixel_colour + &Camera::ray_colour(&r, self.max_depth, world);
                }
                write_colour(&pixel_colour, self.samples_per_pixel);
            }
        }
    }

    fn defocus_disk_sample(&self) -> Vec3 {
        let p = Vec3::random_in_unit_disk();
        return self.camera_center + &(self.defocus_disk_u * p.x()) + &(self.defocus_disk_v * p.y());
    }

    fn get_ray(&self, i: i32, j: i32) -> Ray {
        let pixel_center = self.pixel00_location + &(self.pixel_delta_u * (j as f64)) + &(self.pixel_delta_v * (i as f64));
        let pixel_sample = pixel_center + &self.pixel_sample_square();
        let ray_direction = pixel_sample - &self.camera_center;
        let ray_origin = if self.defocus_angle <= 0.0 { self.camera_center } else { self.defocus_disk_sample() };
        Ray::new( ray_origin, ray_direction )           
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

        let camera_center : Vec3 = opts.look_from;

        // Determine viewport dimensions.
        let theta = degrees_to_radians(opts.vfov);
        let h = (theta/2.0).tan();
        let viewport_height : f64 = 2.0 * h * opts.focus_distance;
        let viewport_width : f64 = (opts.image_width as f64/image_height as f64) * viewport_height;
        
        // Calculate the u,v,w unit basis vectors for the camera coordinate frame.
        let w = (opts.look_from - &opts.look_at).unit_vector();
        let u = opts.vup.cross(&w).unit_vector();
        let v = w.cross(&u);

        // Calculate the vectors across the horizontal and down the vertical viewport edges.
        let viewpoint_u : Vec3 = u * viewport_width;
        let viewpoint_v : Vec3 = (v * viewport_height) * -1.0;
        let pixel_delta_u : Vec3 = &viewpoint_u / (opts.image_width as f64);
        let pixel_delta_v : Vec3 = &viewpoint_v / (image_height as f64);

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