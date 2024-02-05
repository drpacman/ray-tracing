
use ordered_float::OrderedFloat;
use ray_tracing::utils::{ Colour, Sphere, Point3, Ray, Vec3, clamp, Lambertian, Metal, Dielectric };
use ray_tracing::camera::{ Camera, CameraOpts };
use std::rc::Rc;

fn write_colour(pixel_color: &Vec3, samples_per_pixel: i32) {
    let scale = 1.0 / samples_per_pixel as f64;
    let r = (scale * pixel_color.r());
    let g = (scale * pixel_color.g());
    let b = (scale * pixel_color.b());
    let intensity = OrderedFloat(0.0)..OrderedFloat(0.999);
    let x = clamp(r, &intensity);

    println!("{} {} {}", (256.0 * clamp(r, &intensity)) as i32, (256.0 * clamp(g, &intensity)) as i32, (256.0 * clamp(b, &intensity)) as i32);
}

fn main() {
    let material_ground = Rc::new(Lambertian::new(Colour::new(0.8, 0.8, 0.0)));
    let material_center = Rc::new(Lambertian::new(Colour::new(0.1, 0.2, 0.5)));
    let material_left = Rc::new(Dielectric::new(1.5));
    let material_right =  Rc::new(Metal::new(Colour::new(0.8, 0.6, 0.2), 0.1));
    let world = vec! { 
        Sphere::new(Point3::new(0.0, 0.0, -1.0), 0.5, material_center.clone()),
        Sphere::new(Point3::new(-1.0, 0.0, -1.0), 0.5, material_left.clone()),
        Sphere::new(Point3::new(1.0, 0.0, -1.0), 0.5, material_right.clone()),
        Sphere::new(Point3::new(0.0, -100.5, -1.0), 100.0, material_ground.clone())
    };

    let camera_opts = CameraOpts { 
        aspect_ratio: 16.0/9.0, 
        image_width: 400, 
        samples_per_pixel: 5 
    };
    let camera = Camera::initialize(camera_opts);
    camera.render(&world);    
}
