
use rand::{ random, Rng };
use ray_tracing::utils::{ Colour, Sphere, Point3, Vec3 };
use ray_tracing::material::{ CheckerTexture, ImageTexture, Dielectric, Lambertian, Metal };
use ray_tracing::camera::{ Camera, CameraOpts };
use std::rc::Rc;

#[allow(dead_code)]
fn random_spheres() -> (Vec<Sphere>, CameraOpts) {
    let material_ground = Rc::new(Lambertian::from_colour(Colour::new(0.5, 0.5, 0.5)));
    let mut world = Vec::new();
    world.push( Sphere::new(Point3::new(0.0, -1000.0, 0.0), 1000.0, material_ground.clone()) );
    for a in -11..11 {
        for b in -11..11 {
            let center = Point3::new(a as f64 + 0.9 * random::<f64>(), 0.2, b as f64 + 0.9 * random::<f64>());
            if (center - &Point3::new(4.0, 0.2, 0.0)).length() > 0.9 {
                let rand : f64 = random();
                if rand < 0.8 {
                    let albedo = Colour::random() * &Colour::random();
                    world.push(Sphere::new(center, 0.2, Rc::new(Lambertian::from_colour(albedo))));
                } else if rand < 0.95 {
                    let albedo = Vec3::random_in_range(0.5, 1.0);
                    let fuzz = (rand::thread_rng()).gen_range(0.0..0.5);
                    world.push(Sphere::new(center, 0.2, Rc::new(Metal::new(albedo, fuzz))));
                } else {
                    world.push(Sphere::new(center, 0.2, Rc::new(Dielectric::new(1.5))));
                }
            }
        }
    }

    world.push(Sphere::new(Point3::new(-4.0, 1.0, 0.0), 1.0, Rc::new(Lambertian::from_colour(Colour::new(0.4, 0.2, 0.1)))));
    world.push(Sphere::new(Point3::new(0.0, 1.0, 0.0), 1.0, Rc::new(Dielectric::new(1.5))));
    world.push(Sphere::new(Point3::new(4.0, 1.0, 0.0), 1.0, Rc::new(Metal::new(Colour::new(0.7, 0.6, 0.5), 0.0))));

    let camera_opts = CameraOpts { 
        aspect_ratio: 16.0/9.0, 
        image_width: 1200, 
        samples_per_pixel: 200,
        max_depth: 50,
        vfov: 20.0,
        look_from: Point3::new(13.0, 2.0, 3.0),
        look_at: Point3::new(0.0, 0.0, 0.0),
        vup: Vec3::new(0.0, 1.0, 0.0),
        defocus_angle: 0.0,
        focus_distance: 10.0
    };
    (world, camera_opts)
} 

#[allow(dead_code)]
fn exercise() -> (Vec<Sphere>, CameraOpts) {
    let material_ground = Rc::new(Lambertian::from_colour(Colour::new(0.8, 0.8, 0.0)));
    let material_center = Rc::new(Lambertian::from_colour(Colour::new(0.1, 0.2, 0.5)));
    let material_left = Rc::new(Dielectric::new(1.5));
    let material_right =  Rc::new(Metal::new(Colour::new(0.8, 0.6, 0.2), 0.1));
    let world = vec! { 
        Sphere::new(Point3::new(0.0, 0.0, -1.0), 0.5, material_center.clone()),
        Sphere::new(Point3::new(-1.0, 0.0, -1.0), 0.5, material_left.clone()),
        Sphere::new(Point3::new(-1.0, 0.0, -1.0), -0.4, material_left.clone()),
        Sphere::new(Point3::new(1.0, 0.0, -1.0), 0.5, material_right.clone()),
        Sphere::new(Point3::new(0.0, -100.5, -1.0), 100.0, material_ground.clone())
    };

    let camera_opts = CameraOpts { 
        aspect_ratio: 16.0/9.0, 
        image_width: 400, 
        samples_per_pixel: 200,
        max_depth: 10,
        vfov: 20.0,
        look_from: Point3::new(-2.0, 2.0, 1.0),
        look_at: Point3::new(0.0, 0.0, -1.0),
        vup: Vec3::new(0.0, 1.0, 0.0),
        defocus_angle: 0.0,
        focus_distance: 10.0
    };
    (world, camera_opts)
}

#[allow(dead_code)]
fn two_spheres() -> (Vec<Sphere>, CameraOpts) {
    let checker = Rc::new(CheckerTexture::new_with_colours(0.8, Colour::new(0.2, 0.3, 0.1), Colour::new(0.9, 0.9, 0.9)));
    let world = vec! {
        Sphere::new(Point3::new(0.0, -10.0, 0.0), 10.0, Rc::new(Lambertian::new(checker.clone()))),
        Sphere::new(Point3::new(0.0, 10.0, 0.0), 10.0, Rc::new(Lambertian::new(checker.clone())))
    };

    let camera_opts = CameraOpts { 
        aspect_ratio: 16.0/9.0, 
        image_width: 400, 
        samples_per_pixel: 20,
        max_depth: 50,
        vfov: 20.0,
        look_from: Point3::new(13.0,2.0,3.0),
        look_at: Point3::new(0.0, 0.0, 0.0),
        vup: Vec3::new(0.0, 1.0, 0.0),
        defocus_angle: 0.0,
        focus_distance: 10.0
    };
    (world, camera_opts)
}

#[allow(dead_code)]
fn earth() -> (Vec<Sphere>, CameraOpts) {
    let earth_image = image::open("earthmap.jpg").unwrap();
    let earth_texture = Rc::new(ImageTexture::new(earth_image));
    let earth_material = Rc::new(Lambertian::new(earth_texture));

    let material_right =  Rc::new(Metal::new(Colour::new(0.8, 0.6, 0.2), 0.0));
    let world = vec! {
        Sphere::new(Point3::new(0.0, 0.0, 0.0), 2.0, earth_material.clone()),
        Sphere::new(Point3::new(5.0, 2.0, -1.0), 3.0, material_right.clone())
    };

    let camera_opts = CameraOpts { 
        aspect_ratio: 16.0/9.0, 
        image_width: 400, 
        samples_per_pixel: 50,
        max_depth: 10,
        vfov: 30.0,
        look_from: Point3::new(0.0,0.0,12.0),
        look_at: Point3::new(0.0, 0.0, 0.0),
        vup: Vec3::new(0.0, 1.0, 0.0),
        defocus_angle: 0.0,
        focus_distance: 10.0
    };
    (world, camera_opts)
}
fn main() {
    let (world, camera_opts)  = earth();
    let camera = Camera::initialize(camera_opts);
    camera.render(&world);    
}
