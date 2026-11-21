use ray_tracing::utils::{ Colour, Point3, Vec3 };
use ray_tracing::material::Texture;
use ray_tracing::camera::{ Camera, CameraOpts };
use ray_tracing::scene::Scene;
use ray_tracing::rng::Rng;

#[allow(dead_code)]
fn random_spheres(rng: &mut Rng) -> (Scene, CameraOpts) {
    let mut scene = Scene::new();
    let ground = scene.lambertian(Colour::new(0.5, 0.5, 0.5));
    scene.add_sphere(Point3::new(0.0, -1000.0, 0.0), 1000.0, ground);
    let glass = scene.dielectric(1.5);
    for a in -11..11 {
        for b in -11..11 {
            let center = Point3::new(a as f32 + 0.9 * rng.next_f32(), 0.2, b as f32 + 0.9 * rng.next_f32());
            if (center - &Point3::new(4.0, 0.2, 0.0)).length() > 0.9 {
                let choose = rng.next_f32();
                let material = if choose < 0.8 {
                    let albedo = Colour::random(rng) * &Colour::random(rng);
                    scene.lambertian(albedo)
                } else if choose < 0.95 {
                    let albedo = Vec3::random_in_range(rng, 0.5, 1.0);
                    let fuzz = rng.range(0.0, 0.5);
                    scene.metal(albedo, fuzz)
                } else {
                    glass
                };
                scene.add_sphere(center, 0.2, material);
            }
        }
    }

    let brown = scene.lambertian(Colour::new(0.4, 0.2, 0.1));
    let mirror = scene.metal(Colour::new(0.7, 0.6, 0.5), 0.0);
    scene.add_sphere(Point3::new(-4.0, 1.0, 0.0), 1.0, brown);
    scene.add_sphere(Point3::new(0.0, 1.0, 0.0), 1.0, glass);
    scene.add_sphere(Point3::new(4.0, 1.0, 0.0), 1.0, mirror);

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
    (scene, camera_opts)
}

#[allow(dead_code)]
fn exercise(_rng: &mut Rng) -> (Scene, CameraOpts) {
    let mut scene = Scene::new();
    let material_ground = scene.lambertian(Colour::new(0.8, 0.8, 0.0));
    let material_center = scene.lambertian(Colour::new(0.1, 0.2, 0.5));
    let material_left = scene.dielectric(1.5);
    let material_right = scene.metal(Colour::new(0.8, 0.6, 0.2), 0.1);
    scene.add_sphere(Point3::new(0.0, 0.0, -1.0), 0.5, material_center);
    scene.add_sphere(Point3::new(-1.0, 0.0, -1.0), 0.5, material_left);
    scene.add_sphere(Point3::new(-1.0, 0.0, -1.0), -0.4, material_left);
    scene.add_sphere(Point3::new(1.0, 0.0, -1.0), 0.5, material_right);
    scene.add_sphere(Point3::new(0.0, -100.5, -1.0), 100.0, material_ground);

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
    (scene, camera_opts)
}

#[allow(dead_code)]
fn two_spheres(_rng: &mut Rng) -> (Scene, CameraOpts) {
    let mut scene = Scene::new();
    let checker = scene.add_texture(Texture::checker(0.8, Colour::new(0.2, 0.3, 0.1), Colour::new(0.9, 0.9, 0.9)));
    let material = scene.lambertian_textured(checker);
    scene.add_sphere(Point3::new(0.0, -10.0, 0.0), 10.0, material);
    scene.add_sphere(Point3::new(0.0, 10.0, 0.0), 10.0, material);

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
    (scene, camera_opts)
}

#[allow(dead_code)]
fn earth(_rng: &mut Rng) -> (Scene, CameraOpts) {
    let mut scene = Scene::new();
    let earth_image = image::open("earthmap.jpg").unwrap().to_rgb8();
    let earth_texture = scene.image_texture(earth_image);
    let earth_material = scene.lambertian_textured(earth_texture);
    let material_right = scene.metal(Colour::new(0.8, 0.6, 0.2), 0.0);
    scene.add_sphere(Point3::new(0.0, 0.0, 0.0), 2.0, earth_material);
    scene.add_sphere(Point3::new(5.0, 2.0, -1.0), 3.0, material_right);

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
    (scene, camera_opts)
}

// Samples per pixel in each GPU dispatch. Keep each dispatch short (well under
// a second) so macOS doesn't kill it as a hung GPU job.
const SAMPLES_PER_GPU_PASS: u32 = 10;

fn main() {
    // A fixed seed makes the scene (and BVH) the same on every run.
    let mut rng = Rng::new(42);
    let (mut scene, camera_opts) = random_spheres(&mut rng);
    scene.build_bvh();
    let camera = Camera::initialize(camera_opts);
    // GPU by default; pass --cpu to render with the CPU path instead.
    if std::env::args().any(|arg| arg == "--cpu") {
        camera.render(&scene);
    } else {
        camera.render_gpu(&scene, SAMPLES_PER_GPU_PASS);
    }
}
