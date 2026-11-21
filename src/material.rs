use crate::utils::{ Vec3, Ray, HitRecord, Colour };
use crate::scene::Scene;
use crate::rng::Rng;

// Materials are a plain enum rather than trait objects, so they can live in a
// flat array (and later a GPU buffer) and be selected with a match/switch.
#[derive(Copy, Clone)]
pub enum Material {
    /*
    Albedo - represents the _fraction_ of light that is reflected by a body or surface.
    It is commonly used in astronomy to describe the reflective properties of
    planets, satellites, and asteroids.
    */
    Lambertian { albedo: u32 }, // index into Scene::textures
    Metal { albedo: Colour, fuzz: f32 },
    // See https://raytracing.github.io/books/RayTracingInOneWeekend.html#dielectrics
    Dielectric { index_of_refraction: f32 }
}

// Uses Schlick's approximation for reflectance.
// see https://en.wikipedia.org/wiki/Schlick%27s_approximation
fn reflectance(cos_theta : f32, refraction_ratio : f32) -> f32 {
    let r0 = ((1.0 - refraction_ratio) / (1.0 + refraction_ratio)).powi(2);
    r0 + (1.0 - r0) * (1.0 - cos_theta).powi(5)
}

impl Material {
    pub fn scatter(&self, ray: &Ray, hit_record: &HitRecord, scene: &Scene, rng: &mut Rng) -> Option<(Colour, Ray)> {
        match *self {
            Material::Lambertian { albedo } => {
                let mut scatter_direction = hit_record.normal + &Vec3::random_unit_vector(rng);
                if scatter_direction.near_zero() {
                    scatter_direction = hit_record.normal;
                }
                let scattered = Ray::new(hit_record.p, scatter_direction);
                let attenuation = scene.texture_value(albedo, hit_record.u, hit_record.v, &hit_record.p);
                Some((attenuation, scattered))
            }
            Material::Metal { albedo, fuzz } => {
                // For mirrored reflection
                // See https://raytracing.github.io/books/RayTracingInOneWeekend.html#metal/mirroredlightreflection
                let mirror_reflected_ray = ray.dir.unit_vector().reflect(&hit_record.normal);
                // The fuzz parameter is used to simulate the roughness of the metal. A fuzz of zero is a perfect mirror.
                // See https://raytracing.github.io/books/RayTracingInOneWeekend.html#metal/fuzzyreflection
                let fuzzed_ray = mirror_reflected_ray + &(Vec3::random_in_unit_sphere(rng) * fuzz);
                Some((albedo, Ray::new(hit_record.p, fuzzed_ray)))
            }
            Material::Dielectric { index_of_refraction } => {
                // Example refractive indices (typically air = 1.0, glass = 1.3–1.7, diamond = 2.4)
                // assuming front face we are going from air -> material and backface we are going from material to air
                // refraction_ratio is the ratio of the refractive indices of the two materials.
                let refraction_ratio = if hit_record.front_face { 1.0 / index_of_refraction } else { index_of_refraction };
                let unit_direction = ray.dir.unit_vector();
                let cos_theta = f32::min(1.0, unit_direction.dot(&hit_record.normal) * -1.0);
                let sin_theta = (1.0 - cos_theta * cos_theta).sqrt();
                // when the ray is in the material with the higher refractive index, depending on the value of sin theta
                // there can be no real solution to Snell’s law. In those cases there is no refraction possible.
                let cannot_refract = refraction_ratio * sin_theta > 1.0;
                let direction = if cannot_refract || reflectance(cos_theta, refraction_ratio) > rng.next_f32() {
                    unit_direction.reflect(&hit_record.normal)
                } else {
                    unit_direction.refract(&hit_record.normal, refraction_ratio)
                };
                Some((Colour::new(1.0, 1.0, 1.0), Ray::new(hit_record.p, direction)))
            }
        }
    }
}

#[derive(Copy, Clone)]
pub enum Texture {
    Solid(Colour),
    Checker { inverse_scale: f32, odd: Colour, even: Colour },
    Image { image: u32 } // index into Scene::images
}

impl Texture {
    pub fn checker(scale: f32, odd: Colour, even: Colour) -> Self {
        Texture::Checker { inverse_scale: 1.0 / scale, odd, even }
    }

    // Image lookups need the scene's images, so Scene::texture_value calls this.
    pub fn value(&self, u: f32, v: f32, p: &Vec3, images: &[image::RgbImage]) -> Colour {
        match *self {
            Texture::Solid(colour) => colour,
            Texture::Checker { inverse_scale, odd, even } => {
                let x_int = (inverse_scale * p.x()).floor() as i64;
                let y_int = (inverse_scale * p.y()).floor() as i64;
                let z_int = (inverse_scale * p.z()).floor() as i64;
                if (x_int + y_int + z_int) % 2 == 0 { even } else { odd }
            }
            Texture::Image { image } => {
                let image = &images[image as usize];
                // Clamp input texture coordinates to [0,1] x [1,0]
                let u = u.clamp(0.0, 1.0);
                let v = 1.0 - v.clamp(0.0, 1.0);
                let (width, height) = image.dimensions();
                let x = ((u * width as f32) as u32).min(width - 1);
                let y = ((v * height as f32) as u32).min(height - 1);
                let pixel = image.get_pixel(x, y);
                Colour::new(pixel[0] as f32 / 255.0, pixel[1] as f32 / 255.0, pixel[2] as f32 / 255.0)
            }
        }
    }
}
