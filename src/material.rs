use crate::utils::{ Vec3, Ray, HitRecord, Colour };
use rand::Rng;
use std::rc::Rc;
use image::{GenericImageView, Pixel};
pub trait Material {
    fn scatter(&self, ray: &Ray, hit_record: &HitRecord) -> Option<(Colour, Ray)>;
}

/*
Albedo - represents the _fraction_ of light that is reflected by a body or surface. 
It is commonly used in astronomy to describe the reflective properties of 
planets, satellites, and asteroids.
 */
pub struct Lambertian {
    pub albedo: Rc<dyn Texture>
}

impl Lambertian {
    pub fn from_colour(albedo: Colour) -> Self {
        Self::new(Rc::new(SolidColour::new(albedo)))
    }

    pub fn new(albedo: Rc<dyn Texture>) -> Self {
        Self { albedo }
    }
}

impl Material for Lambertian {    
    fn scatter(&self, _ray: &Ray, hit_record: &HitRecord) -> Option<(Colour, Ray)> {
        let mut scatter_direction = hit_record.normal + &Vec3::random_unit_vector();
        if scatter_direction.near_zero() {
            scatter_direction = hit_record.normal.clone();
        }
        let scattered = Ray::new(hit_record.p.clone(), scatter_direction);
        let attenuation = self.albedo.value(hit_record.u, hit_record.v, &hit_record.p);
        Some((attenuation, scattered))
    }
}

pub struct Metal {
    pub albedo: Colour,
    pub fuzz: f64
}

impl Metal {
    pub fn new(albedo: Colour, fuzz: f64) -> Self {
        Self { albedo, fuzz: if fuzz < 1.0 { fuzz } else { 1.0 } }
    }
}

impl Material for Metal {    
    fn scatter(&self, ray: &Ray, hit_record: &HitRecord) -> Option<(Colour, Ray)> {
        // For mirrored reflection
        // See https://raytracing.github.io/books/RayTracingInOneWeekend.html#metal/mirroredlightreflection    
        let mirror_reflected_ray = ray.dir.unit_vector().reflect(&hit_record.normal);
        // The fuzz parameter is used to simulate the roughness of the metal. A fuzz of zero is a perfect mirror.
        // See https://raytracing.github.io/books/RayTracingInOneWeekend.html#metal/fuzzyreflection
        let fuzzed_ray = mirror_reflected_ray + &(Vec3::random_in_unit_sphere() * self.fuzz);
        let scattered = Ray::new(hit_record.p.clone(), fuzzed_ray);
        Some((self.albedo.clone(), scattered))
    }
}

// See https://raytracing.github.io/books/RayTracingInOneWeekend.html#dielectrics
pub struct Dielectric {
    pub index_of_refraction: f64
}

impl Dielectric {
    pub fn new(index_of_refraction: f64) -> Self {
        Self { index_of_refraction }
    }

    // Uses Schlick's approximation for reflectance.
    // see https://en.wikipedia.org/wiki/Schlick%27s_approximation
    fn reflectance(cos_theta : f64, refraction_ratio : f64) -> f64 {
        let r0 = ((1.0 - refraction_ratio) / (1.0 + refraction_ratio)).powi(2);
        r0 + (1.0 - r0) * (1.0 - cos_theta).powi(5)
    }
}

impl Material for Dielectric {
    fn scatter(&self, ray: &Ray, hit_record: &HitRecord) -> Option<(Colour, Ray)> {
        // Example refractive indices (typically air = 1.0, glass = 1.3–1.7, diamond = 2.4)
        // assuming front face we are going from air -> material and backface we are going from material to air
        // refraction_ratio is the ratio of the refractive indices of the two materials.
        let refraction_ratio = if hit_record.front_face { 1.0 / self.index_of_refraction } else { self.index_of_refraction };
        let unit_direction = ray.dir.unit_vector();
        let cos_theta = f64::min(1.0, unit_direction.dot(&hit_record.normal) * -1.0);
        let sin_theta = (1.0 - cos_theta * cos_theta).sqrt();
        // when the ray is in the material with the higher refractive index, depending on the value of sin theta
        // there can be no real solution to Snell’s law. In those cases there is no refraction possible.
        let cannot_refract = refraction_ratio * sin_theta > 1.0;
        let direction = if cannot_refract || Self::reflectance(cos_theta, refraction_ratio) > rand::thread_rng().gen_range(0.0..1.0){
            unit_direction.reflect(&hit_record.normal)
        } else {
            unit_direction.refract(&hit_record.normal, refraction_ratio)
        };
        Some((Colour::new(1.0, 1.0, 1.0), Ray::new(hit_record.p.clone(), direction)))
    }
}

pub trait Texture {
    fn value(&self, u: f64, v: f64, p: &Vec3) -> Colour;
}

pub struct SolidColour {
    pub colour: Colour
}

impl SolidColour {
    pub fn new(colour: Colour) -> Self {
        Self { colour }
    }
}

impl Texture for SolidColour {
    fn value(&self, _u: f64, _v: f64, _p: &Vec3) -> Colour {
        self.colour.clone()
    }
}

pub struct CheckerTexture {
    pub inverse_scale : f64,
    pub odd: Box<dyn Texture>,
    pub even: Box<dyn Texture>
}

impl CheckerTexture {
    pub fn new_with_colours(scale: f64, colour1: Colour, colour2: Colour) -> Self {
        Self { 
            inverse_scale: 1.0/scale, 
            odd : Box::new(SolidColour::new(colour1)), 
            even : Box::new(SolidColour::new(colour2))
        }
    }

    pub fn new_with_textures(scale: f64, texture1: Box<dyn Texture>, texture2: Box<dyn Texture>) -> Self {
        Self { 
            inverse_scale: 1.0/scale, 
            odd : texture1, 
            even : texture2
        }
    }
}

impl Texture for CheckerTexture {
    fn value(&self, u: f64, v: f64, p: &Vec3) -> Colour {
        let x_int : i64 = (self.inverse_scale * p.x()).floor() as i64;
        let y_int : i64 = (self.inverse_scale * p.y()).floor() as i64;
        let z_int : i64 = (self.inverse_scale * p.z()).floor() as i64;
        let is_even = (x_int + y_int + z_int) % 2 == 0;
        if is_even {
            self.even.value(u, v, p)
        } else {
            self.odd.value(u, v, p)
        } 
    }
}

pub struct ImageTexture {
    pub image: image::DynamicImage
}

impl ImageTexture {
    pub fn new(image: image::DynamicImage) -> Self {
        Self { image }
    }
}

impl Texture for ImageTexture {
    fn value(&self, u: f64, v: f64, _p: &Vec3) -> Colour {
        // Clamp input texture coordinates to [0,1] x [1,0]
        let u = f64::min(f64::max(u, 0.0), 1.0);
        let v = 1.0 - f64::min(f64::max(v, 0.0), 1.0);
        let (width, height) = self.image.dimensions();
        let x = (u * width as f64) as u32;
        let y = (v * height as f64) as u32;
        let pixel = self.image.get_pixel(x, y);
        let channels = pixel.channels();
        Colour::new(channels[0] as f64 / 255.0, channels[1] as f64 / 255.0, channels[2] as f64 / 255.0)
    }
}