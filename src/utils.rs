use std::ops::{ Div, Add, Mul, Sub, Index };
use std::ops::Range;
use std::rc::Rc;
use ordered_float::OrderedFloat;
use rand::Rng;

#[derive(Copy, Clone)]
pub struct Vec3 {
    e: [f64; 3]
}

impl Vec3 {
    pub fn new(e0: f64, e1: f64, e2: f64) -> Self {
        Self { e: [e0, e1, e2] }
    }

    pub fn x(&self) -> f64 { self.e[0] }
    pub fn y(&self) -> f64 { self.e[1] }
    pub fn z(&self) -> f64 { self.e[2] }

    pub fn r(&self) -> f64 { self.e[0] }
    pub fn g(&self) -> f64 { self.e[1] }
    pub fn b(&self) -> f64 { self.e[2] }

    pub fn length(&self) -> f64 {
        self.length_squared().sqrt()
    }

    pub fn length_squared(&self) -> f64 {
        self.e[0] * self.e[0] + self.e[1] * self.e[1] + self.e[2] * self.e[2]
    }

    pub fn near_zero(&self) -> bool {
        let s = 1e-8;
        self.e[0].abs() < s && self.e[1].abs() < s && self.e[2].abs() < s
    } 
    
    pub fn dot(&self, rhs: &Self) -> f64 {
        self.e[0] * rhs.e[0] + self.e[1] * rhs.e[1] + self.e[2] * rhs.e[2]
    }

    pub fn cross(&self, rhs: &Self) -> Self {
        Self::new(
            self.e[1] * rhs.e[2] - self.e[2] * rhs.e[1],
            self.e[2] * rhs.e[0] - self.e[0] * rhs.e[2],
            self.e[0] * rhs.e[1] - self.e[1] * rhs.e[0]
        )
    }

    pub fn unit_vector(&self) -> Self {
        self / self.length()
    }
        
    pub fn reflect(&self, normal: &Vec3) -> Self {
        *self - &(*normal * self.dot(normal) * 2.0)
    }

    // See https://raytracing.github.io/books/RayTracingInOneWeekend.html#dielectrics/snell'slaw
    pub fn refract(&self, normal: &Vec3, etai_over_etat: f64) -> Self {
        let cos_theta = f64::min(1.0, (*self * -1.0).dot(normal));
        let ray_out_perpendicular = (*self + &((*normal * cos_theta))) * etai_over_etat;
        let ray_out_parallel = *normal * -1.0*f64::abs(1.0 - ray_out_perpendicular.length_squared()).sqrt();
        ray_out_perpendicular + &ray_out_parallel
    }

    pub fn random() -> Vec3 {
        Vec3::random_in_range(0.0, 1.0)
    }
    
    pub fn random_in_range(min: f64, max:f64) -> Vec3 {
        let mut rng = rand::thread_rng();
        Vec3::new(rng.gen_range(min..max), rng.gen_range(min..max), rng.gen_range(min..max))
    }
    
    pub fn random_in_unit_sphere() -> Vec3 {
        loop {
            let p = Vec3::random_in_range(-1.0, 1.0);
            if p.length_squared() >= 1.0 {
                continue;
            }
            return p;
        }
    }

    pub fn random_unit_vector() -> Vec3 {
        Self::random_in_unit_sphere().unit_vector()
    }

    pub fn random_on_hemisphere(normal: &Vec3) -> Vec3 {
        let on_unit_sphere = Self::random_unit_vector();
        if on_unit_sphere.dot(normal) > 0.0 {
            on_unit_sphere
        } else {
            on_unit_sphere * -1.0
        }
    }

    pub fn random_in_unit_disk() -> Self {
        loop {
            let p = Vec3::new(rand::thread_rng().gen_range(-1.0..1.0), rand::thread_rng().gen_range(-1.0..1.0), 0.0);
            if p.length_squared() < 1.0 {
                return p
            }
        }
    }   
}

impl Add<&Vec3> for Vec3 {
    type Output = Self;

    fn add(self, rhs: &Self) -> Self::Output {
        Self::new(self.e[0] + rhs.e[0], self.e[1] + rhs.e[1], self.e[2] + rhs.e[2])
    }
}

impl <'a> Div<&'a Vec3> for &'a Vec3 {
    type Output = Vec3;

    fn div(self, rhs: Self) -> Self::Output {
        Vec3::new(self.e[0] / rhs.e[0], self.e[1] / rhs.e[1], self.e[2] / rhs.e[2])
    }
}

impl Div<f64> for &Vec3 {
    type Output = Vec3;

    fn div(self, rhs: f64) -> Self::Output {
        Vec3::new(self.e[0] / rhs, self.e[1] / rhs, self.e[2] / rhs)
    }
}

impl Sub<&Vec3> for Vec3 {
    type Output = Self;

    fn sub(self, rhs: &Self) -> Self::Output {
        Self::new(self.e[0] - rhs.e[0], self.e[1] - rhs.e[1], self.e[2] - rhs.e[2])
    }
}

impl Mul<&Vec3> for Vec3 {
    type Output = Self;

    fn mul(self, rhs: &Self) -> Self::Output {
        Self::new(self.e[0] * rhs.e[0], self.e[1] * rhs.e[1], self.e[2] * rhs.e[2])
    }
}

impl Mul<f64> for Vec3 {
    type Output = Self;

    fn mul(self, rhs: f64) -> Self::Output {
        Self::new(self.e[0] * rhs, self.e[1] * rhs, self.e[2] * rhs)
    }
}

impl Index<usize> for Vec3 {
    type Output = f64;

    fn index(&self, idx: usize) -> &Self::Output {
        &self.e[idx]
    }
}

pub type Point3 = Vec3;
pub type Colour = Vec3;

fn linear_to_gamma(x: f64) -> f64 {
    x.sqrt()
}
pub fn write_colour(pixel_color: &Vec3, samples_per_pixel: i32) {
    let scale = 1.0 / samples_per_pixel as f64;
    let r = linear_to_gamma(scale * pixel_color.r());
    let g = linear_to_gamma(scale * pixel_color.g());
    let b = linear_to_gamma(scale * pixel_color.b());
    let intensity = OrderedFloat(0.0)..OrderedFloat(0.999);
    
    println!("{} {} {}", (256.0 * clamp(r, &intensity)) as i32, (256.0 * clamp(g, &intensity)) as i32, (256.0 * clamp(b, &intensity)) as i32);
}

pub fn clamp(t: f64, range : &Range<OrderedFloat<f64>>) -> f64 {
    if t < range.start.into_inner() {
        range.start.into_inner()
    } else if t > range.end.into_inner() {
        range.end.into_inner()
    } else {
        t
    }
}

pub struct Ray {
    pub orig: Point3,
    pub dir: Vec3
}

impl Ray {
    pub fn new(origin: Point3, direction: Vec3) -> Self {
        Self { orig: origin, dir: direction }
    }

    pub fn origin(&self) -> &Point3 {
        &self.orig
    }

    pub fn direction(&self) -> &Vec3 {
        &self.dir
    }

    pub fn at(&self, t: f64) -> Point3 {
        self.orig + &( (self.dir * &Vec3::new(t, t, t)))
    }
}

pub struct HitRecord {
    pub p: Point3,
    pub normal: Vec3,
    pub t: f64,
    pub front_face: bool,
    pub material: Rc<dyn Material>
}

pub trait Material {
    fn scatter(&self, ray: &Ray, hit_record: &HitRecord) -> Option<(Colour, Ray)>;
}

/*
Albedo - represents the _fraction_ of light that is reflected by a body or surface. 
It is commonly used in astronomy to describe the reflective properties of 
planets, satellites, and asteroids.
 */
pub struct Lambertian {
    pub albedo: Colour
}

impl Lambertian {
    pub fn new(albedo: Colour) -> Self {
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
        Some((self.albedo.clone(), scattered))
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

pub trait Hittable {
    fn hit(&self, ray: &Ray, range : &Range<OrderedFloat<f64>>) -> Option<HitRecord>;
}

pub struct Sphere {
    center: Point3,
    radius: f64,
    material: Rc<dyn Material>
}

impl Sphere {
    pub fn new(center: Point3, radius: f64, material: Rc<dyn Material>) -> Self {
        Self { center, radius, material }
    }
}

impl Hittable for Sphere {
    fn  hit(&self, ray: &Ray, range : &Range<OrderedFloat<f64>>) -> Option<HitRecord> {        
        let oc = ray.origin().clone() - &self.center;
        let a = ray.direction().dot(ray.direction());
        let half_b = oc.dot(ray.direction());
        let c = oc.dot(&oc) - &self.radius * &self.radius;
        let discriminant = half_b*half_b - a*c;
        if discriminant < 0.0 {
            return None;
        }

        let sqrtd = discriminant.sqrt();
        let mut root = (-half_b - sqrtd)/a;
        if !range.contains(&OrderedFloat(root)) {
            root = (-half_b + sqrtd)/a;
            if !range.contains(&OrderedFloat(root)) {
                return None;
            }
        }
        let position = ray.at(root);
        let unit_normal_to_sphere = &(position - &self.center) / self.radius;
        // if angle to the normal is between 90 and and 270, we are on the front face.
        // cos is negative for angles between 90 and 270.
        let front_face = ray.direction().dot(&unit_normal_to_sphere) < 0.0;
        let normal = if front_face { unit_normal_to_sphere.clone() } else { unit_normal_to_sphere.clone() * -1.0 };

        Some(HitRecord {
            t : root,
            p : position,
            normal : normal,
            front_face: front_face,
            material: self.material.clone()
        })
    }
}

impl <T> Hittable for Vec<T> where T: Hittable {
    fn hit(&self, ray: &Ray, range: &Range<OrderedFloat<f64>>) -> Option<HitRecord> {
        let mut shortest_range = range.clone();
        let mut hit = None;
        for object in self.into_iter() {
            if let Some(rec) = object.hit(ray, &shortest_range) {
                shortest_range = shortest_range.start..OrderedFloat(rec.t);
                hit = Some(rec);
            }
        }
        return hit;
    }
}
