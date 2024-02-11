use std::ops::{ Div, Add, Mul, Sub, Index };
use std::ops::Range;
use std::rc::Rc;
use ordered_float::OrderedFloat;
use rand::Rng;
use crate::material::Material;
use std::f64::consts::PI;

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
    pub material: Rc<dyn Material>,
    pub u : f64,
    pub v : f64
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

    fn get_sphere_uv(p: &Point3) -> (f64, f64) {
        // p: a given point on the sphere of radius one, centered at the origin.
        // u: returned value [0,1] of angle around the Y axis from X=-1.
        // v: returned value [0,1] of angle from Y=-1 to Y=+1.
        //     <1 0 0> yields <0.50 0.50>       <-1  0  0> yields <0.00 0.50>
        //     <0 1 0> yields <0.50 1.00>       < 0 -1  0> yields <0.50 0.00>
        //     <0 0 1> yields <0.25 0.50>       < 0  0 -1> yields <0.75 0.50>
        let theta = (-p.y()).acos();
        let phi = (-p.z()).atan2(p.x()) + PI;
        (phi / (2.0*PI), theta/PI)
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
        let (u, v) = Sphere::get_sphere_uv(&unit_normal_to_sphere);
        Some(HitRecord {
            t : root,
            p : position,
            normal : normal,
            front_face: front_face,
            material: self.material.clone(),
            u,
            v
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
