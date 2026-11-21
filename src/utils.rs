use std::ops::{ Div, Add, Mul, Sub, Index };
use std::io::{self, Write};
use std::f32::consts::PI;
use crate::aabb::Aabb;
use crate::rng::Rng;

// #[repr(C)] gives a fixed, C-compatible memory layout, ready for copying into GPU buffers.
#[derive(Copy, Clone, Debug)]
#[repr(C)]
pub struct Vec3 {
    e: [f32; 3]
}

impl Vec3 {
    pub fn new(e0: f32, e1: f32, e2: f32) -> Self {
        Self { e: [e0, e1, e2] }
    }

    pub fn x(&self) -> f32 { self.e[0] }
    pub fn y(&self) -> f32 { self.e[1] }
    pub fn z(&self) -> f32 { self.e[2] }

    pub fn r(&self) -> f32 { self.e[0] }
    pub fn g(&self) -> f32 { self.e[1] }
    pub fn b(&self) -> f32 { self.e[2] }

    pub fn length(&self) -> f32 {
        self.length_squared().sqrt()
    }

    pub fn length_squared(&self) -> f32 {
        self.e[0] * self.e[0] + self.e[1] * self.e[1] + self.e[2] * self.e[2]
    }

    pub fn near_zero(&self) -> bool {
        let s = 1e-8;
        self.e[0].abs() < s && self.e[1].abs() < s && self.e[2].abs() < s
    }

    pub fn dot(&self, rhs: &Self) -> f32 {
        self.e[0] * rhs.e[0] + self.e[1] * rhs.e[1] + self.e[2] * rhs.e[2]
    }

    pub fn cross(&self, rhs: &Self) -> Self {
        Self::new(
            self.e[1] * rhs.e[2] - self.e[2] * rhs.e[1],
            self.e[2] * rhs.e[0] - self.e[0] * rhs.e[2],
            self.e[0] * rhs.e[1] - self.e[1] * rhs.e[0]
        )
    }

    pub fn min(&self, rhs: &Self) -> Self {
        Self::new(self.e[0].min(rhs.e[0]), self.e[1].min(rhs.e[1]), self.e[2].min(rhs.e[2]))
    }

    pub fn max(&self, rhs: &Self) -> Self {
        Self::new(self.e[0].max(rhs.e[0]), self.e[1].max(rhs.e[1]), self.e[2].max(rhs.e[2]))
    }

    pub fn unit_vector(&self) -> Self {
        self / self.length()
    }

    pub fn reflect(&self, normal: &Vec3) -> Self {
        *self - &(*normal * self.dot(normal) * 2.0)
    }

    // See https://raytracing.github.io/books/RayTracingInOneWeekend.html#dielectrics/snell'slaw
    pub fn refract(&self, normal: &Vec3, etai_over_etat: f32) -> Self {
        let cos_theta = f32::min(1.0, (*self * -1.0).dot(normal));
        let ray_out_perpendicular = (*self + &((*normal * cos_theta))) * etai_over_etat;
        let ray_out_parallel = *normal * -1.0*f32::abs(1.0 - ray_out_perpendicular.length_squared()).sqrt();
        ray_out_perpendicular + &ray_out_parallel
    }

    pub fn random(rng: &mut Rng) -> Vec3 {
        Vec3::random_in_range(rng, 0.0, 1.0)
    }

    pub fn random_in_range(rng: &mut Rng, min: f32, max: f32) -> Vec3 {
        Vec3::new(rng.range(min, max), rng.range(min, max), rng.range(min, max))
    }

    pub fn random_in_unit_sphere(rng: &mut Rng) -> Vec3 {
        loop {
            let p = Vec3::random_in_range(rng, -1.0, 1.0);
            if p.length_squared() < 1.0 {
                return p;
            }
        }
    }

    pub fn random_unit_vector(rng: &mut Rng) -> Vec3 {
        Self::random_in_unit_sphere(rng).unit_vector()
    }

    pub fn random_in_unit_disk(rng: &mut Rng) -> Self {
        loop {
            let p = Vec3::new(rng.range(-1.0, 1.0), rng.range(-1.0, 1.0), 0.0);
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

impl Div<f32> for &Vec3 {
    type Output = Vec3;

    fn div(self, rhs: f32) -> Self::Output {
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

impl Mul<f32> for Vec3 {
    type Output = Self;

    fn mul(self, rhs: f32) -> Self::Output {
        Self::new(self.e[0] * rhs, self.e[1] * rhs, self.e[2] * rhs)
    }
}

impl Index<usize> for Vec3 {
    type Output = f32;

    fn index(&self, idx: usize) -> &Self::Output {
        &self.e[idx]
    }
}

pub type Point3 = Vec3;
pub type Colour = Vec3;

fn linear_to_gamma(x: f32) -> f32 {
    x.sqrt()
}

pub fn write_colour(out: &mut impl Write, pixel_color: &Vec3, samples_per_pixel: i32) -> io::Result<()> {
    let scale = 1.0 / samples_per_pixel as f32;
    let to_byte = |c: f32| (256.0 * linear_to_gamma(scale * c).clamp(0.0, 0.999)) as i32;
    writeln!(out, "{} {} {}", to_byte(pixel_color.r()), to_byte(pixel_color.g()), to_byte(pixel_color.b()))
}

#[derive(Copy, Clone)]
pub struct Ray {
    pub orig: Point3,
    pub dir: Vec3
}

impl Ray {
    pub fn new(origin: Point3, direction: Vec3) -> Self {
        Self { orig: origin, dir: direction }
    }

    pub fn at(&self, t: f32) -> Point3 {
        self.orig + &(self.dir * t)
    }
}

pub struct HitRecord {
    pub p: Point3,
    pub normal: Vec3,
    pub t: f32,
    pub front_face: bool,
    pub material: u32, // index into Scene::materials
    pub u : f32,
    pub v : f32
}

#[derive(Copy, Clone)]
#[repr(C)]
pub struct Sphere {
    pub center: Point3,
    pub radius: f32,
    pub material: u32 // index into Scene::materials
}

impl Sphere {
    pub fn new(center: Point3, radius: f32, material: u32) -> Self {
        Self { center, radius, material }
    }

    fn get_sphere_uv(p: &Point3) -> (f32, f32) {
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

    pub fn bounding_box(&self) -> Aabb {
        // abs() because a negative radius is used to model hollow glass spheres.
        let r = self.radius.abs();
        let rvec = Vec3::new(r, r, r);
        Aabb::new(self.center - &rvec, self.center + &rvec)
    }

    // Returns only the distance t to the nearest hit in [t_min, t_max).
    // Building the full HitRecord (normal, uv, ...) is deferred until we know
    // this is the closest object, see hit_record().
    pub fn intersect(&self, ray: &Ray, t_min: f32, t_max: f32) -> Option<f32> {
        let oc = ray.orig - &self.center;
        let a = ray.dir.length_squared();
        let half_b = oc.dot(&ray.dir);
        // The textbook discriminant half_b² - a·c subtracts two large, nearly
        // equal numbers for big spheres (the radius-1000 ground), which loses
        // too much precision in f32. This equivalent form measures the distance
        // from the centre to the ray's closest point instead.
        // See Ray Tracing Gems, ch. 7 "Precision Improvements for Ray/Sphere Intersection".
        let l = oc - &(ray.dir * (half_b / a));
        let discriminant = a * (self.radius * self.radius - l.length_squared());
        if discriminant < 0.0 {
            return None;
        }

        let sqrtd = discriminant.sqrt();
        let mut root = (-half_b - sqrtd) / a;
        if root < t_min || root >= t_max {
            root = (-half_b + sqrtd) / a;
            if root < t_min || root >= t_max {
                return None;
            }
        }
        Some(root)
    }

    pub fn hit_record(&self, ray: &Ray, t: f32) -> HitRecord {
        let position = ray.at(t);
        let unit_normal_to_sphere = &(position - &self.center) / self.radius;
        // if angle to the normal is between 90 and and 270, we are on the front face.
        // cos is negative for angles between 90 and 270.
        let front_face = ray.dir.dot(&unit_normal_to_sphere) < 0.0;
        let normal = if front_face { unit_normal_to_sphere } else { unit_normal_to_sphere * -1.0 };
        let (u, v) = Sphere::get_sphere_uv(&unit_normal_to_sphere);
        HitRecord { t, p: position, normal, front_face, material: self.material, u, v }
    }
}
