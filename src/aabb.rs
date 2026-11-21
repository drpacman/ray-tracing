use crate::utils::{ Point3, Ray };

// Axis-aligned bounding box.
// See https://raytracing.github.io/books/RayTracingTheNextWeek.html#boundingvolumehierarchies
#[derive(Copy, Clone)]
#[repr(C)]
pub struct Aabb {
    pub min: Point3,
    pub max: Point3
}

impl Aabb {
    // Treat the two points a and b as extrema for the bounding box, so we don't
    // require a particular minimum/maximum coordinate order.
    pub fn new(a: Point3, b: Point3) -> Self {
        Self { min: a.min(&b), max: a.max(&b) }
    }

    // The smallest box that encloses both boxes.
    pub fn surrounding(box0: &Aabb, box1: &Aabb) -> Self {
        Self { min: box0.min.min(&box1.min), max: box0.max.max(&box1.max) }
    }

    // Index (0 = x, 1 = y, 2 = z) of the box's longest side.
    pub fn longest_axis(&self) -> usize {
        let size = self.max - &self.min;
        if size.x() > size.y() {
            if size.x() > size.z() { 0 } else { 2 }
        } else {
            if size.y() > size.z() { 1 } else { 2 }
        }
    }

    // The "slab" method: for each axis, compute the t interval over which the
    // ray is between the two planes bounding the box on that axis. The ray hits
    // the box only if those intervals all overlap (with each other and the
    // incoming range).
    pub fn hit(&self, ray: &Ray, t_min: f32, t_max: f32) -> bool {
        let mut t_min = t_min;
        let mut t_max = t_max;
        for axis in 0..3 {
            let inv_d = 1.0 / ray.dir[axis];
            let t0 = (self.min[axis] - ray.orig[axis]) * inv_d;
            let t1 = (self.max[axis] - ray.orig[axis]) * inv_d;
            // A negative direction component flips which plane is entered first.
            let (t0, t1) = if t0 < t1 { (t0, t1) } else { (t1, t0) };
            if t0 > t_min { t_min = t0; }
            if t1 < t_max { t_max = t1; }
            if t_max <= t_min {
                return false;
            }
        }
        true
    }
}
