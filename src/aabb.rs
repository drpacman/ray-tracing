use crate::utils::{ Point3, Ray };
use ordered_float::OrderedFloat;
use std::ops::Range;

// Axis-aligned bounding box.
// See https://raytracing.github.io/books/RayTracingTheNextWeek.html#boundingvolumehierarchies
#[derive(Copy, Clone)]
pub struct Aabb {
    pub min: Point3,
    pub max: Point3
}

impl Aabb {
    // Treat the two points a and b as extrema for the bounding box, so we don't
    // require a particular minimum/maximum coordinate order.
    pub fn new(a: Point3, b: Point3) -> Self {
        Self {
            min: Point3::new(a.x().min(b.x()), a.y().min(b.y()), a.z().min(b.z())),
            max: Point3::new(a.x().max(b.x()), a.y().max(b.y()), a.z().max(b.z()))
        }
    }

    // The smallest box that encloses both boxes.
    pub fn surrounding(box0: &Aabb, box1: &Aabb) -> Self {
        Self {
            min: Point3::new(box0.min.x().min(box1.min.x()), box0.min.y().min(box1.min.y()), box0.min.z().min(box1.min.z())),
            max: Point3::new(box0.max.x().max(box1.max.x()), box0.max.y().max(box1.max.y()), box0.max.z().max(box1.max.z()))
        }
    }

    // The "slab" method: for each axis, compute the t interval over which the
    // ray is between the two planes bounding the box on that axis. The ray hits
    // the box only if those intervals all overlap (with each other and the
    // incoming range).
    pub fn hit(&self, ray: &Ray, range: &Range<OrderedFloat<f64>>) -> bool {
        let mut t_min = range.start.into_inner();
        let mut t_max = range.end.into_inner();
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
