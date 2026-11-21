use crate::aabb::Aabb;
use crate::utils::{ HitRecord, Hittable, Ray };
use ordered_float::OrderedFloat;
use rand::Rng;
use std::ops::Range;
use std::sync::Arc;

// A node in a bounding volume hierarchy: a binary tree where each node's box
// encloses everything below it, so a ray that misses the box can skip the
// whole subtree.
// See https://raytracing.github.io/books/RayTracingTheNextWeek.html#boundingvolumehierarchies
pub struct BvhNode {
    left: Arc<dyn Hittable>,
    right: Arc<dyn Hittable>,
    bbox: Aabb
}

impl BvhNode {
    pub fn from_list<T: Hittable + 'static>(objects: Vec<T>) -> Self {
        let mut objects: Vec<Arc<dyn Hittable>> = objects
            .into_iter()
            .map(|o| Arc::new(o) as Arc<dyn Hittable>)
            .collect();
        Self::new(&mut objects)
    }

    pub fn new(objects: &mut [Arc<dyn Hittable>]) -> Self {
        assert!(!objects.is_empty(), "cannot build a BVH from an empty list");

        // The simple approach from the book: split along a randomly chosen axis.
        let axis = rand::thread_rng().gen_range(0..3);

        let (left, right): (Arc<dyn Hittable>, Arc<dyn Hittable>) = match objects.len() {
            // A single object goes in both children, so every node has two
            // children and we don't need to check for null.
            1 => (objects[0].clone(), objects[0].clone()),
            2 => (objects[0].clone(), objects[1].clone()),
            n => {
                objects.sort_by(|a, b| {
                    a.bounding_box().min[axis].total_cmp(&b.bounding_box().min[axis])
                });
                let (lower, upper) = objects.split_at_mut(n / 2);
                (Arc::new(BvhNode::new(lower)), Arc::new(BvhNode::new(upper)))
            }
        };

        let bbox = Aabb::surrounding(&left.bounding_box(), &right.bounding_box());
        Self { left, right, bbox }
    }
}

impl Hittable for BvhNode {
    fn hit(&self, ray: &Ray, range: &Range<OrderedFloat<f64>>) -> Option<HitRecord> {
        if !self.bbox.hit(ray, range) {
            return None;
        }

        let hit_left = self.left.hit(ray, range);
        // Only accept a right-hand hit if it is closer than the left-hand one.
        let right_range = match &hit_left {
            Some(rec) => range.start..OrderedFloat(rec.t),
            None => range.clone()
        };
        let hit_right = self.right.hit(ray, &right_range);
        hit_right.or(hit_left)
    }

    fn bounding_box(&self) -> Aabb {
        self.bbox
    }
}
