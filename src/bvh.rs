use crate::aabb::Aabb;
use crate::rng::Rng;
use crate::utils::Sphere;

// A bounding volume hierarchy stored as a flat array of nodes rather than a
// tree of pointers, so it can be copied straight into a GPU buffer.
// Node 0 is the root.
//   - Leaf (count > 0): holds spheres[first .. first + count].
//   - Interior (count == 0): its two children are nodes[first] and nodes[first + 1].
// See https://raytracing.github.io/books/RayTracingTheNextWeek.html#boundingvolumehierarchies
#[derive(Copy, Clone)]
#[repr(C)]
pub struct BvhNode {
    pub bbox: Aabb,
    pub first: u32,
    pub count: u32
}

impl BvhNode {
    pub fn is_leaf(&self) -> bool {
        self.count > 0
    }
}

// Builds the BVH, reordering `spheres` so that each leaf's spheres are contiguous.
pub fn build(spheres: &mut [Sphere], rng: &mut Rng) -> Vec<BvhNode> {
    assert!(!spheres.is_empty(), "cannot build a BVH from an empty list");
    let placeholder = BvhNode { bbox: spheres[0].bounding_box(), first: 0, count: 0 };
    let mut nodes = vec![placeholder];
    subdivide(&mut nodes, 0, spheres, 0, rng);
    nodes
}

// Fills in nodes[node_index] for the spheres in `spheres`, which start at
// index `offset` in the full sphere array.
fn subdivide(nodes: &mut Vec<BvhNode>, node_index: usize, spheres: &mut [Sphere], offset: usize, rng: &mut Rng) {
    let bbox = spheres.iter()
        .map(|s| s.bounding_box())
        .reduce(|a, b| Aabb::surrounding(&a, &b))
        .unwrap();

    // The book stops at one or two objects per leaf.
    if spheres.len() <= 2 {
        nodes[node_index] = BvhNode { bbox, first: offset as u32, count: spheres.len() as u32 };
        return;
    }

    // The simple approach from the book: sort along a random axis and split in half.
    let axis = rng.below(3) as usize;
    spheres.sort_by(|a, b| a.bounding_box().min[axis].total_cmp(&b.bounding_box().min[axis]));
    let mid = spheres.len() / 2;

    // Allocate both children next to each other.
    let left = nodes.len();
    nodes.push(BvhNode { bbox, first: 0, count: 0 });
    nodes.push(BvhNode { bbox, first: 0, count: 0 });
    nodes[node_index] = BvhNode { bbox, first: left as u32, count: 0 };

    let (lower, upper) = spheres.split_at_mut(mid);
    subdivide(nodes, left, lower, offset, rng);
    subdivide(nodes, left + 1, upper, offset + mid, rng);
}
