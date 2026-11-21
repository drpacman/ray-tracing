// GPU port of the CPU path tracer. One invocation renders one pixel:
// it traces `samples_per_pass` samples and adds them into `accum`.
// The host runs this repeatedly (one "pass" at a time, each with a new seed)
// and divides by the total sample count at the end.
//
// Each function mirrors the Rust function of the same name; see the Rust
// sources for the detailed comments and links to the book.

const PI: f32 = 3.14159265358979;
const T_MIN: f32 = 0.001;
const T_INFINITY: f32 = 3.4e38;
// Russian roulette (see ray_colour). Must match the constants in src/camera.rs.
const RR_MIN_DEPTH: u32 = 3u;
const RR_MAX_SURVIVAL: f32 = 0.95;

// --- Data shared with the host (layouts must match src/gpu.rs) ---------------

struct Params {
    pixel00: vec3f,
    samples_per_pass: u32,
    delta_u: vec3f,
    width: u32,
    delta_v: vec3f,
    height: u32,
    center: vec3f,
    max_depth: u32,
    defocus_u: vec3f,
    defocus_angle: f32,
    defocus_v: vec3f,
    seed: u32,
}

struct Sphere {
    center: vec3f,
    radius: f32,
    material: u32,
}

const MATERIAL_LAMBERTIAN: u32 = 0u;
const MATERIAL_METAL: u32 = 1u;
const MATERIAL_DIELECTRIC: u32 = 2u;

struct Material {
    albedo: vec3f,  // Metal
    fuzz: f32,      // Metal
    kind: u32,
    texture: u32,   // Lambertian: index into textures
    ior: f32,       // Dielectric
}

const TEXTURE_SOLID: u32 = 0u;
const TEXTURE_CHECKER: u32 = 1u;
const TEXTURE_IMAGE: u32 = 2u;

struct Texture {
    colour0: vec3f,     // Solid colour, or Checker odd
    inv_scale: f32,     // Checker
    colour1: vec3f,     // Checker even
    kind: u32,
    image_offset: u32,  // Image: first pixel in image_pixels
    image_width: u32,
    image_height: u32,
}

// Leaf (count > 0): spheres[first .. first + count].
// Interior (count == 0): children are bvh[first] and bvh[first + 1].
struct BvhNode {
    min: vec3f,
    first: u32,
    max: vec3f,
    count: u32,
}

@group(0) @binding(0) var<uniform> params: Params;
@group(0) @binding(1) var<storage, read> spheres: array<Sphere>;
@group(0) @binding(2) var<storage, read> materials: array<Material>;
@group(0) @binding(3) var<storage, read> textures: array<Texture>;
@group(0) @binding(4) var<storage, read> bvh: array<BvhNode>;
// All images' pixels, packed as RGBA8 in one u32 each.
@group(0) @binding(5) var<storage, read> image_pixels: array<u32>;
// Running sum of sample colours per pixel (w is unused).
@group(0) @binding(6) var<storage, read_write> accum: array<vec4f>;

// --- Random numbers (src/rng.rs) ---------------------------------------------

fn pcg_hash(input: u32) -> u32 {
    let state = input * 747796405u + 2891336453u;
    let word = ((state >> ((state >> 28u) + 4u)) ^ state) * 277803737u;
    return (word >> 22u) ^ word;
}

fn rand_u32(state: ptr<function, u32>) -> u32 {
    let s = *state;
    *state = s * 747796405u + 2891336453u;
    let word = ((s >> ((s >> 28u) + 4u)) ^ s) * 277803737u;
    return (word >> 22u) ^ word;
}

// Uniform in [0, 1).
fn rand_f32(state: ptr<function, u32>) -> f32 {
    return f32(rand_u32(state) >> 8u) * (1.0 / 16777216.0);
}

fn rand_range(state: ptr<function, u32>, min: f32, max: f32) -> f32 {
    return min + (max - min) * rand_f32(state);
}

fn random_in_unit_sphere(state: ptr<function, u32>) -> vec3f {
    var p: vec3f;
    loop {
        p = vec3f(rand_range(state, -1.0, 1.0), rand_range(state, -1.0, 1.0), rand_range(state, -1.0, 1.0));
        if (dot(p, p) < 1.0) {
            break;
        }
    }
    return p;
}

fn random_unit_vector(state: ptr<function, u32>) -> vec3f {
    return normalize(random_in_unit_sphere(state));
}

fn random_in_unit_disk(state: ptr<function, u32>) -> vec3f {
    var p: vec3f;
    loop {
        p = vec3f(rand_range(state, -1.0, 1.0), rand_range(state, -1.0, 1.0), 0.0);
        if (dot(p, p) < 1.0) {
            break;
        }
    }
    return p;
}

// --- Geometry (src/utils.rs, src/aabb.rs, src/scene.rs) -----------------------

struct Ray {
    orig: vec3f,
    dir: vec3f,
}

struct HitRecord {
    p: vec3f,
    normal: vec3f,
    t: f32,
    front_face: bool,
    material: u32,
    u: f32,
    v: f32,
}

// Returns t of the nearest hit in [t_min, t_max), or -1 for a miss.
fn sphere_intersect(s: Sphere, r: Ray, t_min: f32, t_max: f32) -> f32 {
    let oc = r.orig - s.center;
    let a = dot(r.dir, r.dir);
    let half_b = dot(oc, r.dir);
    // Precision-friendly discriminant, see Sphere::intersect.
    let l = oc - r.dir * (half_b / a);
    let discriminant = a * (s.radius * s.radius - dot(l, l));
    if (discriminant < 0.0) {
        return -1.0;
    }
    let sqrtd = sqrt(discriminant);
    var root = (-half_b - sqrtd) / a;
    if (root < t_min || root >= t_max) {
        root = (-half_b + sqrtd) / a;
        if (root < t_min || root >= t_max) {
            return -1.0;
        }
    }
    return root;
}

fn sphere_hit_record(s: Sphere, r: Ray, t: f32) -> HitRecord {
    let p = r.orig + r.dir * t;
    let outward_normal = (p - s.center) / s.radius;
    let front_face = dot(r.dir, outward_normal) < 0.0;
    let normal = select(-outward_normal, outward_normal, front_face);
    let theta = acos(-outward_normal.y);
    let phi = atan2(-outward_normal.z, outward_normal.x) + PI;
    return HitRecord(p, normal, t, front_face, s.material, phi / (2.0 * PI), theta / PI);
}

// The slab test, done on all three axes at once with vector min/max.
fn aabb_hit(node: BvhNode, r: Ray, inv_dir: vec3f, t_min: f32, t_max: f32) -> bool {
    let t0 = (node.min - r.orig) * inv_dir;
    let t1 = (node.max - r.orig) * inv_dir;
    let t_near = min(t0, t1);
    let t_far = max(t0, t1);
    let lo = max(t_min, max(t_near.x, max(t_near.y, t_near.z)));
    let hi = min(t_max, min(t_far.x, min(t_far.y, t_far.z)));
    return hi > lo;
}

// Closest hit via the BVH. Returns the sphere index (or -1) and its t.
struct Hit {
    t: f32,
    sphere: i32,
}

fn scene_hit(r: Ray, t_min: f32, t_max: f32) -> Hit {
    var stack: array<u32, 64>;
    var stack_len = 1u; // stack[0] = 0, the root
    var closest_t = t_max;
    var closest_sphere = -1;
    let inv_dir = 1.0 / r.dir;

    while (stack_len > 0u) {
        stack_len -= 1u;
        let node = bvh[stack[stack_len]];
        if (!aabb_hit(node, r, inv_dir, t_min, closest_t)) {
            continue;
        }
        if (node.count > 0u) {
            for (var i = node.first; i < node.first + node.count; i++) {
                let t = sphere_intersect(spheres[i], r, t_min, closest_t);
                if (t >= 0.0) {
                    closest_t = t;
                    closest_sphere = i32(i);
                }
            }
        } else {
            stack[stack_len] = node.first + 1u;
            stack[stack_len + 1u] = node.first;
            stack_len += 2u;
        }
    }
    return Hit(closest_t, closest_sphere);
}

// --- Textures and materials (src/material.rs) ---------------------------------

fn texture_value(index: u32, u: f32, v: f32, p: vec3f) -> vec3f {
    let tex = textures[index];
    if (tex.kind == TEXTURE_SOLID) {
        return tex.colour0;
    }
    if (tex.kind == TEXTURE_CHECKER) {
        let cell = vec3i(floor(tex.inv_scale * p));
        if ((cell.x + cell.y + cell.z) % 2 == 0) {
            return tex.colour1;
        }
        return tex.colour0;
    }
    // TEXTURE_IMAGE: nearest-pixel lookup, as on the CPU.
    let uu = clamp(u, 0.0, 1.0);
    let vv = 1.0 - clamp(v, 0.0, 1.0);
    let x = min(u32(uu * f32(tex.image_width)), tex.image_width - 1u);
    let y = min(u32(vv * f32(tex.image_height)), tex.image_height - 1u);
    let packed = image_pixels[tex.image_offset + y * tex.image_width + x];
    return unpack4x8unorm(packed).rgb;
}

fn reflectance(cos_theta: f32, refraction_ratio: f32) -> f32 {
    var r0 = (1.0 - refraction_ratio) / (1.0 + refraction_ratio);
    r0 = r0 * r0;
    return r0 + (1.0 - r0) * pow(1.0 - cos_theta, 5.0);
}

fn refract_dir(uv: vec3f, n: vec3f, etai_over_etat: f32) -> vec3f {
    let cos_theta = min(1.0, dot(-uv, n));
    let r_out_perp = (uv + n * cos_theta) * etai_over_etat;
    let r_out_parallel = n * -sqrt(abs(1.0 - dot(r_out_perp, r_out_perp)));
    return r_out_perp + r_out_parallel;
}

struct Scatter {
    attenuation: vec3f,
    ray: Ray,
}

// None of the current materials absorb a ray outright, so scatter always
// produces a ray (the CPU version's None case never happens).
fn scatter(r: Ray, rec: HitRecord, state: ptr<function, u32>) -> Scatter {
    let m = materials[rec.material];
    if (m.kind == MATERIAL_LAMBERTIAN) {
        var direction = rec.normal + random_unit_vector(state);
        if (all(abs(direction) < vec3f(1e-8))) {
            direction = rec.normal;
        }
        return Scatter(texture_value(m.texture, rec.u, rec.v, rec.p), Ray(rec.p, direction));
    }
    if (m.kind == MATERIAL_METAL) {
        let reflected = reflect(normalize(r.dir), rec.normal);
        let fuzzed = reflected + random_in_unit_sphere(state) * m.fuzz;
        return Scatter(m.albedo, Ray(rec.p, fuzzed));
    }
    // MATERIAL_DIELECTRIC
    let refraction_ratio = select(m.ior, 1.0 / m.ior, rec.front_face);
    let unit_direction = normalize(r.dir);
    let cos_theta = min(1.0, -dot(unit_direction, rec.normal));
    let sin_theta = sqrt(1.0 - cos_theta * cos_theta);
    let cannot_refract = refraction_ratio * sin_theta > 1.0;
    var direction: vec3f;
    if (cannot_refract || reflectance(cos_theta, refraction_ratio) > rand_f32(state)) {
        direction = reflect(unit_direction, rec.normal);
    } else {
        direction = refract_dir(unit_direction, rec.normal, refraction_ratio);
    }
    return Scatter(vec3f(1.0), Ray(rec.p, direction));
}

// --- Camera (src/camera.rs) ---------------------------------------------------

fn get_ray(i: u32, j: u32, state: ptr<function, u32>) -> Ray {
    let pixel_center = params.pixel00 + params.delta_u * f32(j) + params.delta_v * f32(i);
    let px = -0.5 + rand_f32(state);
    let py = -0.5 + rand_f32(state);
    let pixel_sample = pixel_center + params.delta_u * px + params.delta_v * py;
    var origin = params.center;
    if (params.defocus_angle > 0.0) {
        let d = random_in_unit_disk(state);
        origin = params.center + params.defocus_u * d.x + params.defocus_v * d.y;
    }
    return Ray(origin, pixel_sample - origin);
}

fn ray_colour(initial: Ray, state: ptr<function, u32>) -> vec3f {
    var r = initial;
    var throughput = vec3f(1.0);
    for (var depth = 0u; depth < params.max_depth; depth++) {
        let hit = scene_hit(r, T_MIN, T_INFINITY);
        if (hit.sphere < 0) {
            let unit_direction = normalize(r.dir);
            let t = 0.5 * (unit_direction.y + 1.0);
            let background = vec3f(1.0) * (1.0 - t) + vec3f(0.5, 0.7, 1.0) * t;
            return throughput * background;
        }
        let rec = sphere_hit_record(spheres[hit.sphere], r, hit.t);
        let s = scatter(r, rec, state);
        throughput *= s.attenuation;
        r = s.ray;

        // Russian roulette: randomly end low-throughput paths, scaling up the
        // survivors by 1/p so the average is unchanged (see camera.rs).
        if (depth >= RR_MIN_DEPTH) {
            let p = min(max(throughput.r, max(throughput.g, throughput.b)), RR_MAX_SURVIVAL);
            if (rand_f32(state) >= p) {
                return vec3f(0.0);
            }
            throughput /= p;
        }
    }
    return vec3f(0.0);
}

@compute @workgroup_size(8, 8)
fn main(@builtin(global_invocation_id) gid: vec3u) {
    // x = column j, y = row i. The grid is rounded up to whole workgroups.
    if (gid.x >= params.width || gid.y >= params.height) {
        return;
    }
    let index = gid.y * params.width + gid.x;
    // Same seeding as Rng::for_pixel.
    var state = pcg_hash(index ^ pcg_hash(params.seed + 0x9E3779B9u));
    var colour = vec3f(0.0);
    for (var s = 0u; s < params.samples_per_pass; s++) {
        colour += ray_colour(get_ray(gid.y, gid.x, &state), &state);
    }
    accum[index] += vec4f(colour, 0.0);
}
