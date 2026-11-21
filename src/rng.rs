// A tiny PCG random number generator (32-bit state, RXS-M-XS output).
// It is fast, has no global/thread-local state, and is trivial to port to a
// GPU shader, where each pixel gets its own generator seeded from its index.
// See https://www.pcg-random.org and
// https://www.reedbeta.com/blog/hash-functions-for-gpu-rendering/
#[derive(Clone)]
pub struct Rng {
    state: u32
}

// Stateless PCG hash: scrambles an integer into a well-distributed one.
pub fn pcg_hash(input: u32) -> u32 {
    let state = input.wrapping_mul(747796405).wrapping_add(2891336453);
    let word = ((state >> ((state >> 28) + 4)) ^ state).wrapping_mul(277803737);
    (word >> 22) ^ word
}

impl Rng {
    pub fn new(seed: u32) -> Self {
        Self { state: pcg_hash(seed) }
    }

    // A generator for one pixel. Different seeds give independent images,
    // e.g. for accumulating several passes.
    pub fn for_pixel(pixel_index: u32, seed: u32) -> Self {
        Self::new(pixel_index ^ pcg_hash(seed.wrapping_add(0x9E3779B9)))
    }

    pub fn next_u32(&mut self) -> u32 {
        let state = self.state;
        self.state = state.wrapping_mul(747796405).wrapping_add(2891336453);
        let word = ((state >> ((state >> 28) + 4)) ^ state).wrapping_mul(277803737);
        (word >> 22) ^ word
    }

    // Uniform in [0, 1): use the top 24 bits, which is all the precision an f32 mantissa holds.
    pub fn next_f32(&mut self) -> f32 {
        (self.next_u32() >> 8) as f32 * (1.0 / 16_777_216.0)
    }

    // Uniform in [min, max).
    pub fn range(&mut self, min: f32, max: f32) -> f32 {
        min + (max - min) * self.next_f32()
    }

    // Uniform integer in [0, n).
    pub fn below(&mut self, n: u32) -> u32 {
        self.next_u32() % n
    }
}
