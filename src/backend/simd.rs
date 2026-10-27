//! Batch products with explicit SIMD instructions, chosen when the program runs.

use std::any::Any;
use std::ops::Range;

use fearless_simd::{Level, Simd, SimdBase, SimdFloat, SimdFloatElement, dispatch};

use super::Backend;
use crate::batch::{Batch, Inputs, Product};
use crate::blade::BladeSet;
use crate::coefficient::Coefficient;
use crate::metric::Metric;
use crate::plan::{BySign, Entry, ProductPlan, SandwichPlan};

/// Runs batch products with explicit SIMD instructions, on the widest instruction set the CPU
/// supports: AVX-512 or AVX2 with FMA, falling back to SSE4.2 or SSE2, on x86; NEON on ARM; and
/// SIMD128 on WebAssembly built with it. With the `rayon` feature it can also split batches
/// across all cores.
///
/// Its buffers are plain [`Batch`] values, as for [`CpuBackend`](crate::CpuBackend), and it runs
/// the same cached plans. Each output coefficient is summed in registers for several vectors of
/// elements at once, with fused multiply-adds where the instruction set has them, so results
/// agree with the CPU backend to within a few rounding steps of the terms they add up rather than
/// bit for bit. Batches of other coefficient types than `f32` and `f64` run the CPU backend's
/// loops.
///
/// # Examples
///
/// ```
/// use seagrass_math::prelude::*;
///
/// let motor = Motor::<f32, Pga3>::new_translation_motor(1.0, 0.0, 0.0);
/// let motors = Batch::from_values(Pga3::default(), &[motor]);
/// let points = Batch::<f32, Pga3>::from_points(&[Vector3::new(0.0, 2.0, 0.0)]);
/// let simd = SimdBackend::default();
/// let moved = simd.sandwich(&motors, &points);
/// assert_eq!(moved.points_to_euclidean(), [Vector3::new(1.0, 2.0, 0.0)]);
/// ```
#[derive(Debug, Clone, Copy)]
pub struct SimdBackend {
    level: Level,
    parallel: bool,
}

impl SimdBackend {
    /// A backend that runs on the calling thread.
    pub fn single_threaded() -> Self {
        Self {
            level: Level::new(),
            parallel: false,
        }
    }

    /// A backend that splits batches across all cores. Without the `rayon` feature it runs on the
    /// calling thread.
    pub fn parallel() -> Self {
        Self {
            level: Level::new(),
            parallel: true,
        }
    }

    /// The name of the instruction set the backend runs on, such as `"AVX2"` or `"NEON"`.
    pub fn instruction_set(&self) -> &'static str {
        match self.level {
            #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
            Level::Avx512(_) => "AVX-512",
            #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
            Level::Avx2(_) => "AVX2",
            #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
            Level::Sse4_2(_) => "SSE4.2",
            #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
            Level::Sse2(_) => "SSE2",
            #[cfg(target_arch = "aarch64")]
            Level::Neon(_) => "NEON",
            #[cfg(all(target_arch = "wasm32", target_feature = "simd128"))]
            Level::WasmSimd128(_) => "SIMD128",
            _ => "none",
        }
    }
}

impl Default for SimdBackend {
    fn default() -> Self {
        Self::single_threaded()
    }
}

impl<T: Coefficient> Backend<T> for SimdBackend {
    type Buffer<M: Metric> = Batch<T, M>;

    fn upload<M: Metric>(&self, batch: &Batch<T, M>) -> Batch<T, M> {
        batch.clone()
    }

    fn download<M: Metric>(&self, buffer: &Batch<T, M>) -> Batch<T, M> {
        buffer.clone()
    }

    fn product<M: Metric>(&self, op: Product, lhs: &Batch<T, M>, rhs: &Batch<T, M>) -> Batch<T, M> {
        let mut out = Batch::zeros(lhs.metric(), M::Bits::EMPTY, 0);
        self.product_into(op, lhs, rhs, &mut out);
        out
    }

    fn sandwich<M: Metric>(&self, versors: &Batch<T, M>, operands: &Batch<T, M>) -> Batch<T, M> {
        let mut out = Batch::zeros(versors.metric(), M::Bits::EMPTY, 0);
        self.sandwich_into(versors, operands, &mut out);
        out
    }

    fn product_into<M: Metric>(
        &self,
        op: Product,
        lhs: &Batch<T, M>,
        rhs: &Batch<T, M>,
        out: &mut Batch<T, M>,
    ) {
        if let (Some(lhs), Some(rhs), Some(out)) = (cast(lhs), cast(rhs), cast_mut(out)) {
            self.simd_product_into::<f32, M>(op, lhs, rhs, out);
        } else if let (Some(lhs), Some(rhs), Some(out)) = (cast(lhs), cast(rhs), cast_mut(out)) {
            self.simd_product_into::<f64, M>(op, lhs, rhs, out);
        } else {
            lhs.product_into(op, rhs, out, self.parallel);
        }
    }

    fn sandwich_into<M: Metric>(
        &self,
        versors: &Batch<T, M>,
        operands: &Batch<T, M>,
        out: &mut Batch<T, M>,
    ) {
        if let (Some(m), Some(x), Some(out)) = (cast(versors), cast(operands), cast_mut(out)) {
            self.simd_sandwich_into::<f32, M>(m, x, out);
        } else if let (Some(m), Some(x), Some(out)) = (cast(versors), cast(operands), cast_mut(out))
        {
            self.simd_sandwich_into::<f64, M>(m, x, out);
        } else {
            versors.sandwich_into(operands, out, self.parallel);
        }
    }

    fn finish(&self) {}
}

impl SimdBackend {
    fn simd_product_into<T: Lanes, M: Metric>(
        &self,
        op: Product,
        lhs: &Batch<T, M>,
        rhs: &Batch<T, M>,
        out: &mut Batch<T, M>,
    ) {
        let level = self.level;
        lhs.product_into_with(op, rhs, out, self.parallel, |plan, inputs, range, out, offset| {
            dispatch!(level, simd => product_block::<_, T>(simd, plan, inputs, range, out, offset));
        });
    }

    fn simd_sandwich_into<T: Lanes, M: Metric>(
        &self,
        versors: &Batch<T, M>,
        operands: &Batch<T, M>,
        out: &mut Batch<T, M>,
    ) {
        let level = self.level;
        versors.sandwich_into_with(operands, out, self.parallel, |plan, inputs, range, out, offset| {
            dispatch!(level, simd => sandwich_block::<_, T>(simd, plan, inputs, range, out, offset));
        });
    }
}

/// `batch` as a batch of `U`, if `T` is `U`. The type check is resolved at compile time.
fn cast<T: Coefficient, U: Coefficient, M: Metric>(batch: &Batch<T, M>) -> Option<&Batch<U, M>> {
    (batch as &dyn Any).downcast_ref()
}

/// `batch` as a batch of `U`, if `T` is `U`.
fn cast_mut<T: Coefficient, U: Coefficient, M: Metric>(
    batch: &mut Batch<T, M>,
) -> Option<&mut Batch<U, M>> {
    (batch as &mut dyn Any).downcast_mut()
}

/// A coefficient type with SIMD kernels, and its native-width vector at each instruction set.
trait Lanes: Coefficient + SimdFloatElement {
    type Vector<S: Simd>: SimdFloat<S, Element = Self>;
}

impl Lanes for f32 {
    type Vector<S: Simd> = S::f32s;
}

impl Lanes for f64 {
    type Vector<S: Simd> = S::f64s;
}

/// Vectors of elements a product computes together: each output keeps one sum per vector, so
/// each entry's plane lookups are shared by this many vectors.
const PRODUCT_GROUP: usize = 8;

/// Vectors of elements a sandwich computes together. Each output keeps two sums per vector, one
/// per sign, which is as many registers as a product uses; that measured faster for sandwiches
/// than one sum per vector over twice as many vectors.
const SANDWICH_GROUP: usize = 4;

// Every function below that takes a `Simd` token is `#[inline(always)]`, so it is compiled into
// the `dispatch!` call with that instruction set's target features enabled.

/// One block of an elementwise bilinear product. Each output is computed for the whole block
/// before the next, a group of vectors at a time, so the block's inputs stay in the cache and each
/// output plane is written in one run.
#[inline(always)]
fn product_block<S: Simd, T: Lanes>(
    simd: S,
    plan: &ProductPlan,
    inputs: Inputs<'_, T>,
    range: Range<usize>,
    out: &mut [&mut [T]],
    offset: usize,
) {
    let plan = plan.by_sign();
    let width = PRODUCT_GROUP * T::Vector::<S>::LEN;
    let whole = range.len() / width * width;
    let groups = range.start..range.start + whole;
    product_groups::<S, T>(simd, plan, inputs, groups, out, offset);
    if whole < range.len() {
        let mut tail = Tail::new(inputs, range.start + whole..range.end, out.len(), width);
        let (inputs, mut group) = tail.split();
        product_groups::<S, T>(simd, plan, inputs, 0..width, &mut group, 0);
        tail.copy_to(out, offset + whole);
    }
}

/// One block of an elementwise sandwich, output by output as for products.
#[inline(always)]
fn sandwich_block<S: Simd, T: Lanes>(
    simd: S,
    plan: &SandwichPlan,
    inputs: Inputs<'_, T>,
    range: Range<usize>,
    out: &mut [&mut [T]],
    offset: usize,
) {
    let width = SANDWICH_GROUP * T::Vector::<S>::LEN;
    let whole = range.len() / width * width;
    let groups = range.start..range.start + whole;
    sandwich_groups::<S, T>(simd, plan, inputs, groups, out, offset);
    if whole < range.len() {
        let mut tail = Tail::new(inputs, range.start + whole..range.end, out.len(), width);
        let (inputs, mut group) = tail.split();
        sandwich_groups::<S, T>(simd, plan, inputs, 0..width, &mut group, 0);
        tail.copy_to(out, offset + whole);
    }
}

/// The last group of a batch whose length is not a whole number of groups, copied into planes of
/// a whole group's width padded with zeros, so it runs the same kernel as every other group.
struct Tail<T> {
    lhs: Vec<T>,
    rhs: Vec<T>,
    out: Vec<T>,
    len: usize,
    width: usize,
}

impl<T: Coefficient> Tail<T> {
    fn new(
        (lhs, rhs, len): Inputs<'_, T>,
        range: Range<usize>,
        outputs: usize,
        width: usize,
    ) -> Self {
        let padded = |planes: &[T]| {
            let mut padded = vec![T::zero(); planes.len() / len * width];
            for (plane, padded) in planes.chunks(len).zip(padded.chunks_mut(width)) {
                padded[..range.len()].copy_from_slice(&plane[range.clone()]);
            }
            padded
        };
        Self {
            lhs: padded(lhs),
            rhs: padded(rhs),
            out: vec![T::zero(); outputs * width],
            len: range.len(),
            width,
        }
    }

    /// The padded inputs, and the output planes to compute into.
    fn split(&mut self) -> (Inputs<'_, T>, Vec<&mut [T]>) {
        let inputs = (self.lhs.as_slice(), self.rhs.as_slice(), self.width);
        (inputs, self.out.chunks_mut(self.width).collect())
    }

    /// Copies the computed elements into `out` from `at`.
    fn copy_to(&self, out: &mut [&mut [T]], at: usize) {
        for (out, group) in out.iter_mut().zip(self.out.chunks(self.width)) {
            out[at..at + self.len].copy_from_slice(&group[..self.len]);
        }
    }
}

/// The groups of an elementwise bilinear product whose elements are `elements` in the input
/// planes, written from `at` in the output planes. Each output sums its negative entries,
/// negates the sum, and adds its positive entries, so no entry branches on its sign.
#[inline(always)]
fn product_groups<S: Simd, T: Lanes>(
    simd: S,
    plan: &BySign,
    inputs: Inputs<'_, T>,
    elements: Range<usize>,
    out: &mut [&mut [T]],
    at: usize,
) {
    let lanes = T::Vector::<S>::LEN;
    let width = PRODUCT_GROUP * lanes;
    for ([negative, positive], out) in plan.outputs().zip(out.iter_mut()) {
        let out = &mut out[at..][..elements.len()];
        for (start, out) in elements
            .clone()
            .step_by(width)
            .zip(out.chunks_exact_mut(width))
        {
            let mut sums = [T::Vector::<S>::splat(simd, T::zero()); PRODUCT_GROUP];
            add_products(simd, &mut sums, negative, inputs, start);
            for sum in &mut sums {
                *sum = -*sum;
            }
            add_products(simd, &mut sums, positive, inputs, start);
            for (sum, out) in sums.iter().zip(out.chunks_exact_mut(lanes)) {
                sum.store_slice(out);
            }
        }
    }
}

/// Adds `lhs[entry.lhs] * rhs[entry.rhs]` for every entry to every vector of the group.
#[inline(always)]
fn add_products<S: Simd, T: Lanes>(
    simd: S,
    sums: &mut [T::Vector<S>; PRODUCT_GROUP],
    entries: &[Entry],
    (lhs, rhs, len): Inputs<'_, T>,
    start: usize,
) {
    let width = PRODUCT_GROUP * T::Vector::<S>::LEN;
    for entry in entries {
        let a = &lhs[usize::from(entry.lhs) * len + start..][..width];
        let b = &rhs[usize::from(entry.rhs) * len + start..][..width];
        for (g, sum) in sums.iter_mut().enumerate() {
            let (a, b) = (load::<S, T>(simd, a, g), load::<S, T>(simd, b, g));
            *sum = a.mul_add(b, *sum);
        }
    }
}

/// The groups of an elementwise sandwich whose elements are `elements` in the input planes,
/// written from `at` in the output planes.
#[inline(always)]
fn sandwich_groups<S: Simd, T: Lanes>(
    simd: S,
    plan: &SandwichPlan,
    (m, x, len): Inputs<'_, T>,
    elements: Range<usize>,
    out: &mut [&mut [T]],
    at: usize,
) {
    let lanes = T::Vector::<S>::LEN;
    let width = SANDWICH_GROUP * lanes;
    let zero = T::Vector::<S>::splat(simd, T::zero());
    for (terms, out) in plan.outputs().zip(out.iter_mut()) {
        let out = &mut out[at..][..elements.len()];
        for (start, out) in elements
            .clone()
            .step_by(width)
            .zip(out.chunks_exact_mut(width))
        {
            let (mut plus, mut minus) = ([zero; SANDWICH_GROUP], [zero; SANDWICH_GROUP]);
            for term in terms {
                let factors = (
                    &m[usize::from(term.first) * len + start..][..width],
                    &x[usize::from(term.operand) * len + start..][..width],
                    &m[usize::from(term.second) * len + start..][..width],
                );
                // A merged term sums at most two triples: the coefficient is -2, -1, 1 or 2.
                match term.coefficient {
                    1 => add_triple::<S, T, false>(simd, &mut plus, factors),
                    2 => add_triple::<S, T, true>(simd, &mut plus, factors),
                    -1 => add_triple::<S, T, false>(simd, &mut minus, factors),
                    _ => add_triple::<S, T, true>(simd, &mut minus, factors),
                }
            }
            for ((plus, minus), out) in plus.iter().zip(&minus).zip(out.chunks_exact_mut(lanes)) {
                (*plus - *minus).store_slice(out);
            }
        }
    }
}

/// Adds `f * x * s`, doubled if `DOUBLE`, to every vector of the group.
#[inline(always)]
fn add_triple<S: Simd, T: Lanes, const DOUBLE: bool>(
    simd: S,
    sums: &mut [T::Vector<S>; SANDWICH_GROUP],
    (f, x, s): (&[T], &[T], &[T]),
) {
    for (g, sum) in sums.iter_mut().enumerate() {
        let fx = load::<S, T>(simd, f, g) * load::<S, T>(simd, x, g);
        // Doubling is exact, so it adds no rounding.
        let fx = if DOUBLE { fx + fx } else { fx };
        *sum = fx.mul_add(load::<S, T>(simd, s, g), *sum);
    }
}

/// Vector `g` of a group.
#[inline(always)]
fn load<S: Simd, T: Lanes>(simd: S, group: &[T], g: usize) -> T::Vector<S> {
    let lanes = T::Vector::<S>::LEN;
    T::Vector::<S>::from_slice(simd, &group[g * lanes..][..lanes])
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::metric::Dynamic;

    /// Every instruction set this CPU supports, from the baseline up.
    fn levels() -> Vec<Level> {
        #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
        {
            let level = Level::new();
            let mut levels: Vec<Level> = level.as_sse2().map(Level::Sse2).into_iter().collect();
            levels.extend(level.as_sse4_2().map(Level::Sse4_2));
            levels.extend(level.as_avx2().map(Level::Avx2));
            levels.extend(level.as_avx512().map(Level::Avx512));
            levels
        }
        #[cfg(not(any(target_arch = "x86", target_arch = "x86_64")))]
        vec![Level::new()]
    }

    /// A deterministic generator of blade sets and of values in [-2, 2).
    struct Values(u64);

    impl Values {
        fn next(&mut self) -> u64 {
            self.0 ^= self.0 << 13;
            self.0 ^= self.0 >> 7;
            self.0 ^= self.0 << 17;
            self.0
        }

        fn batch<T: Coefficient, M: Metric>(
            &mut self,
            metric: M,
            blades: M::Bits,
            len: usize,
        ) -> Batch<T, M> {
            let planes = (0..blades.count_ones() as usize * len)
                .map(|_| T::from((self.next() % 4001) as f64 / 1000.0 - 2.0).unwrap())
                .collect();
            Batch::from_planes(metric, blades, len, planes).unwrap()
        }
    }

    /// The sum of the absolute coefficients of element `i`.
    fn l1<T: Coefficient, M: Metric>(batch: &Batch<T, M>, i: usize) -> T {
        let planes = batch.blades().count_ones() as usize;
        (0..planes)
            .map(|k| batch.planes()[k * batch.len() + i].abs())
            .fold(T::zero(), |sum, c| sum + c)
    }

    /// Each coefficient agrees to `4 * epsilon` relative to `scale(i)`, a bound on the sum of
    /// the absolute values of the terms it adds up. Summing in a different order, or with fused
    /// multiply-adds, changes a sum by a few rounding steps of that size, however much the terms
    /// cancel.
    #[track_caller]
    fn assert_close<T: Coefficient, M: Metric>(
        simd: &Batch<T, M>,
        cpu: &Batch<T, M>,
        scale: impl Fn(usize) -> T,
        what: &str,
    ) {
        assert_eq!(simd.blades(), cpu.blades(), "{what}: stored blades");
        let tolerance = T::epsilon() * T::from(4).unwrap();
        let planes = cpu.blades().count_ones() as usize;
        for i in 0..cpu.len() {
            let scale = scale(i).max(T::one());
            for k in 0..planes {
                let (s, c) = (
                    simd.planes()[k * cpu.len() + i],
                    cpu.planes()[k * cpu.len() + i],
                );
                assert!(
                    (s - c).abs() <= tolerance * scale,
                    "{what}: element {i}, slot {k}: {s:?} != {c:?}"
                );
            }
        }
    }

    const PRODUCTS: [Product; 5] = [
        Product::Geometric,
        Product::Outer,
        Product::Inner,
        Product::LeftContraction,
        Product::RightContraction,
    ];

    fn random_structures<T: Coefficient, B: BladeSet>(backend: SimdBackend, seed: u64) {
        let mut values = Values(seed);
        let dim = B::DIM;
        for p in (0..=dim).step_by(2) {
            let metric =
                Dynamic::<B>::new(p, dim - p - u8::from(p < dim), u8::from(p < dim)).unwrap();
            // A batch shorter than one group, one that ends partway through a group, and one
            // longer than a rayon chunk.
            for len in [5, 2 * 256 + 17, 4_099] {
                let (lhs, rhs) = (B::from_u64(values.next()), B::from_u64(values.next()));
                let a: Batch<T, _> = values.batch(metric, lhs, len);
                let b: Batch<T, _> = values.batch(metric, rhs, len);
                let what = |op: &str| {
                    format!(
                        "{} {op} in {metric} {lhs:?} x {rhs:?}, {len} elements",
                        backend.instruction_set()
                    )
                };
                for op in PRODUCTS {
                    let scale = |i| l1(&a, i) * l1(&b, i);
                    let simd = backend.product(op, &a, &b);
                    assert_close(&simd, &a.product(op, &b), scale, &what(&format!("{op:?}")));
                }
                let scale = |i| l1(&a, i) * l1(&a, i) * l1(&b, i);
                let simd = backend.sandwich(&a, &b);
                assert_close(&simd, &a.sandwich(&b), scale, &what("sandwich"));
            }
        }
    }

    #[test]
    fn every_instruction_set_matches_the_cpu() {
        for level in levels() {
            for parallel in [false, true] {
                let backend = SimdBackend { level, parallel };
                random_structures::<f32, u8>(backend, 1);
                random_structures::<f32, u16>(backend, 2);
                random_structures::<f64, u16>(backend, 3);
                random_structures::<f32, u32>(backend, 4);
                random_structures::<f64, u64>(backend, 5);
            }
        }
    }

    #[test]
    fn empty_batches_and_empty_structures() {
        let backend = SimdBackend::default();
        let metric = Dynamic::<u16>::new(3, 0, 1).unwrap();
        let empty = Batch::<f32, _>::zeros(metric, 0x9669, 0);
        assert!(backend.sandwich(&empty, &empty).is_empty());
        // No stored blades: every product has an empty structure, at every length.
        let none = Batch::<f64, _>::zeros(metric, 0, 37);
        let product = backend.product(Product::Geometric, &none, &none);
        assert_eq!((product.len(), product.blades()), (37, 0));
    }
}
