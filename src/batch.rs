//! Many multivectors that share one set of stored blades, kept as one array per blade.

use std::ops::Range;

#[cfg(feature = "nalgebra")]
use nalgebra::{Matrix4, Vector3};

use crate::blade::{BladeIndex, BladeSet};
use crate::coefficient::Coefficient;
use crate::error::BladeError;
use crate::metric::Metric;
#[cfg(feature = "nalgebra")]
use crate::metric::Pga3;
use crate::multivector::Multivector;
use crate::plan::{Op, ProductPlan, SandwichPlan, with_product_plan, with_sandwich_plan};
#[cfg(feature = "nalgebra")]
use crate::typed::Motor;

/// The bilinear products a batch can compute elementwise.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Product {
    /// The geometric product, `a * b`.
    Geometric,
    /// The outer product, `a ^ b`.
    Outer,
    /// The inner product, `a | b`.
    Inner,
    /// The left contraction, `a << b`.
    LeftContraction,
    /// The right contraction, `a >> b`.
    RightContraction,
}

impl Product {
    pub(crate) fn op(self) -> Op {
        match self {
            Self::Geometric => Op::Geometric,
            Self::Outer => Op::Outer,
            Self::Inner => Op::Inner,
            Self::LeftContraction => Op::LeftContraction,
            Self::RightContraction => Op::RightContraction,
        }
    }
}

/// Elements processed together, so that a block of every plane stays in the CPU cache.
const BLOCK: usize = 256;

/// The inputs of a block kernel: the left and right operands' planes, and their length.
pub(crate) type Inputs<'a, T> = (&'a [T], &'a [T], usize);

/// Many multivectors of one algebra that store the same blades, kept as one array per blade
/// ("structure of arrays"): the coefficient of blade slot k of element i is
/// `planes()[k * len() + i]`.
///
/// Elementwise products run one cached plan over whole arrays, so they compute exactly what the
/// single-value products compute, in the same order, with no per-element plan lookup. A batch has
/// one fixed structure: values that store fewer blades are padded with zeros, and products keep
/// every blade their plan can produce, including the blades a single-value sandwich would drop as
/// zero.
///
/// # Examples
///
/// ```
/// use seagrass_math::prelude::*;
///
/// let motor = Motor::<f64, Pga3>::new_translation_motor(0.0, 0.0, 1.0);
/// let points = [
///     Multivector::<f64, Pga3>::new_point(1.0, 0.0, 0.0),
///     Multivector::<f64, Pga3>::new_point(0.0, 2.0, 0.0),
/// ];
/// let motors = Batch::from_values(Pga3::default(), &[motor.clone(), motor]);
/// let points = Batch::from_values(Pga3::default(), &points);
/// let moved = motors.sandwich(&points);
/// assert_eq!(moved.get(1).point_to_euclidean(), Vector3::new(0.0, 2.0, 1.0));
/// ```
#[derive(Debug, Clone, PartialEq)]
pub struct Batch<T, M: Metric> {
    metric: M,
    blades: M::Bits,
    len: usize,
    planes: Vec<T>,
}

impl<T: Coefficient, M: Metric> Batch<T, M> {
    /// A batch of `len` zeros that store `blades`.
    pub fn zeros(metric: M, blades: M::Bits, len: usize) -> Self {
        let planes = vec![T::zero(); blades.count_ones() as usize * len];
        Self {
            metric,
            blades,
            len,
            planes,
        }
    }

    /// The values in order, storing the union of their blades; a value that does not store a
    /// blade of the union gets a zero there.
    ///
    /// # Panics
    ///
    /// Panics if a value's runtime metric differs from `metric`.
    pub fn from_values<V: AsRef<Multivector<T, M>>>(metric: M, values: &[V]) -> Self {
        let blades = values
            .iter()
            .fold(M::Bits::EMPTY, |union, v| union | v.as_ref().blades());
        let mut batch = Self::zeros(metric, blades, values.len());
        for (i, value) in values.iter().enumerate() {
            batch
                .set(i, value.as_ref())
                .expect("every value's blades are in the union");
        }
        batch
    }

    /// The batch with the given blades and planes, `planes[k * len + i]` holding blade slot k of
    /// element i.
    ///
    /// # Errors
    ///
    /// [`BladeError::CoefficientCount`] if `planes` does not hold one coefficient per blade and
    /// element.
    pub fn from_planes(
        metric: M,
        blades: M::Bits,
        len: usize,
        planes: Vec<T>,
    ) -> Result<Self, BladeError> {
        let expected = blades.count_ones() as usize * len;
        if planes.len() != expected {
            return Err(BladeError::CoefficientCount {
                expected,
                found: planes.len(),
            });
        }
        Ok(Self {
            metric,
            blades,
            len,
            planes,
        })
    }

    /// The metric.
    pub fn metric(&self) -> M {
        self.metric
    }

    /// The blades every element stores.
    pub fn blades(&self) -> M::Bits {
        self.blades
    }

    /// The number of elements.
    pub fn len(&self) -> usize {
        self.len
    }

    /// Whether the batch has no elements.
    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// All coefficients, plane by plane.
    pub fn planes(&self) -> &[T] {
        &self.planes
    }

    /// All coefficients, plane by plane, mutably.
    pub fn planes_mut(&mut self) -> &mut [T] {
        &mut self.planes
    }

    /// The coefficients of blade `index` for every element, or `None` if the batch does not store
    /// that blade.
    pub fn plane(&self, index: BladeIndex) -> Option<&[T]> {
        self.blades.has(index).then(|| {
            let slot = self.blades.count_below(index) as usize;
            &self.planes[slot * self.len..(slot + 1) * self.len]
        })
    }

    /// Element `i` as a multivector that stores every blade of the batch.
    ///
    /// # Panics
    ///
    /// Panics if `i` is out of bounds.
    pub fn get(&self, i: usize) -> Multivector<T, M> {
        assert!(
            i < self.len,
            "index {i} out of bounds for a batch of {}",
            self.len
        );
        let coefficients = (0..self.blades.count_ones() as usize)
            .map(|k| self.planes[k * self.len + i])
            .collect();
        Multivector::from_raw(self.metric, self.blades, coefficients)
    }

    /// Replaces element `i` with `value`, writing zeros to the blades it does not store.
    ///
    /// # Errors
    ///
    /// [`BladeError::NotStored`] if `value` stores a blade the batch does not.
    ///
    /// # Panics
    ///
    /// Panics if `i` is out of bounds or the metrics differ.
    pub fn set(&mut self, i: usize, value: &Multivector<T, M>) -> Result<(), BladeError> {
        assert!(
            i < self.len,
            "index {i} out of bounds for a batch of {}",
            self.len
        );
        assert_eq!(
            self.metric,
            value.metric(),
            "the value has a different metric"
        );
        if let Some(index) = (value.blades() & !self.blades).iter().next() {
            return Err(BladeError::NotStored { index });
        }
        for (slot, index) in self.blades.iter().enumerate() {
            self.planes[slot * self.len + i] = value.get(index);
        }
        Ok(())
    }

    /// Every element as a multivector.
    pub fn to_values(&self) -> Vec<Multivector<T, M>> {
        (0..self.len).map(|i| self.get(i)).collect()
    }

    /// The elementwise product `op(self[i], rhs[i])`.
    ///
    /// # Panics
    ///
    /// Panics if the batches differ in length or runtime metric.
    pub fn product(&self, op: Product, rhs: &Self) -> Self {
        self.product_with(op, rhs, false)
    }

    /// The elementwise geometric product.
    ///
    /// # Panics
    ///
    /// Panics if the batches differ in length or runtime metric.
    pub fn geometric_product(&self, rhs: &Self) -> Self {
        self.product(Product::Geometric, rhs)
    }

    /// The elementwise sandwich `self[i] * x[i] * reverse(self[i])`.
    ///
    /// # Panics
    ///
    /// Panics if the batches differ in length or runtime metric.
    pub fn sandwich(&self, x: &Self) -> Self {
        self.sandwich_with(x, false)
    }

    pub(crate) fn assert_matches(&self, rhs: &Self) {
        assert_eq!(self.metric, rhs.metric, "batches have different metrics");
        assert_eq!(self.len, rhs.len, "batches have different lengths");
    }

    pub(crate) fn product_with(&self, op: Product, rhs: &Self, parallel: bool) -> Self {
        let mut out = Self::zeros(self.metric, M::Bits::EMPTY, 0);
        self.product_into(op, rhs, &mut out, parallel);
        out
    }

    pub(crate) fn sandwich_with(&self, x: &Self, parallel: bool) -> Self {
        let mut out = Self::zeros(self.metric, M::Bits::EMPTY, 0);
        self.sandwich_into(x, &mut out, parallel);
        out
    }

    /// Gives `out` these blades and length, reusing its allocation. The coefficients are left
    /// as they are; the kernels overwrite every one.
    fn reset(out: &mut Self, metric: M, blades: M::Bits, len: usize) {
        out.metric = metric;
        out.blades = blades;
        out.len = len;
        out.planes
            .resize(blades.count_ones() as usize * len, T::zero());
    }

    /// The elementwise product written into `out`, whose allocation is reused. Writing into a
    /// batch that is already the right size avoids allocating (and, for large batches, the page
    /// faults of touching fresh memory).
    pub(crate) fn product_into(&self, op: Product, rhs: &Self, out: &mut Self, parallel: bool) {
        self.product_into_with(op, rhs, out, parallel, product_block);
    }

    /// Like [`product_into`](Self::product_into), with `kernel` computing each block.
    pub(crate) fn product_into_with(
        &self,
        op: Product,
        rhs: &Self,
        out: &mut Self,
        parallel: bool,
        kernel: impl Fn(&ProductPlan, Inputs<'_, T>, Range<usize>, &mut [&mut [T]], usize) + Sync,
    ) {
        self.assert_matches(rhs);
        with_product_plan(op.op(), &self.metric, self.blades, rhs.blades, |plan| {
            Self::reset(out, self.metric, plan.output(), self.len);
            let inputs = (self.planes.as_slice(), rhs.planes.as_slice(), self.len);
            run_blocks(&mut out.planes, self.len, parallel, |range, out, offset| {
                kernel(plan, inputs, range, out, offset);
            });
        });
    }

    /// The elementwise sandwich written into `out`, whose allocation is reused.
    pub(crate) fn sandwich_into(&self, x: &Self, out: &mut Self, parallel: bool) {
        self.sandwich_into_with(x, out, parallel, sandwich_block);
    }

    /// Like [`sandwich_into`](Self::sandwich_into), with `kernel` computing each block.
    pub(crate) fn sandwich_into_with(
        &self,
        x: &Self,
        out: &mut Self,
        parallel: bool,
        kernel: impl Fn(&SandwichPlan, Inputs<'_, T>, Range<usize>, &mut [&mut [T]], usize) + Sync,
    ) {
        self.assert_matches(x);
        with_sandwich_plan(&self.metric, self.blades, x.blades, |plan| {
            Self::reset(out, self.metric, plan.output(), self.len);
            let inputs = (self.planes.as_slice(), x.planes.as_slice(), self.len);
            run_blocks(&mut out.planes, self.len, parallel, |range, out, offset| {
                kernel(plan, inputs, range, out, offset);
            });
        });
    }
}

/// Splits the elements into chunks, one for the calling thread or a few per core with the
/// `rayon` feature and `parallel` set, and runs `kernel` over each chunk in blocks of [`BLOCK`]
/// elements. `kernel` gets the block's element range, the chunk's slice of every output plane,
/// and the block's offset in those slices.
fn run_blocks<T: Coefficient>(
    out: &mut [T],
    len: usize,
    parallel: bool,
    kernel: impl Fn(Range<usize>, &mut [&mut [T]], usize) + Sync,
) {
    if len == 0 {
        return;
    }
    let chunk = chunk_len(len, parallel);
    let mut chunks: Vec<Vec<&mut [T]>> = (0..len.div_ceil(chunk)).map(|_| Vec::new()).collect();
    for plane in out.chunks_mut(len) {
        for (c, piece) in plane.chunks_mut(chunk).enumerate() {
            chunks[c].push(piece);
        }
    }
    let run = |(c, mut planes): (usize, Vec<&mut [T]>)| {
        let (start, end) = (c * chunk, ((c + 1) * chunk).min(len));
        let mut block = start;
        while block < end {
            let block_end = (block + BLOCK).min(end);
            kernel(block..block_end, &mut planes, block - start);
            block = block_end;
        }
    };
    #[cfg(feature = "rayon")]
    if parallel {
        use rayon::prelude::*;
        chunks.into_par_iter().enumerate().for_each(run);
        return;
    }
    chunks.into_iter().enumerate().for_each(run);
}

/// Elements per chunk: all of them on one thread, or about four chunks per core, in whole
/// blocks.
fn chunk_len(len: usize, parallel: bool) -> usize {
    #[cfg(feature = "rayon")]
    if parallel {
        let chunks = 4 * rayon::current_num_threads();
        return len.div_ceil(chunks).div_ceil(BLOCK).max(1) * BLOCK;
    }
    let _ = parallel;
    len
}

/// One block of an elementwise bilinear product. Each element sums its plan entries in plan
/// order, as the single-value product does, so the results agree bit for bit.
fn product_block<T: Coefficient>(
    plan: &ProductPlan,
    (lhs, rhs, len): Inputs<'_, T>,
    range: Range<usize>,
    out: &mut [&mut [T]],
    offset: usize,
) {
    let count = range.len();
    for (entries, out) in plan.outputs().zip(out.iter_mut()) {
        let out = &mut out[offset..offset + count];
        let Some((first, rest)) = entries.split_first() else {
            out.fill(T::zero());
            continue;
        };
        let slot = |planes, slot| block_of(planes, slot, len, &range);
        // The first entry starts the sum from zero, as the single-value product does.
        multiply_start(
            out,
            slot(lhs, first.lhs),
            slot(rhs, first.rhs),
            first.negative,
        );
        for entry in rest {
            let (a, b) = (slot(lhs, entry.lhs), slot(rhs, entry.rhs));
            if entry.negative {
                multiply_subtract(out, a, b);
            } else {
                multiply_add(out, a, b);
            }
        }
    }
}

/// One block of an elementwise sandwich, term by term as the single-value plan runs.
fn sandwich_block<T: Coefficient>(
    plan: &SandwichPlan,
    (m, x, len): Inputs<'_, T>,
    range: Range<usize>,
    out: &mut [&mut [T]],
    offset: usize,
) {
    let count = range.len();
    for (terms, out) in plan.outputs().zip(out.iter_mut()) {
        let out = &mut out[offset..offset + count];
        out.fill(T::zero());
        for term in terms {
            let first = &m[usize::from(term.first) * len + range.start..][..count];
            let operand = &x[usize::from(term.operand) * len + range.start..][..count];
            let second = &m[usize::from(term.second) * len + range.start..][..count];
            let c = T::from(term.coefficient).expect("small integers are representable");
            triple_multiply_add(out, (first, operand, second), term.coefficient, c);
        }
    }
}

// The loops below take their arrays as separate reference parameters, which tells the compiler
// they do not overlap, so it can vectorize them.

/// The coefficients of plane `slot` for the elements in `range`.
#[inline]
fn block_of<'a, T>(planes: &'a [T], slot: u8, len: usize, range: &Range<usize>) -> &'a [T] {
    &planes[usize::from(slot) * len + range.start..][..range.len()]
}

/// `out[i] = 0 + a[i] * b[i]`, or `0 - a[i] * b[i]` if `negative`: the first step of a sum.
#[inline]
fn multiply_start<T: Coefficient>(out: &mut [T], a: &[T], b: &[T], negative: bool) {
    let zero = T::zero();
    if negative {
        for ((o, a), b) in out.iter_mut().zip(a).zip(b) {
            *o = zero - *a * *b;
        }
    } else {
        for ((o, a), b) in out.iter_mut().zip(a).zip(b) {
            *o = zero + *a * *b;
        }
    }
}

/// `out[i] += a[i] * b[i]`.
#[inline]
fn multiply_add<T: Coefficient>(out: &mut [T], a: &[T], b: &[T]) {
    for ((o, a), b) in out.iter_mut().zip(a).zip(b) {
        *o += *a * *b;
    }
}

/// `out[i] -= a[i] * b[i]`.
#[inline]
fn multiply_subtract<T: Coefficient>(out: &mut [T], a: &[T], b: &[T]) {
    for ((o, a), b) in out.iter_mut().zip(a).zip(b) {
        *o -= *a * *b;
    }
}

/// `out[i] += c * (f[i] * x[i] * s[i])`, with the products for `c` of 1 and -1 added or
/// subtracted as the single-value plan does.
#[inline]
fn triple_multiply_add<T: Coefficient>(
    out: &mut [T],
    (f, x, s): (&[T], &[T], &[T]),
    coefficient: i8,
    c: T,
) {
    let terms = out.iter_mut().zip(f).zip(x).zip(s);
    match coefficient {
        1 => terms.for_each(|(((o, f), x), s)| *o += *f * *x * *s),
        -1 => terms.for_each(|(((o, f), x), s)| *o += -(*f * *x * *s)),
        _ => terms.for_each(|(((o, f), x), s)| *o += c * (*f * *x * *s)),
    }
}

#[cfg(feature = "nalgebra")]
impl<T: Coefficient> Batch<T, Pga3> {
    /// The points (x, y, z), as PGA3D trivectors.
    ///
    /// # Examples
    ///
    /// ```
    /// use seagrass_math::prelude::*;
    ///
    /// let points = [Vector3::new(1.0, 2.0, 3.0), Vector3::new(0.0, 0.0, 0.0)];
    /// let batch = Batch::<f64, Pga3>::from_points(&points);
    /// assert_eq!(batch.points_to_euclidean(), points);
    /// ```
    pub fn from_points(points: &[Vector3<T>]) -> Self {
        let len = points.len();
        let mut planes = Vec::with_capacity(4 * len);
        // Blade slots e123, e124, e134, e234 hold w, -z, y, -x.
        planes.extend(points.iter().map(|_| T::one()));
        planes.extend(points.iter().map(|p| -p.z));
        planes.extend(points.iter().map(|p| p.y));
        planes.extend(points.iter().map(|p| -p.x));
        Self {
            metric: Pga3::default(),
            blades: <u16 as BladeSet>::GRADE[3],
            len,
            planes,
        }
    }

    /// The Euclidean coordinates of every element read as a point.
    ///
    /// # Panics
    ///
    /// Panics if the batch does not store the four trivector blades of a point.
    pub fn points_to_euclidean(&self) -> Vec<Vector3<T>> {
        let plane = |index| self.plane(index).expect("the batch stores points");
        let (w, z, y, x) = (
            plane(BladeIndex::E123),
            plane(BladeIndex::E124),
            plane(BladeIndex::E134),
            plane(BladeIndex::E234),
        );
        (0..self.len)
            .map(|i| Vector3::new(-x[i] / w[i], y[i] / w[i], -z[i] / w[i]))
            .collect()
    }

    /// The rigid motions of homogeneous matrices, as motors that store every even blade.
    pub fn from_matrices(matrices: &[Matrix4<T>]) -> Self {
        let motors: Vec<Motor<T, Pga3>> = matrices.iter().map(Motor::from).collect();
        let mut batch = Self::zeros(Pga3::default(), <u16 as BladeSet>::EVEN, matrices.len());
        for (i, motor) in motors.iter().enumerate() {
            batch
                .set(i, motor.as_multivector())
                .expect("motors store even blades");
        }
        batch
    }

    /// The homogeneous matrix of every element read as a unit motor.
    pub fn motors_to_matrices(&self) -> Vec<Matrix4<T>> {
        (0..self.len)
            .map(|i| {
                let motor = Motor::new(self.get(i)).expect("the batch stores motors");
                motor.to_homogeneous_transformation_matrices().0
            })
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::metric::Dynamic;
    use crate::traits::{
        GeometricProduct, InnerProduct, LeftContraction, OuterProduct, RightContraction,
        SandwichProduct,
    };

    /// A single-value product, to compare a batch product against.
    type SingleProduct<B> = fn(
        &Multivector<f64, Dynamic<B>>,
        &Multivector<f64, Dynamic<B>>,
    ) -> Multivector<f64, Dynamic<B>>;

    /// Random values with random structures, from a small deterministic generator.
    fn values<B: BladeSet>(
        metric: Dynamic<B>,
        count: usize,
        seed: u64,
    ) -> Vec<Multivector<f64, Dynamic<B>>> {
        let mut state = seed;
        let mut next = move || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };
        (0..count)
            .map(|_| {
                let blades = B::from_u64(next() & next());
                let coefficients = blades
                    .iter()
                    .map(|_| (next() % 2001) as f64 / 500.0 - 2.0)
                    .collect();
                Multivector::from_raw(metric, blades, coefficients)
            })
            .collect()
    }

    fn check<B: BladeSet>(p: u8, q: u8, r: u8) {
        let metric = Dynamic::<B>::new(p, q, r).unwrap();
        // More elements than one block, so the block split is exercised.
        let (a, b) = (
            values(metric, 2 * BLOCK + 17, 1),
            values(metric, 2 * BLOCK + 17, 2),
        );
        let (batch_a, batch_b) = (
            Batch::from_values(metric, &a),
            Batch::from_values(metric, &b),
        );
        let padded = |v: &Multivector<f64, Dynamic<B>>, blades| {
            let mut out = Batch::zeros(metric, blades, 1);
            out.set(0, v).unwrap();
            out.get(0)
        };
        let single: [(Product, SingleProduct<B>); 5] = [
            (Product::Geometric, |a, b| a.geometric_product(b)),
            (Product::Outer, |a, b| a.outer_product(b)),
            (Product::Inner, |a, b| a.inner_product(b)),
            (Product::LeftContraction, |a, b| a.left_contraction(b)),
            (Product::RightContraction, |a, b| a.right_contraction(b)),
        ];
        let (pa, pb) = (
            |i: usize| padded(&a[i], batch_a.blades()),
            |i: usize| padded(&b[i], batch_b.blades()),
        );
        for (op, f) in single {
            let batch = batch_a.product(op, &batch_b);
            assert_eq!(
                batch_a.product_with(op, &batch_b, true),
                batch,
                "{op:?} in parallel"
            );
            for i in [0, 1, BLOCK, 2 * BLOCK + 16] {
                assert_eq!(
                    batch.get(i),
                    f(&pa(i), &pb(i)),
                    "{op:?} element {i} in ({p},{q},{r})"
                );
            }
        }
        let batch = batch_a.sandwich(&batch_b);
        for i in [0, BLOCK + 3] {
            let mut expected = pa(i).sandwich_product(&pb(i));
            let mut got = batch.get(i);
            expected.remove_zeros();
            got.remove_zeros();
            assert_eq!(got, expected, "sandwich element {i} in ({p},{q},{r})");
        }
    }

    #[test]
    fn elementwise_products_equal_single_value_products_bit_for_bit() {
        check::<u8>(3, 0, 0);
        check::<u16>(3, 0, 1);
        check::<u16>(1, 3, 0);
        check::<u32>(4, 1, 0);
        check::<u64>(3, 2, 1);
    }

    #[test]
    fn set_rejects_blades_outside_the_batch() {
        let metric = Dynamic::<u8>::new(3, 0, 0).unwrap();
        let mut batch = Batch::<f64, _>::zeros(metric, 0b11, 1);
        let value = Multivector::basis(metric, BladeIndex::E2);
        assert_eq!(
            batch.set(0, &value),
            Err(BladeError::NotStored {
                index: BladeIndex::E2
            })
        );
    }
}
