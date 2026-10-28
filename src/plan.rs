//! Product plans: the structure of a product, worked out once per pair of operand blade sets.
//!
//! The bitset tricks in [`BladeSet`] find, for each output blade, which blade pairs land on it and
//! with which sign. A plan records that result as a list of multiply-adds between coefficient
//! slots, so a product with the same operand blade sets skips the bit work, the popcounts and
//! the mask building, and runs one tight loop. Plans are cached per thread.

use std::cell::RefCell;
use std::collections::HashMap;
use std::hash::{BuildHasherDefault, Hash, Hasher};
use std::rc::Rc;

use smallvec::SmallVec;

use crate::blade::{BladeIndex, BladeSet};
use crate::coefficient::Coefficient;
use crate::metric::{Metric, ProductConstants};

/// Coefficient storage: up to 16 coefficients inline, so values up to 4D never allocate.
pub(crate) type Coefficients<T> = SmallVec<[T; 16]>;

/// The products that have plans. The regressive product is the outer product of duals, which
/// are cheap sign flips, and the scalar product touches at most one pair per blade.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub(crate) enum Op {
    Geometric,
    Outer,
    Inner,
    LeftContraction,
    RightContraction,
    Sandwich,
}

/// What a plan depends on: the operation, the signature and both operand blade sets.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct PlanKey {
    op: Op,
    signature: [u8; 3],
    lhs: u64,
    rhs: u64,
}

impl PlanKey {
    fn new<M: Metric>(op: Op, metric: &M, lhs: M::Bits, rhs: M::Bits) -> Self {
        Self {
            op,
            signature: [metric.positive(), metric.negative(), metric.zero()],
            lhs: lhs.to_u64(),
            rhs: rhs.to_u64(),
        }
    }
}

impl Hash for PlanKey {
    fn hash<H: Hasher>(&self, state: &mut H) {
        let [p, q, r] = self.signature;
        state.write_u64(
            (self.op as u64) << 24 | u64::from(p) << 16 | u64::from(q) << 8 | u64::from(r),
        );
        state.write_u64(self.lhs);
        state.write_u64(self.rhs);
    }
}

/// The Fx hash from rustc: one rotate, xor and multiply per word. Plan keys are a few machine
/// words chosen by the program, not by an adversary, so a fast non-cryptographic hash is enough.
#[derive(Default)]
struct FxHasher(u64);

impl Hasher for FxHasher {
    fn write(&mut self, bytes: &[u8]) {
        for &byte in bytes {
            self.write_u64(u64::from(byte));
        }
    }

    #[inline]
    fn write_u64(&mut self, word: u64) {
        self.0 = (self.0.rotate_left(5) ^ word).wrapping_mul(0x517c_c1b7_2722_0a95);
    }

    #[inline]
    fn finish(&self) -> u64 {
        self.0
    }
}

/// Plans kept per thread. A program that keeps meeting new blade-set pairs, such as one fed
/// random structures, is bounded by clearing the cache when it is full.
const CACHE_CAPACITY: usize = 4096;

/// A per-thread plan cache: a hash map, plus the last plan used for each operation, so a loop
/// that repeats one product compares four words instead of hashing.
struct PlanCache<P> {
    plans: HashMap<PlanKey, Rc<P>, BuildHasherDefault<FxHasher>>,
    recent: [Option<(PlanKey, Rc<P>)>; 6],
}

impl<P> Default for PlanCache<P> {
    fn default() -> Self {
        Self {
            plans: HashMap::default(),
            recent: Default::default(),
        }
    }
}

impl<P> PlanCache<P> {
    /// The plan for `key`, built by `build` the first time.
    #[inline]
    fn get(&mut self, key: PlanKey, build: impl FnOnce() -> P) -> &P {
        let slot = key.op as usize;
        if !matches!(&self.recent[slot], Some((recent, _)) if *recent == key) {
            if self.plans.len() >= CACHE_CAPACITY && !self.plans.contains_key(&key) {
                self.plans.clear();
            }
            let plan = Rc::clone(self.plans.entry(key).or_insert_with(|| Rc::new(build())));
            self.recent[slot] = Some((key, plan));
        }
        match &self.recent[slot] {
            Some((_, plan)) => plan,
            None => unreachable!("the slot was just filled"),
        }
    }
}

thread_local! {
    static PRODUCT_PLANS: RefCell<PlanCache<ProductPlan>> = RefCell::default();
    static SANDWICH_PLANS: RefCell<PlanCache<SandwichPlan>> = RefCell::default();
}

/// Runs `run` with the cached plan of `op` for these operand blade sets, building it first if
/// needed.
#[inline]
pub(crate) fn with_product_plan<M: Metric, R>(
    op: Op,
    metric: &M,
    lhs: M::Bits,
    rhs: M::Bits,
    run: impl FnOnce(&ProductPlan) -> R,
) -> R {
    let key = PlanKey::new(op, metric, lhs, rhs);
    PRODUCT_PLANS.with(|cache| {
        let mut cache = cache.borrow_mut();
        run(cache.get(key, || {
            ProductPlan::build(op, ProductConstants::of(metric), lhs, rhs)
        }))
    })
}

/// Runs `run` with the cached sandwich plan for a versor and an operand with these blade sets.
#[inline]
pub(crate) fn with_sandwich_plan<M: Metric, R>(
    metric: &M,
    versor: M::Bits,
    operand: M::Bits,
    run: impl FnOnce(&SandwichPlan) -> R,
) -> R {
    let key = PlanKey::new(Op::Sandwich, metric, versor, operand);
    SANDWICH_PLANS.with(|cache| {
        let mut cache = cache.borrow_mut();
        run(cache.get(key, || {
            SandwichPlan::build(ProductConstants::of(metric), versor, operand)
        }))
    })
}

/// One multiply-add of a bilinear product: `± lhs[lhs] * rhs[rhs]`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct Entry {
    pub(crate) lhs: u8,
    pub(crate) rhs: u8,
    pub(crate) negative: bool,
}

/// A bilinear product with fixed operand blade sets, as multiply-adds between coefficient slots.
///
/// Every output coefficient is a fixed sum that no other output touches: the entries of output
/// slot k are `entries[offsets[k]..offsets[k + 1]]`.
#[derive(Debug, Clone)]
pub(crate) struct ProductPlan {
    output: u64,
    offsets: Box<[u16]>,
    entries: Box<[Entry]>,
    #[cfg(feature = "simd")]
    by_sign: BySign,
}

impl ProductPlan {
    /// Works out the product structure with the bit tricks of [`BladeSet`]. Each output's
    /// entries are listed in the order the per-call loop visits its pairs, so a plan sums in the
    /// same order and gives the same bits.
    pub(crate) fn build<B: BladeSet>(op: Op, constants: ProductConstants, lhs: B, rhs: B) -> Self {
        let output = match op {
            Op::Geometric => lhs.geometric_structure(rhs, constants.non_null_count),
            Op::Outer => lhs.outer_structure(rhs),
            Op::Inner => lhs.inner_structure(rhs, constants.non_null_count),
            Op::LeftContraction => lhs.left_contraction_structure(rhs, constants.non_null_count),
            Op::RightContraction => lhs.right_contraction_structure(rhs, constants.non_null_count),
            Op::Sandwich => unreachable!("sandwiches have their own plan type"),
        };
        let reversed_rhs = rhs.reverse_bits();

        let mut offsets = Vec::with_capacity(output.count_ones() as usize + 1);
        let mut entries = Vec::new();
        offsets.push(0);
        for o in output.iter() {
            let mut push_pairs = |pairs: B, signs: B| {
                for x in pairs.iter() {
                    entries.push(Entry {
                        lhs: lhs.count_below(x) as u8,
                        rhs: rhs.count_below(x ^ o) as u8,
                        negative: signs.has(x),
                    });
                }
            };
            let signs = B::product_signs(o, constants.negative);
            match op {
                Op::Geometric => {
                    push_pairs(lhs & rhs.geometric_pairs(o, constants.non_null), signs)
                }
                Op::Outer => push_pairs(lhs & reversed_rhs.outer_pairs(o), B::reorder_signs(o)),
                Op::LeftContraction => {
                    push_pairs(
                        lhs & rhs.left_contraction_pairs(o, constants.non_null),
                        signs,
                    );
                }
                Op::RightContraction => {
                    push_pairs(
                        lhs & rhs.right_contraction_pairs(o, constants.non_null),
                        signs,
                    );
                }
                Op::Inner => {
                    push_pairs(
                        lhs & rhs.left_contraction_pairs(o, constants.non_null),
                        signs,
                    );
                    // Scalar pairs are counted once, by the left contraction.
                    if o != BladeIndex::SCALAR {
                        push_pairs(
                            lhs & rhs.right_contraction_pairs(o, constants.non_null),
                            signs,
                        );
                    }
                }
                Op::Sandwich => unreachable!("sandwiches have their own plan type"),
            }
            offsets.push(entries.len() as u16);
        }
        Self {
            output: output.to_u64(),
            #[cfg(feature = "simd")]
            by_sign: BySign::new(&offsets, &entries),
            offsets: offsets.into(),
            entries: entries.into(),
        }
    }

    /// The output blade set.
    #[inline]
    pub(crate) fn output<B: BladeSet>(&self) -> B {
        B::from_u64(self.output)
    }

    /// The entries of each output slot, in slot order.
    pub(crate) fn outputs(&self) -> impl ExactSizeIterator<Item = &[Entry]> + '_ {
        self.offsets
            .windows(2)
            .map(|range| &self.entries[usize::from(range[0])..usize::from(range[1])])
    }

    /// The same entries with each output's negative ones first.
    #[cfg(feature = "simd")]
    pub(crate) fn by_sign(&self) -> &BySign {
        &self.by_sign
    }

    /// The output coefficients for operands with the plan's blade sets.
    #[inline]
    pub(crate) fn run<T: Coefficient>(&self, lhs: &[T], rhs: &[T]) -> Coefficients<T> {
        let mut out = Coefficients::from_elem(T::zero(), self.offsets.len() - 1);
        for (slot, range) in out.iter_mut().zip(self.offsets.windows(2)) {
            let mut sum = T::zero();
            for entry in &self.entries[usize::from(range[0])..usize::from(range[1])] {
                let product = lhs[usize::from(entry.lhs)] * rhs[usize::from(entry.rhs)];
                if entry.negative {
                    sum -= product;
                } else {
                    sum += product;
                }
            }
            *slot = sum;
        }
        out
    }
}

/// One term of a sandwich: `coefficient * m[first] * x[operand] * m[second]`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct Term {
    pub(crate) first: u8,
    pub(crate) operand: u8,
    pub(crate) second: u8,
    pub(crate) coefficient: i8,
}

/// The sandwich `M X reverse(M)` with fixed blade sets, flattened to trilinear terms.
///
/// Built from the two geometric products of the per-call sandwich, so the reverse of M is still
/// folded into the second product's sign mask. Expanding both products gives terms
/// `± m[i] x[j] m[k]`; since coefficients commute, the terms for (i, j, k) and (k, j, i) are
/// merged, and terms that cancel are dropped, along with output blades whose terms all cancel.
/// For a PGA3D motor applied to a point that leaves 40 terms where the two products have 68
/// pairs. Offsets are `u32` because two dense 6D values give up to 64³ triples, more than a
/// `u16` can count.
#[derive(Debug, Clone)]
pub(crate) struct SandwichPlan {
    output: u64,
    offsets: Box<[u32]>,
    terms: Box<[Term]>,
}

impl SandwichPlan {
    pub(crate) fn build<B: BladeSet>(constants: ProductConstants, versor: B, operand: B) -> Self {
        let first = ProductPlan::build(Op::Geometric, constants, versor, operand);
        let first_blades: B = first.output();
        let output = first_blades.geometric_structure(versor, constants.non_null_count);

        // (output blade, i, j, k) with i <= k -> summed integer coefficient, in sorted order.
        let mut sums = std::collections::BTreeMap::<(u8, u8, u8, u8), i32>::new();
        for o in output.iter() {
            let pairs = first_blades & versor.geometric_pairs(o, constants.non_null);
            let signs = B::product_signs(o, constants.negative)
                ^ B::REVERSE_SIGNS.geometric_pairs(o, constants.non_null);
            for t in pairs.iter() {
                let k = versor.count_below(t ^ o) as u8;
                let outer_sign = if signs.has(t) { -1 } else { 1 };
                let slot = first_blades.count_below(t) as usize;
                let range = usize::from(first.offsets[slot])..usize::from(first.offsets[slot + 1]);
                for entry in &first.entries[range] {
                    let inner_sign = if entry.negative { -1 } else { 1 };
                    let (lo, hi) = (entry.lhs.min(k), entry.lhs.max(k));
                    *sums.entry((o.bits(), lo, entry.rhs, hi)).or_default() +=
                        outer_sign * inner_sign;
                }
            }
        }

        let mut output_blades = B::EMPTY;
        for (&(o, ..), &c) in &sums {
            if c != 0 {
                output_blades |= B::from_index(BladeIndex::new(o));
            }
        }
        let mut offsets = vec![0_u32];
        let mut terms = Vec::new();
        let mut current = None;
        for ((o, first, operand, second), coefficient) in sums {
            if coefficient == 0 {
                continue;
            }
            if current.is_some_and(|c| c != o) {
                offsets.push(terms.len() as u32);
            }
            current = Some(o);
            // Each ordered (i, j, k) lands on exactly one output blade, so a merged term sums at
            // most two of them: the coefficient is -2, -1, 1 or 2.
            let coefficient = i8::try_from(coefficient).expect("a merged term sums two triples");
            terms.push(Term {
                first,
                operand,
                second,
                coefficient,
            });
        }
        if current.is_some() {
            offsets.push(terms.len() as u32);
        }
        Self {
            output: output_blades.to_u64(),
            offsets: offsets.into(),
            terms: terms.into(),
        }
    }

    /// The output blade set: the blades with at least one term that does not cancel.
    #[inline]
    pub(crate) fn output<B: BladeSet>(&self) -> B {
        B::from_u64(self.output)
    }

    /// The terms of each output slot, in slot order.
    pub(crate) fn outputs(&self) -> impl ExactSizeIterator<Item = &[Term]> + '_ {
        self.offsets
            .windows(2)
            .map(|range| &self.terms[range[0] as usize..range[1] as usize])
    }

    /// The output coefficients for a versor `m` and operand `x` with the plan's blade sets.
    #[inline]
    pub(crate) fn run<T: Coefficient>(&self, m: &[T], x: &[T]) -> Coefficients<T> {
        let mut out = Coefficients::from_elem(T::zero(), self.offsets.len().saturating_sub(1));
        for (slot, range) in out.iter_mut().zip(self.offsets.windows(2)) {
            let mut sum = T::zero();
            for term in &self.terms[range[0] as usize..range[1] as usize] {
                let product = m[usize::from(term.first)]
                    * x[usize::from(term.operand)]
                    * m[usize::from(term.second)];
                sum += match term.coefficient {
                    1 => product,
                    -1 => -product,
                    c => T::from(c).expect("small integers are representable") * product,
                };
            }
            *slot = sum;
        }
        out
    }
}

/// A product plan's entries with each output's negative entries first, keeping plan order within
/// a sign, for kernels that sum each sign without branching on it: output k's negative entries
/// are `entries[bounds[2 k]..bounds[2 k + 1]]` and its positive ones
/// `entries[bounds[2 k + 1]..bounds[2 k + 2]]`.
#[cfg(feature = "simd")]
#[derive(Debug, Clone)]
pub(crate) struct BySign {
    bounds: Box<[u16]>,
    entries: Box<[Entry]>,
}

#[cfg(feature = "simd")]
impl BySign {
    fn new(offsets: &[u16], entries: &[Entry]) -> Self {
        let mut bounds = Vec::with_capacity(2 * offsets.len() - 1);
        let mut grouped = Vec::with_capacity(entries.len());
        for range in offsets.windows(2) {
            let output = &entries[usize::from(range[0])..usize::from(range[1])];
            bounds.push(grouped.len() as u16);
            grouped.extend(output.iter().filter(|entry| entry.negative));
            bounds.push(grouped.len() as u16);
            grouped.extend(output.iter().filter(|entry| !entry.negative));
        }
        bounds.push(grouped.len() as u16);
        Self {
            bounds: bounds.into(),
            entries: grouped.into(),
        }
    }

    /// The `[negative, positive]` entries of each output slot, in slot order.
    pub(crate) fn outputs(&self) -> impl ExactSizeIterator<Item = [&[Entry]; 2]> + '_ {
        self.bounds.windows(3).step_by(2).map(|bounds| {
            let [start, middle, end] = [bounds[0], bounds[1], bounds[2]].map(usize::from);
            [&self.entries[start..middle], &self.entries[middle..end]]
        })
    }
}
