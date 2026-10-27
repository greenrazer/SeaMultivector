//! Sets of basis blades as bitsets, and the bit tricks that compute product structure.

use std::fmt::Debug;
use std::hash::Hash;
use std::ops::{
    BitAnd, BitAndAssign, BitOr, BitOrAssign, BitXor, BitXorAssign, Not, Shl, ShlAssign, Shr,
    ShrAssign,
};

use super::BladeIndex;
use super::tables;
use crate::typed::{Even, Kind, Odd};

/// A set of basis blades of one algebra, stored as a bitset with one bit per blade.
///
/// Implemented for `u8` (3D, 8 blades), `u16` (4D), `u32` (5D) and `u64` (6D). Bit i is set when
/// the set holds the blade with [`BladeIndex`] i. Besides the set operations, the trait carries
/// the sign and membership masks of the algebra and the bit tricks that find which blade pairs
/// of two multivectors land on each output blade.
///
/// # Examples
///
/// ```
/// use seagrass_math::{BladeIndex, BladeSet};
///
/// // {1, e12, e3} in 3D, then every blade multiplied by e2 at once.
/// let set = 0b0001_1001_u8;
/// assert_eq!(set.times_vector(2), 0b0100_0110); // {e1, e2, e23}
/// assert!(set.has(BladeIndex::E12));
/// ```
pub trait BladeSet:
    Copy
    + Eq
    + Ord
    + Hash
    + Debug
    + Default
    + Send
    + Sync
    + 'static
    + BitAnd<Output = Self>
    + BitOr<Output = Self>
    + BitXor<Output = Self>
    + Not<Output = Self>
    + Shl<u32, Output = Self>
    + Shr<u32, Output = Self>
    + BitAndAssign
    + BitOrAssign
    + BitXorAssign
    + ShlAssign<u32>
    + ShrAssign<u32>
{
    /// The number of basis vectors, 3 to 6.
    const DIM: u8;
    /// The number of basis blades, 2^`DIM`.
    const BLADES: u32;
    /// The empty set.
    const EMPTY: Self;
    /// The set holding only the scalar blade.
    const SCALAR: Self;
    /// Every blade.
    const ALL: Self;
    /// Blades of even grade.
    const EVEN: Self;
    /// Blades of odd grade.
    const ODD: Self;
    /// Blades that reversion negates.
    const REVERSE_SIGNS: Self;
    /// Blades that grade involution negates.
    const INVOLUTE_SIGNS: Self;
    /// Blades that Clifford conjugation negates.
    const CONJUGATE_SIGNS: Self;
    /// Output blades the right complement negates after bit reversal.
    const RIGHT_COMPLEMENT_SIGNS: Self;
    /// Output blades the left complement negates after bit reversal.
    const LEFT_COMPLEMENT_SIGNS: Self;
    /// Output blades the dual negates after bit reversal, for a metric without negative or null
    /// vectors.
    const DUAL_SIGNS: Self;
    /// `GRADE[g]` holds the blades of grade g.
    const GRADE: [Self; 7];
    /// `CONTAINS[k]` holds the blades that contain e(k); `CONTAINS[0]` is every blade.
    const CONTAINS: [Self; 7];
    /// `EXCLUDES[k]` holds the blades without e(k).
    const EXCLUDES: [Self; 7];
    /// `REORDER[k]` holds the blades that contain a vector numbered above k.
    const REORDER: [Self; 7];

    /// The grade parity the dual and the regressive product add: [`Even`] in even dimensions,
    /// [`Odd`] in odd ones.
    type Parity: Kind;

    /// The number of blades in the set.
    fn count_ones(self) -> u32;
    /// The index of the lowest blade in the set; `BLADES` when the set is empty.
    fn trailing_zeros(self) -> u32;
    /// Reverses the bits, which sends blade i to its complement 2^n - 1 - i.
    fn reverse_bits(self) -> Self;
    /// The low `BLADES` bits of a 6D mask.
    fn from_u64(bits: u64) -> Self;
    /// The bits as a `u64`.
    fn to_u64(self) -> u64;
    /// The blades with index below `index`.
    fn below(index: BladeIndex) -> Self;

    /// The set holding only `index`.
    #[inline]
    fn from_index(index: BladeIndex) -> Self {
        Self::SCALAR << u32::from(index.bits())
    }

    /// Whether the set holds `index`.
    #[inline]
    fn has(self, index: BladeIndex) -> bool {
        self & Self::from_index(index) != Self::EMPTY
    }

    /// Whether the set is empty.
    #[inline]
    fn is_empty(self) -> bool {
        self == Self::EMPTY
    }

    /// The number of blades in the set with index below `index`, which is the coefficient slot of
    /// `index` in a multivector with these blades.
    #[inline]
    fn count_below(self, index: BladeIndex) -> u32 {
        (self & Self::below(index)).count_ones()
    }

    /// The blades in ascending order.
    #[inline]
    fn iter(self) -> SetBlades<Self> {
        SetBlades(self)
    }

    /// Multiplies every blade in the set by basis vector e(`vector`), ignoring signs.
    ///
    /// Multiplying by e(k) flips bit k - 1 of every index, which swaps neighbouring groups of
    /// 2^(k-1) bits of the set: blades that contain e(k) move down and lose it, the rest move up
    /// and gain it. For example `{1, e12, e3} * e2` swaps groups of 2 bits:
    /// `0b0001_1001` becomes `0b0100_0110`, which is `{e1, e2, e23}`.
    ///
    /// # Panics
    ///
    /// Panics if `vector` is 0.
    #[inline]
    fn times_vector(self, vector: u8) -> Self {
        assert_ne!(vector, 0, "basis vectors are numbered from 1");
        let contains = Self::CONTAINS[vector as usize];
        let shift = 1 << (vector - 1);
        ((self & !contains) << shift) | ((self & contains) >> shift)
    }

    /// The sub-blades of `index`: the blades whose vectors all lie in `index`.
    ///
    /// Built by doubling: start from {1} and, for each vector at bit k of `index`, add a copy of
    /// the set shifted by 2^k, which adds that vector to every blade so far.
    #[inline]
    fn sub_blades(index: BladeIndex) -> Self {
        let mut out = Self::SCALAR;
        for position in index.positions() {
            out ^= out << (1_u32 << position);
        }
        out
    }

    /// The blades that contain none of the vectors of `index`.
    #[inline]
    fn without_any(index: BladeIndex) -> Self {
        let mut out = Self::ALL;
        for position in index.positions() {
            out &= Self::EXCLUDES[position as usize + 1];
        }
        out
    }

    /// The blades that contain every vector of `index`.
    #[inline]
    fn with_all(index: BladeIndex) -> Self {
        let mut out = Self::ALL;
        for position in index.positions() {
            out &= Self::CONTAINS[position as usize + 1];
        }
        out
    }

    /// For output blade `output`, the left blades x whose pair (x, x XOR `output`) has a
    /// negative reordering sign.
    #[inline]
    fn reorder_signs(output: BladeIndex) -> Self {
        let mut out = Self::REVERSE_SIGNS;
        for position in output.positions() {
            out ^= Self::REORDER[position as usize + 1];
        }
        out
    }

    /// For output blade `output`, the left blades x whose pair gets a factor -1 from the metric:
    /// every negative vector outside `output` lies in both blades of the pair when it lies in x,
    /// and squares to -1.
    #[inline]
    fn metric_signs(output: BladeIndex, negative_vectors: BladeIndex) -> Self {
        let squared = !output & negative_vectors;
        let mut out = Self::EMPTY;
        for position in squared.positions() {
            out ^= Self::CONTAINS[position as usize + 1];
        }
        out
    }

    /// The full sign mask of the pairs landing on `output`: reordering and metric signs.
    #[inline]
    fn product_signs(output: BladeIndex, negative_vectors: BladeIndex) -> Self {
        Self::reorder_signs(output) ^ Self::metric_signs(output, negative_vectors)
    }

    /// The sign mask for the dual (`inverse == false`) or undual in a metric without null vectors.
    #[inline]
    fn dual_signs(
        inverse: bool,
        pseudoscalar_inverse_negative: bool,
        negative_vectors: BladeIndex,
    ) -> Self {
        let mut out = Self::DUAL_SIGNS;
        if !inverse && pseudoscalar_inverse_negative {
            out = !out;
        }
        for position in negative_vectors.positions() {
            out ^= !Self::CONTAINS[position as usize + 1];
        }
        out
    }

    /// The blades of the product structure of a blade set with one blade `index` on its left,
    /// skipping pairs that share a null vector (one numbered above `non_null`).
    #[inline]
    fn geometric_structure_blade(self, index: BladeIndex, non_null: u8) -> Self {
        let mut out = self;
        for position in index.positions() {
            let vector = position + 1;
            out = out.times_vector(vector);
            if vector > non_null {
                out &= Self::CONTAINS[vector as usize];
            }
        }
        out
    }

    /// The blades the geometric product of two sets can produce.
    #[inline]
    fn geometric_structure(self, rhs: Self, non_null: u8) -> Self {
        let mut out = Self::EMPTY;
        for index in self.iter() {
            out |= rhs.geometric_structure_blade(index, non_null);
        }
        out
    }

    /// The blades `index ⌋ y` can produce for y in this set.
    #[inline]
    fn left_contraction_structure_blade(self, index: BladeIndex, non_null: u8) -> Self {
        if index.has_vector_above(non_null) {
            return Self::EMPTY;
        }
        let mut out = self >> u32::from(index.bits());
        for position in index.positions() {
            out &= !Self::CONTAINS[position as usize + 1];
        }
        out
    }

    /// The blades the left contraction of two sets can produce.
    #[inline]
    fn left_contraction_structure(self, rhs: Self, non_null: u8) -> Self {
        let mut out = Self::EMPTY;
        for index in self.iter() {
            out |= rhs.left_contraction_structure_blade(index, non_null);
        }
        out
    }

    /// The blades the right contraction of two sets can produce.
    #[inline]
    fn right_contraction_structure(self, rhs: Self, non_null: u8) -> Self {
        rhs.left_contraction_structure(self, non_null)
    }

    /// The blades the inner product (left plus right contraction) of two sets can produce.
    #[inline]
    fn inner_structure(self, rhs: Self, non_null: u8) -> Self {
        self.left_contraction_structure(rhs, non_null)
            | self.right_contraction_structure(rhs, non_null)
    }

    /// The blades `index ∧ y` can produce for y in this set.
    #[inline]
    fn outer_structure_blade(self, index: BladeIndex) -> Self {
        let mut out = self << u32::from(index.bits());
        for position in index.positions() {
            out &= Self::CONTAINS[position as usize + 1];
        }
        out
    }

    /// The blades the outer product of two sets can produce.
    #[inline]
    fn outer_structure(self, rhs: Self) -> Self {
        let mut out = Self::EMPTY;
        for index in self.iter() {
            out |= rhs.outer_structure_blade(index);
        }
        out
    }

    /// For a right operand with these blades, the left blades x whose partner `x XOR output` is
    /// in the set and whose pair does not vanish on a null vector: the pairs of a geometric
    /// product that land on `output`. `non_null` holds the vectors that are not null.
    ///
    /// Multiplying the set by `output` moves each partner onto x; ANDing with the sub-blades of
    /// `output | non_null` drops pairs that share a null vector.
    #[inline]
    fn geometric_pairs(self, output: BladeIndex, non_null: BladeIndex) -> Self {
        let mut out = self;
        for position in output.positions() {
            out = out.times_vector(position + 1);
        }
        out & Self::sub_blades(output | non_null)
    }

    /// The left blades x whose pair `x ⌋ (x + output)` is in the set and does not vanish.
    ///
    /// Every partner is `x + output`, so one right shift by `output` lines them all up.
    #[inline]
    fn left_contraction_pairs(self, output: BladeIndex, non_null: BladeIndex) -> Self {
        let shifted = (self >> u32::from(output.bits())) & Self::without_any(output);
        shifted & Self::sub_blades(output | non_null)
    }

    /// The left blades x whose pair `x ⌊ (x - output)` is in the set and does not vanish.
    #[inline]
    fn right_contraction_pairs(self, output: BladeIndex, non_null: BladeIndex) -> Self {
        let shifted = (self << u32::from(output.bits())) & Self::with_all(output);
        shifted & Self::sub_blades(output | non_null)
    }

    /// For the bit-reversed blades of a right operand, the left blades x whose pair
    /// `x ∧ (output - x)` is in the set.
    ///
    /// Reversal sends y to `MAX - y`, which turns the reflection `y = output - x` into a right
    /// shift by `MAX - output`.
    #[inline]
    fn outer_pairs(self, output: BladeIndex) -> Self {
        let max = Self::BLADES - 1;
        (self >> (max - u32::from(output.bits()))) & Self::sub_blades(output)
    }
}

/// The low bits of seven 6D masks as a narrower bitset type, optionally complemented.
macro_rules! truncate {
    ($ty:ty, $masks:expr, $complement:expr) => {{
        let mut out = [0 as $ty; 7];
        let mut k = 0;
        while k < 7 {
            out[k] = if $complement {
                !($masks[k] as $ty)
            } else {
                $masks[k] as $ty
            };
            k += 1;
        }
        out
    }};
}

macro_rules! impl_blade_set {
    ($($ty:ty => $dim:literal, $parity:ty;)*) => {$(
        impl BladeSet for $ty {
            const DIM: u8 = $dim;
            const BLADES: u32 = <$ty>::BITS;
            const EMPTY: Self = 0;
            const SCALAR: Self = 1;
            const ALL: Self = <$ty>::MAX;
            const EVEN: Self = tables::EVEN as $ty;
            const ODD: Self = !tables::EVEN as $ty;
            const REVERSE_SIGNS: Self = tables::REVERSE_SIGNS as $ty;
            const INVOLUTE_SIGNS: Self = !tables::EVEN as $ty;
            const CONJUGATE_SIGNS: Self = tables::CONJUGATE_SIGNS as $ty;
            const RIGHT_COMPLEMENT_SIGNS: Self = tables::right_complement_signs($dim) as $ty;
            const LEFT_COMPLEMENT_SIGNS: Self = tables::left_complement_signs($dim) as $ty;
            const DUAL_SIGNS: Self = tables::dual_signs($dim) as $ty;
            const GRADE: [Self; 7] = truncate!($ty, tables::GRADE, false);
            const CONTAINS: [Self; 7] = truncate!($ty, tables::CONTAINS, false);
            const EXCLUDES: [Self; 7] = truncate!($ty, tables::CONTAINS, true);
            const REORDER: [Self; 7] = truncate!($ty, tables::REORDER, false);

            type Parity = $parity;

            #[inline]
            fn count_ones(self) -> u32 {
                <$ty>::count_ones(self)
            }

            #[inline]
            fn trailing_zeros(self) -> u32 {
                <$ty>::trailing_zeros(self)
            }

            #[inline]
            fn reverse_bits(self) -> Self {
                <$ty>::reverse_bits(self)
            }

            #[inline]
            fn from_u64(bits: u64) -> Self {
                bits as $ty
            }

            #[inline]
            fn to_u64(self) -> u64 {
                u64::from(self)
            }

            #[inline]
            fn below(index: BladeIndex) -> Self {
                (1 as $ty).checked_shl(u32::from(index.bits())).unwrap_or(0).wrapping_sub(1)
            }
        }

    )*};
}

impl_blade_set! {
    u8 => 3, Odd;
    u16 => 4, Even;
    u32 => 5, Odd;
    u64 => 6, Even;
}

/// Iterator over the blades of a [`BladeSet`], lowest index first.
///
/// # Examples
///
/// ```
/// use seagrass_math::{BladeIndex, BladeSet};
///
/// let blades: Vec<BladeIndex> = 0b0100_1001_u8.iter().collect();
/// assert_eq!(blades, [BladeIndex::SCALAR, BladeIndex::E12, BladeIndex::E23]);
/// ```
#[derive(Debug, Clone)]
pub struct SetBlades<B>(B);

impl<B: BladeSet> Iterator for SetBlades<B> {
    type Item = BladeIndex;

    #[inline]
    fn next(&mut self) -> Option<BladeIndex> {
        if self.0.is_empty() {
            return None;
        }
        let index = self.0.trailing_zeros() as u8;
        self.0 ^= B::SCALAR << u32::from(index);
        Some(BladeIndex::new(index))
    }

    #[inline]
    fn size_hint(&self) -> (usize, Option<usize>) {
        let n = self.0.count_ones() as usize;
        (n, Some(n))
    }
}

impl<B: BladeSet> ExactSizeIterator for SetBlades<B> {}
