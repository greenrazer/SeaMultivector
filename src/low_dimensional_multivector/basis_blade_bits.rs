use std::ops::{Deref, DerefMut, AddAssign, BitOrAssign, Shl, Add, BitOr, BitAndAssign, BitAnd, Not, ShlAssign, ShrAssign, BitXor, Shr, BitXorAssign};
use std::str::FromStr;

use funty::Unsigned;

use crate::errors::ParseMultivectorError;
use crate::multivector_product_traits::{GeometricProduct, LeftContraction, RightContraction, InnerProduct, OuterProduct};

use super::basis_blade_index::SeaBasisBladeIndexLow;


pub trait SeaBasisBladeBitsInner: Unsigned{ 
	const DIMENSIONS: u8;
	const NUM_BASIS_BLADES: u8 = Self::BITS as u8;
	
	/// All bits of an even grade,
	/// or all bits with an index with an even number of bits set
	/// or even binary weight.
	/// 
	/// e.g 10010110
	/// because for each index
	/// 0 -> 0 bits set -> even
	/// 1 -> 1 bits set -> odd
	/// 2 -> 1 bits set -> odd
	/// 3 -> 2 bits set -> even
	/// 4 -> 1 bits set -> odd
	/// 5 -> 2 bits set -> odd
	/// 6 -> 2 bits set -> odd
	/// 7 -> 3 bits set -> odd
	const EVEN_GRADE_BITS: Self;

	/// All bits of an odd grade,
	/// or all bits with an index with an odd number of bits set
	/// or odd binary weight.
	const ODD_GRADE_BITS: Self;

	/// All indexes that would be negitive after
	/// being cannonically reordered with themselves,
	/// or all indexes that would be negitive after reversing all
	/// it's basis vectors.
	/// 
	/// if 1/2 * binary weight of index * (binary weight of index-1) is even -> 0
	/// if odd -> 1
	const REVERSION_NEGATION_BITS: Self;

	/// xor between reversion_negation_bits and odd_grade_bits
	const CONJUGATION_NEGATION_BITS: Self;

	/// right compliment negation bits
	/// such that x * compliment(x) = I(pseudo scalar)
	/// complement = reverse(x)*these_bits
	const RIGHT_COMPLIMENT_NEGATION_BITS: Self;

	/// left compliment negation bits
	/// such that compliment(x) * x  = I(pseudo scalar)
	/// complement = reverse(x)*these_bits
	const LEFT_COMPLIMENT_NEGATION_BITS: Self;

	/// negation bits corresponding to the sign of
	/// sign(x << I(pseudo scalar))
	/// for positive-only metrics
	const DUAL_NEGATION_BITS: Self;
	
	/// for each grade, the bits in that grade
	/// for example for grade 0 the only bit set is always in the 0's place
	/// for grade 1 every index with a binary weight of 1 is 1
	/// for grade 2 every index with a binary weight of 2 is 1
	/// ...
	const GRADE_BITS: [Self; 7];
	
	/// for each basis_vector(1->0, e1->1, e2->2, e3->3)
	/// 1 at each index that contains that basis vector
	/// e.g. 
	/// basis vector 1
	/// 01010101
	/// basis vector 2
	/// 00110011
	const CONTAINS_BASIS_VEC_BITS: [Self; 7];
	
	/// opposite of contains basis vector bits
	const EXCLUDES_BASIS_VEC_BITS: [Self; 7];

	/// a shorter name for CONTAINS_BASIS_VEC_BITS
	const K: [Self; 7] = Self::CONTAINS_BASIS_VEC_BITS;

	/// for each basis vector in an output index after some product
	/// this negates the reversion grade bits to get the current
	/// cannonical reordering sign bits
	/// for example for each product output index
	/// the cannonical reordering sign mask is
	/// 1 -> Reversion bits
	/// e1 -> Reversion bits ^ output reversion negation bits[1]
	/// e2 -> Reversion bits ^ output reversion negation bits[2]
	/// e12 -> Reversion bits ^ output reversion negation bits[1] ^ output reversion negation bits[2]
	const OUTPUT_REVERSION_NEGATION_BITS: [Self; 7];

	/// sets all least signifigant bits up to n
	/// e.g. 5
	/// ...0011111
	/// 3
	/// ...0000111
	fn set_lsb_up_to_n(n: u8) -> Self {
		let one = Self::ONE;
		one.checked_shl(n as u32).unwrap_or(Self::ZERO).overflowing_sub(one).0
	}

	/// swaps adjacent bit groups with length 2^n
	/// e.g. 2 1111011001111001
	/// 2^2 = 4
	/// 1111 0110 0111 1001
	/// 0110 1111 1001 0111
	/// 0110111110010111
	fn swap_adjacent_pow2_n_bit_groups(&mut self, n:u8) {
		assert_ne!(n, 0);
		let bitmask = Self::CONTAINS_BASIS_VEC_BITS[n as usize];
		let shift_amnt = Self::ONE << (n-1);
		let msb = *self & bitmask;
		*self &= !bitmask;
		*self <<= shift_amnt;
		*self |= msb >> shift_amnt;
	}
}
impl SeaBasisBladeBitsInner for u8 {
	const DIMENSIONS: u8 = 3;
	// 0110 1001
	const EVEN_GRADE_BITS: Self = 0x69;
	// 1001 0110
	const ODD_GRADE_BITS: Self = !Self::EVEN_GRADE_BITS;
	// 1110 1000
	const REVERSION_NEGATION_BITS: Self = 0xE8;
	// 0111 1110
	const CONJUGATION_NEGATION_BITS: Self = 0x7E;
	// 0010 0100
	const RIGHT_COMPLIMENT_NEGATION_BITS: Self = 0x24;
	// 0010 0100
	const LEFT_COMPLIMENT_NEGATION_BITS: Self = 0x24;
	// 0011 0011
	const DUAL_NEGATION_BITS: Self = 0x33;

	const GRADE_BITS: [Self;7] = [
		0x01,
		0x16,
		0x68,
		0x80,
		0,
		0,
		0
	];
	const CONTAINS_BASIS_VEC_BITS: [Self; 7] = [
		0xFF,
		0xAA,
		0xCC,
		0xF0,
		0,
		0,
		0
	];
	const EXCLUDES_BASIS_VEC_BITS: [Self; 7] = [
		!Self::CONTAINS_BASIS_VEC_BITS[0],
		!Self::CONTAINS_BASIS_VEC_BITS[1],
		!Self::CONTAINS_BASIS_VEC_BITS[2],
		!Self::CONTAINS_BASIS_VEC_BITS[3],
		!Self::CONTAINS_BASIS_VEC_BITS[4],
		!Self::CONTAINS_BASIS_VEC_BITS[5],
		!Self::CONTAINS_BASIS_VEC_BITS[6],
	];
	const OUTPUT_REVERSION_NEGATION_BITS: [Self; 7] = [
		Self::K[1] ^ Self::K[2] ^ Self::K[3] ^ Self::K[4] ^ Self::K[5] ^ Self::K[6],
					 Self::K[2] ^ Self::K[3] ^ Self::K[4] ^ Self::K[5] ^ Self::K[6],
					 			  Self::K[3] ^ Self::K[4] ^ Self::K[5] ^ Self::K[6],
					 			  			   Self::K[4] ^ Self::K[5] ^ Self::K[6],
					 			  			                Self::K[5] ^ Self::K[6],
					 			  			                			 Self::K[6],
					 			  			                			 0,
	];
}
impl SeaBasisBladeBitsInner for u16 {
	const DIMENSIONS: u8 = 4;
	// 1001 0110  0110 1001
	const EVEN_GRADE_BITS: Self = 0x96_69;
	// 0110 1001  1001 0110
	const ODD_GRADE_BITS: Self = !Self::EVEN_GRADE_BITS;
	// 0111 1110  1110 1000
	const REVERSION_NEGATION_BITS: Self = 0x7E_E8;
	// 0001 0111  0111 1110
	const CONJUGATION_NEGATION_BITS: Self = 0x17_7E;
	// 0100 1101 0010 0100
	const RIGHT_COMPLIMENT_NEGATION_BITS: Self = 0x4D_24;
	// 0010 0100 1011 0010
	const LEFT_COMPLIMENT_NEGATION_BITS: Self = 0x24_B2;
	// 0011 0011 1100 1100
	const DUAL_NEGATION_BITS: Self = 0x33_CC;

	#[allow(clippy::mistyped_literal_suffixes)]
	const GRADE_BITS: [Self;7] = [
		0x00_01,
		0x01_16,
		0x16_68,
		0x68_80,
		0x80_00,
		0,
		0
	];
	const CONTAINS_BASIS_VEC_BITS: [Self; 7] = [
		0xFF_FF,
		0xAA_AA,
		0xCC_CC,
		0xF0_F0,
		0xFF_00,
		0,
		0
	];
	const EXCLUDES_BASIS_VEC_BITS: [Self; 7] = [
		!Self::CONTAINS_BASIS_VEC_BITS[0],
		!Self::CONTAINS_BASIS_VEC_BITS[1],
		!Self::CONTAINS_BASIS_VEC_BITS[2],
		!Self::CONTAINS_BASIS_VEC_BITS[3],
		!Self::CONTAINS_BASIS_VEC_BITS[4],
		!Self::CONTAINS_BASIS_VEC_BITS[5],
		!Self::CONTAINS_BASIS_VEC_BITS[6],
	];
	const OUTPUT_REVERSION_NEGATION_BITS: [Self; 7] = [
		Self::K[1] ^ Self::K[2] ^ Self::K[3] ^ Self::K[4] ^ Self::K[5] ^ Self::K[6],
					 Self::K[2] ^ Self::K[3] ^ Self::K[4] ^ Self::K[5] ^ Self::K[6],
					 			  Self::K[3] ^ Self::K[4] ^ Self::K[5] ^ Self::K[6],
					 			  			   Self::K[4] ^ Self::K[5] ^ Self::K[6],
					 			  			                Self::K[5] ^ Self::K[6],
					 			  			                			 Self::K[6],
					 			  			                			 0,
	];
}
impl SeaBasisBladeBitsInner for u32 {
	const DIMENSIONS: u8 = 5;
	// 0110 1001  1001 0110  1001 0110  0110 1001
	const EVEN_GRADE_BITS: Self = 0x69969669;
	// 1001 0110  0110 1001  0110 1001  1001 0110
	const ODD_GRADE_BITS: Self = !Self::EVEN_GRADE_BITS;
	// 0001 0111  0111 1110  0111 1110  1110 1000
	const REVERSION_NEGATION_BITS: Self = 0x177E7EE8;
	// 1000 0001  0001 0111  0001 0111  0111 1110
	const CONJUGATION_NEGATION_BITS: Self = 0x8117177E;
	// 0010 0100 1011 0010 0100 1101 0010 0100
	const RIGHT_COMPLIMENT_NEGATION_BITS: Self = 0x24_B2_4D_24;
	// 0010 0100 1011 0010 0100 1101 0010 0100
	const LEFT_COMPLIMENT_NEGATION_BITS: Self = 0x24_B2_4D_24;
	// 0011 0011 1100 1100 0011 0011 1100 1100
	const DUAL_NEGATION_BITS: Self = 0x33_CC_33_CC;

	#[allow(clippy::mistyped_literal_suffixes)]
	const GRADE_BITS: [Self; 7] = [
		0x00_00_00_01,
		0x00_01_01_16,
		0x01_16_16_68,
		0x16_68_68_80,
		0x68_80_80_00,
		0x80_00_00_00,
		0
	];
	const CONTAINS_BASIS_VEC_BITS: [Self; 7] = [
		0xFF_FF_FF_FF,
		0xAA_AA_AA_AA,
		0xCC_CC_CC_CC,
		0xF0_F0_F0_F0,
		0xFF_00_FF_00,
		0xFF_FF_00_00,
		0
	];
	const EXCLUDES_BASIS_VEC_BITS: [Self; 7] = [
		!Self::CONTAINS_BASIS_VEC_BITS[0],
		!Self::CONTAINS_BASIS_VEC_BITS[1],
		!Self::CONTAINS_BASIS_VEC_BITS[2],
		!Self::CONTAINS_BASIS_VEC_BITS[3],
		!Self::CONTAINS_BASIS_VEC_BITS[4],
		!Self::CONTAINS_BASIS_VEC_BITS[5],
		!Self::CONTAINS_BASIS_VEC_BITS[6],
	];
	const OUTPUT_REVERSION_NEGATION_BITS: [Self; 7] = [
		Self::K[1] ^ Self::K[2] ^ Self::K[3] ^ Self::K[4] ^ Self::K[5] ^ Self::K[6],
					 Self::K[2] ^ Self::K[3] ^ Self::K[4] ^ Self::K[5] ^ Self::K[6],
					 			  Self::K[3] ^ Self::K[4] ^ Self::K[5] ^ Self::K[6],
					 			  			   Self::K[4] ^ Self::K[5] ^ Self::K[6],
					 			  			                Self::K[5] ^ Self::K[6],
					 			  			                			 Self::K[6],
					 			  			                			 0,
	];
}
impl SeaBasisBladeBitsInner for u64 {
	const DIMENSIONS: u8 = 6;
	// 1001 0110  0110 1001  0110 1001  1001 0110  0110 1001  1001 0110  1001 0110  0110 1001
	const EVEN_GRADE_BITS: Self = 0x96_69_69_96_69_96_96_69;
	// 0110 1001  1001 0110  1001 0110  0110 1001  1001 0110  0110 1001  0110 1001  1001 0110
	const ODD_GRADE_BITS: Self = !Self::EVEN_GRADE_BITS;
	// 1000 0001  0001 0111  0001 0111  0111 1110  0001 0111  0111 1110  0111 1110  1110 1000
	const REVERSION_NEGATION_BITS: Self = 0x81_17_17_7E_17_7E_7E_E8;
	// 1110 1000  1000 0001  1000 0001  0001 0111  1000 0001  0001 0111  0001 0111  0111 1110
	const CONJUGATION_NEGATION_BITS: Self = 0xE8_81_81_17_81_17_17_7E;
	// 0100 1011 0010 0100 1101 1011 0100 1101 0010 0100 1011 0010 0100 1101 0010 0100
	const RIGHT_COMPLIMENT_NEGATION_BITS: Self = 0x4D_24_DB_4D_24_B2_4D_24;
	// 0010 0100 1011 0010 0100 1101 0010 0100 1011 0010 1101 1011 0010 0100 1101 0010
	const LEFT_COMPLIMENT_NEGATION_BITS: Self = 0x24_D2_4D_24_B2_DB_24_D2;
	// 0011 0011 1100 1100 0011 0011 1100 1100 1100 1100 0011 0011 1100 1100 0011 0011
	const DUAL_NEGATION_BITS: Self = 0x33_CC_33_CC_CC_33_CC_33;

	const GRADE_BITS: [Self; 7] = [
		0x00_00_00_00_00_00_00_01,
		0x00_00_00_01_00_01_01_16,
		0x00_01_01_16_01_16_16_68,
		0x01_16_16_68_16_68_68_80,
		0x16_68_68_80_68_80_80_00,
		0x68_80_80_00_80_00_00_00,
		0x80_00_00_00_00_00_00_00
	];

	const CONTAINS_BASIS_VEC_BITS: [Self; 7] = [
		0xFF_FF_FF_FF_FF_FF_FF_FF,
		0xAA_AA_AA_AA_AA_AA_AA_AA,
		0xCC_CC_CC_CC_CC_CC_CC_CC,
		0xF0_F0_F0_F0_F0_F0_F0_F0,
		0xFF_00_FF_00_FF_00_FF_00,
		0xFF_FF_00_00_FF_FF_00_00,
		0xFF_FF_FF_FF_00_00_00_00
	];
	const EXCLUDES_BASIS_VEC_BITS: [Self; 7] = [
		!Self::CONTAINS_BASIS_VEC_BITS[0],
		!Self::CONTAINS_BASIS_VEC_BITS[1],
		!Self::CONTAINS_BASIS_VEC_BITS[2],
		!Self::CONTAINS_BASIS_VEC_BITS[3],
		!Self::CONTAINS_BASIS_VEC_BITS[4],
		!Self::CONTAINS_BASIS_VEC_BITS[5],
		!Self::CONTAINS_BASIS_VEC_BITS[6],
	];
	const OUTPUT_REVERSION_NEGATION_BITS: [Self; 7] = [
		Self::K[1] ^ Self::K[2] ^ Self::K[3] ^ Self::K[4] ^ Self::K[5] ^ Self::K[6],
						Self::K[2] ^ Self::K[3] ^ Self::K[4] ^ Self::K[5] ^ Self::K[6],
									Self::K[3] ^ Self::K[4] ^ Self::K[5] ^ Self::K[6],
													Self::K[4] ^ Self::K[5] ^ Self::K[6],
																Self::K[5] ^ Self::K[6],
																			Self::K[6],
																			0,
	];
}


pub struct SeaBasisBladeBitsIter<I:SeaBasisBladeBitsInner>(I);
impl<I:SeaBasisBladeBitsInner> Iterator for SeaBasisBladeBitsIter<I> {
	type Item = SeaBasisBladeIndexLow;

	fn next(&mut self) -> Option<Self::Item> {
		if self.0 == I::ZERO {
			return None;
		}
		let out = self.0.trailing_zeros();
		self.0 ^= I::ONE << out;
		Some(SeaBasisBladeIndexLow(out as u8))
	}

	fn size_hint(&self) -> (usize, Option<usize>) {
		(0, Some(self.0.count_ones() as usize))
	}

	fn count(self) -> usize
		where
			Self: Sized, {
		self.0.count_ones() as usize
	} 
}


#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SeaBasisBladeBits<I: SeaBasisBladeBitsInner>(pub I);
impl<I: SeaBasisBladeBitsInner> SeaBasisBladeBits<I> {
	pub const DIMENSIONS: u8 = I::DIMENSIONS;
	pub const NUM_BASIS_BLADES: SeaBasisBladeIndexLow = SeaBasisBladeIndexLow(I::NUM_BASIS_BLADES);
	pub const ZERO: Self = SeaBasisBladeBits(I::ZERO);
	pub const ONE: Self = SeaBasisBladeBits(I::ONE);
	pub const EVEN_GRADES: Self = SeaBasisBladeBits(I::EVEN_GRADE_BITS);
	pub const ODD_GRADES: Self = SeaBasisBladeBits(I::ODD_GRADE_BITS);
	pub const REVERSION_NEGATIONS: Self = SeaBasisBladeBits(I::REVERSION_NEGATION_BITS);
	pub const CONJUGATION_NEGATIONS: Self = SeaBasisBladeBits(I::CONJUGATION_NEGATION_BITS);
	pub const RIGHT_COMPLIMENT_NEGATIONS: Self = SeaBasisBladeBits(I::RIGHT_COMPLIMENT_NEGATION_BITS);
	pub const LEFT_COMPLIMENT_NEGATIONS: Self = SeaBasisBladeBits(I::LEFT_COMPLIMENT_NEGATION_BITS);
	pub const DUAL_NEGATIONS: Self = SeaBasisBladeBits(I::DUAL_NEGATION_BITS);
	pub const GRADES: [Self; 7] = [
		SeaBasisBladeBits(I::GRADE_BITS[0]),
		SeaBasisBladeBits(I::GRADE_BITS[1]),
		SeaBasisBladeBits(I::GRADE_BITS[2]),
		SeaBasisBladeBits(I::GRADE_BITS[3]),
		SeaBasisBladeBits(I::GRADE_BITS[4]),
		SeaBasisBladeBits(I::GRADE_BITS[5]),
		SeaBasisBladeBits(I::GRADE_BITS[6]),
	];
	pub const CONTAINS_BASIS_VECS: [Self; 7] = [
		SeaBasisBladeBits(I::CONTAINS_BASIS_VEC_BITS[0]),
		SeaBasisBladeBits(I::CONTAINS_BASIS_VEC_BITS[1]),
		SeaBasisBladeBits(I::CONTAINS_BASIS_VEC_BITS[2]),
		SeaBasisBladeBits(I::CONTAINS_BASIS_VEC_BITS[3]),
		SeaBasisBladeBits(I::CONTAINS_BASIS_VEC_BITS[4]),
		SeaBasisBladeBits(I::CONTAINS_BASIS_VEC_BITS[5]),
		SeaBasisBladeBits(I::CONTAINS_BASIS_VEC_BITS[6]),
	];
	pub const EXCLUDE_BASIS_VECS: [Self; 7] = [
		SeaBasisBladeBits(I::EXCLUDES_BASIS_VEC_BITS[0]),
		SeaBasisBladeBits(I::EXCLUDES_BASIS_VEC_BITS[1]),
		SeaBasisBladeBits(I::EXCLUDES_BASIS_VEC_BITS[2]),
		SeaBasisBladeBits(I::EXCLUDES_BASIS_VEC_BITS[3]),
		SeaBasisBladeBits(I::EXCLUDES_BASIS_VEC_BITS[4]),
		SeaBasisBladeBits(I::EXCLUDES_BASIS_VEC_BITS[5]),
		SeaBasisBladeBits(I::EXCLUDES_BASIS_VEC_BITS[6]),
	];
	pub const OUTPUT_REVERSION_NEGATIONS: [Self; 7] = [
		SeaBasisBladeBits(I::OUTPUT_REVERSION_NEGATION_BITS[0]),
		SeaBasisBladeBits(I::OUTPUT_REVERSION_NEGATION_BITS[1]),
		SeaBasisBladeBits(I::OUTPUT_REVERSION_NEGATION_BITS[2]),
		SeaBasisBladeBits(I::OUTPUT_REVERSION_NEGATION_BITS[3]),
		SeaBasisBladeBits(I::OUTPUT_REVERSION_NEGATION_BITS[4]),
		SeaBasisBladeBits(I::OUTPUT_REVERSION_NEGATION_BITS[5]),
		SeaBasisBladeBits(I::OUTPUT_REVERSION_NEGATION_BITS[6]),
	];

	pub fn zero_mask_for_output_index(index: &SeaBasisBladeIndexLow) -> Self {
		let mut out = I::ONE;
		for bv in index {
			out ^= out << (1 << bv);
		}
		Self(out)
	}

	pub fn basis_blades_less_than(index: &SeaBasisBladeIndexLow) -> Self {
		Self(I::set_lsb_up_to_n(**index))
	}

	pub fn cannonical_reordering_neg_mask(index: &SeaBasisBladeIndexLow) -> Self {
		let mut out = I::REVERSION_NEGATION_BITS;
		for i in index {
			let temp = I::OUTPUT_REVERSION_NEGATION_BITS[(i+1) as usize];
			out ^= temp;
		}
		Self(out)
	}

	pub fn metric_neg_mask(index: &SeaBasisBladeIndexLow, metric_negs: &SeaBasisBladeIndexLow) -> Self {
		let curr_met_negs = &!index & metric_negs;
		let mut out = I::ZERO;
		for i in &curr_met_negs {
			let temp = I::CONTAINS_BASIS_VEC_BITS[(i+1) as usize];
			out ^= temp;
		}
		Self(out)
	}

	pub fn dual_metric_neg_mask(index: &SeaBasisBladeIndexLow) -> Self {
		let mut out = I::ZERO;
		for i in index {
			let temp = !I::CONTAINS_BASIS_VEC_BITS[(i+1) as usize];
			out ^= temp;
		}
		Self(out)
	}

	pub fn remove_basis_blades_with_any(&mut self, index: &SeaBasisBladeIndexLow) {
		for i in index {
			*self &= &Self::EXCLUDE_BASIS_VECS[(i+1) as usize];
		}
	}

	pub fn remove_basis_blades_without_all(&mut self, index: &SeaBasisBladeIndexLow) {
		for i in index {
			*self &= &Self::CONTAINS_BASIS_VECS[(i+1) as usize];
		}
	}
	
	pub fn is_zero(&self) -> bool {
		*self == Self::ZERO
	}

	pub fn has(&self, index: &SeaBasisBladeIndexLow) -> bool{
		!(self & index).is_zero()		
	}
	
	pub fn set(&mut self, index: &SeaBasisBladeIndexLow, value:bool) {
		if value {
			*self |= index
		}
		else {
			*self &= &Self::from(&(!index))
		}
	}

	pub fn toggle(&mut self, index: &SeaBasisBladeIndexLow) {
		*self ^= &Self::from(index);
	}

	pub fn count_basis_blades_less_than(&self, idx:&SeaBasisBladeIndexLow) -> u32{
		let mask = I::set_lsb_up_to_n(**idx);
		(**self & mask).count_ones()
	}

	pub fn reverse_bits(&self) -> Self {
		Self(self.deref().reverse_bits())
	}

	pub fn geometric_product_swap(&self, swap_by:&SeaBasisBladeIndexLow) -> Self{
		let mut curr = self.0;
		for i in swap_by {
			let bv_idx = i+1;
			curr.swap_adjacent_pow2_n_bit_groups(bv_idx);
		}
		SeaBasisBladeBits(curr)
	}

	pub fn dual(&self) -> Self {
		Self(self.0.reverse_bits())
	}
}

impl<I: SeaBasisBladeBitsInner> Deref for SeaBasisBladeBits<I> {
	type Target = I;
	fn deref(&self) -> &Self::Target {
		&self.0
	}
}
impl<I: SeaBasisBladeBitsInner> DerefMut for SeaBasisBladeBits<I> {
	fn deref_mut(&mut self) -> &mut Self::Target {
		&mut self.0
	}
}

impl<I: SeaBasisBladeBitsInner> IntoIterator for &SeaBasisBladeBits<I> {
	type IntoIter = SeaBasisBladeBitsIter<I>;
	type Item = SeaBasisBladeIndexLow;
	fn into_iter(self) -> Self::IntoIter {
		SeaBasisBladeBitsIter(self.0)
	}
}

impl<I: SeaBasisBladeBitsInner> From<&SeaBasisBladeIndexLow> for SeaBasisBladeBits<I> {
	fn from(value: &SeaBasisBladeIndexLow) -> Self {
		let z = I::ONE << **value;
		SeaBasisBladeBits(z)
	}
}

impl<I: SeaBasisBladeBitsInner> Not for SeaBasisBladeBits<I> {
	type Output = Self;
	fn not(self) -> Self::Output {
		Self(!*self)
	}
}
impl<I: SeaBasisBladeBitsInner> Not for &SeaBasisBladeBits<I> {
	type Output = SeaBasisBladeBits<I>;
	fn not(self) -> Self::Output {
		SeaBasisBladeBits(!**self)
	}
}

impl<I: SeaBasisBladeBitsInner> BitOrAssign<&SeaBasisBladeIndexLow> for SeaBasisBladeBits<I> {
	fn bitor_assign(&mut self, rhs: &SeaBasisBladeIndexLow) {
		*self |= &Self::from(rhs);
	}
}
impl<I: SeaBasisBladeBitsInner> BitOrAssign<&SeaBasisBladeBits<I>> for SeaBasisBladeBits<I> {
	fn bitor_assign(&mut self, rhs: &SeaBasisBladeBits<I>) {
		**self |= **rhs;
	}
}
impl<I: SeaBasisBladeBitsInner> BitOr<&SeaBasisBladeBits<I>> for &SeaBasisBladeBits<I> {
	type Output = SeaBasisBladeBits<I>;
	fn bitor(self, rhs: &SeaBasisBladeBits<I>) -> Self::Output {
		SeaBasisBladeBits(**self | **rhs)
	}
}

impl<I: SeaBasisBladeBitsInner> BitAndAssign<&SeaBasisBladeBits<I>> for SeaBasisBladeBits<I> {
	fn bitand_assign(&mut self, rhs: &SeaBasisBladeBits<I>) {
		**self &= **rhs;
	}
}
impl<I: SeaBasisBladeBitsInner> BitAnd<&SeaBasisBladeIndexLow> for &SeaBasisBladeBits<I> {
	type Output = SeaBasisBladeBits<I>;
	fn bitand(self, rhs: &SeaBasisBladeIndexLow) -> Self::Output {
		self & &SeaBasisBladeBits::<I>::from(rhs)
	}
}
impl<I: SeaBasisBladeBitsInner> BitAnd<&SeaBasisBladeBits<I>> for &SeaBasisBladeBits<I> {
	type Output = SeaBasisBladeBits<I>;
	fn bitand(self, rhs: &SeaBasisBladeBits<I>) -> Self::Output {
		SeaBasisBladeBits(**self & **rhs)
	}
}

impl<I: SeaBasisBladeBitsInner> BitXorAssign<&SeaBasisBladeBits<I>> for SeaBasisBladeBits<I> {
	fn bitxor_assign(&mut self, rhs: &SeaBasisBladeBits<I>) {
		**self ^= **rhs;
		
	}
}
impl<I: SeaBasisBladeBitsInner> BitXor<&SeaBasisBladeBits<I>> for &SeaBasisBladeBits<I> {
	type Output = SeaBasisBladeBits<I>;
	fn bitxor(self, rhs: &SeaBasisBladeBits<I>) -> Self::Output {
		SeaBasisBladeBits(**self ^ **rhs)
	}
}

impl<I: SeaBasisBladeBitsInner> ShlAssign<u8> for SeaBasisBladeBits<I>{
	fn shl_assign(&mut self, rhs: u8) {
		**self <<= rhs;
	}
}
impl<I: SeaBasisBladeBitsInner> Shl<u8> for &SeaBasisBladeBits<I>{
	type Output = SeaBasisBladeBits<I>;
	fn shl(self, rhs: u8) -> Self::Output {
		SeaBasisBladeBits(**self << rhs)
	}
}
impl<I: SeaBasisBladeBitsInner> Shl<&SeaBasisBladeIndexLow> for &SeaBasisBladeBits<I>{
	type Output = SeaBasisBladeBits<I>;
	fn shl(self, rhs: &SeaBasisBladeIndexLow) -> Self::Output {
		SeaBasisBladeBits(**self << **rhs)
	}
}

impl<I: SeaBasisBladeBitsInner> ShrAssign<u8> for SeaBasisBladeBits<I>{
	fn shr_assign(&mut self, rhs: u8) {
		**self >>= rhs;
	}
}
impl<I: SeaBasisBladeBitsInner> Shr<&SeaBasisBladeIndexLow> for &SeaBasisBladeBits<I>{
	type Output = SeaBasisBladeBits<I>;
	fn shr(self, rhs: &SeaBasisBladeIndexLow) -> Self::Output {
		SeaBasisBladeBits(**self >> **rhs)
	}
}

impl<I: SeaBasisBladeBitsInner> LeftContraction<(&SeaBasisBladeIndexLow, u8)> for SeaBasisBladeBits<I> {
	type Output = Self;
	fn left_contraction(&self, (index, num_metric_non_zeros): (&SeaBasisBladeIndexLow, u8)) -> Self::Output {
		if index.contains_basis_vectors_greater_than(num_metric_non_zeros) {
			return SeaBasisBladeBits::<I>::ZERO;
		}
		let mut curr = self.0 >> index.0;
		for i in index {
			let bv_idx = i+1;
			curr &= !I::CONTAINS_BASIS_VEC_BITS[bv_idx as usize];
		}
		SeaBasisBladeBits(curr)
	}
}
impl<I: SeaBasisBladeBitsInner> LeftContraction<(&Self, u8)> for SeaBasisBladeBits<I> {
	type Output = Self;
	fn left_contraction(&self, (rhs, num_metric_non_zeros): (&Self, u8)) -> Self::Output {
		let mut curr = Self::ZERO;
		for i in self {
			curr |= &rhs.left_contraction((&i, num_metric_non_zeros));
		}
		curr
	}
}

impl<I: SeaBasisBladeBitsInner> RightContraction<(&SeaBasisBladeIndexLow, u8)> for SeaBasisBladeBits<I> {
	type Output = Self;
	fn right_contraction(&self, (index, num_metric_non_zeros): (&SeaBasisBladeIndexLow, u8)) -> Self::Output {
		let mut curr = self.0.reverse_bits() >> (I::NUM_BASIS_BLADES - 1 - index.0);
		for i in 0..I::DIMENSIONS {
			let bv_idx = i+1;
			if index.has(bv_idx) {
				if bv_idx > num_metric_non_zeros {
					curr &= I::CONTAINS_BASIS_VEC_BITS[bv_idx as usize];
				}
			}
			else {
				curr &= !I::CONTAINS_BASIS_VEC_BITS[bv_idx as usize];
			}
		}		
		SeaBasisBladeBits(curr)
	}
}
impl<I: SeaBasisBladeBitsInner> RightContraction<(&Self, u8)> for SeaBasisBladeBits<I> {
	type Output = Self;
	fn right_contraction(&self, (rhs, num_metric_non_zeros): (&Self, u8)) -> Self::Output {
		rhs.left_contraction((self, num_metric_non_zeros))
	}
}

impl<I: SeaBasisBladeBitsInner> InnerProduct<(&SeaBasisBladeIndexLow, u8)> for SeaBasisBladeBits<I> {
	type Output = Self;
	fn inner_product(&self, (index, num_metric_non_zeros): (&SeaBasisBladeIndexLow, u8)) -> Self::Output {
		&self.left_contraction((index, num_metric_non_zeros)) | &self.right_contraction((index, num_metric_non_zeros))
	}
}
impl<I: SeaBasisBladeBitsInner> InnerProduct<(&Self, u8)> for SeaBasisBladeBits<I> {
	type Output = Self;
	fn inner_product(&self, (rhs, num_metric_non_zeros): (&Self, u8)) -> Self::Output {
		&self.left_contraction((rhs, num_metric_non_zeros)) | &self.right_contraction((rhs, num_metric_non_zeros))
	}
}

impl<I: SeaBasisBladeBitsInner> OuterProduct<&SeaBasisBladeIndexLow> for SeaBasisBladeBits<I> {
	type Output = Self;
	fn outer_product(&self, index: &SeaBasisBladeIndexLow) -> Self::Output {
		let mut curr = self.0 << index.0;
		for i in index {
			let bv_idx = i+1;
			curr &= I::CONTAINS_BASIS_VEC_BITS[bv_idx as usize];
		}
		SeaBasisBladeBits(curr)
	}
}
impl<I: SeaBasisBladeBitsInner> OuterProduct<&Self> for SeaBasisBladeBits<I> {
	type Output = Self;
	fn outer_product(&self, rhs: &Self) -> Self::Output {
		let mut curr = Self::ZERO;
		for i in self {
			curr |= &rhs.outer_product(&i);
		}
		curr
	}
}

impl<I: SeaBasisBladeBitsInner> GeometricProduct<(&SeaBasisBladeIndexLow, u8)> for SeaBasisBladeBits<I> {
	type Output = Self;
	fn geometric_product(&self, (index, num_metric_non_zeros): (&SeaBasisBladeIndexLow, u8)) -> Self::Output {
		let mut curr = self.0;
		for i in index {
			let bv_idx = i+1;
			curr.swap_adjacent_pow2_n_bit_groups(bv_idx);
			if bv_idx > num_metric_non_zeros {
				curr &= I::CONTAINS_BASIS_VEC_BITS[bv_idx as usize];
			}
		}
		SeaBasisBladeBits(curr)
	}
}
impl<I: SeaBasisBladeBitsInner> GeometricProduct<(&Self, u8)> for SeaBasisBladeBits<I> {
	type Output = Self;
	fn geometric_product(&self, (rhs, num_metric_non_zeros): (&Self, u8)) -> Self::Output {
		let mut curr = Self::ZERO;
		for i in self {
			curr |= &rhs.geometric_product((&i, num_metric_non_zeros));
		}
		curr
	}
}

#[cfg(test)]
mod tests {
    use std::time::Instant;

    use super::*;

	#[test]
	fn thing() {

		// for i in 0..32{
		// 	let a = SeaBasisBladeBits::<u32>::cannonical_reordering_neg_mask(&SeaBasisBladeIndexLow(i));
		// 	println!("{:#x}", a.0);
		// }
	
		let a = SeaBasisBladeBits::<u32>::geometric_product(&SeaBasisBladeBits(513 as u32), (&SeaBasisBladeBits(26752 as u32), 3 as u8));
		println!("{:?}", a);
	}

	#[test]
	fn left_contraction() {
		// 1 + e1 + e2 + e13 + e123
		let a = SeaBasisBladeBits(167_u8);
		// e12
		let b = SeaBasisBladeBits(8_u8);
		// e1 + e2 + e12
		let out = SeaBasisBladeBits(14_u8);
		assert_eq!(
			a.left_contraction((&b, 3)),
			out
		);

		// 1 + e13 + e123
		let a = SeaBasisBladeBits(161_u8);
		// e1 + e12 + e13 + e23 
		let b = SeaBasisBladeBits(106_u8);
		// 1 + e1 + e12 + e13 + e23
		let out = SeaBasisBladeBits(107_u8);
		assert_eq!(
			a.left_contraction((&b, 3)),
			out
		);

		// e1 + e2
		let a = SeaBasisBladeBits(6_u8);
		// e2 + e3
		let b = SeaBasisBladeBits(20_u8);
		// 1
		let out = SeaBasisBladeBits(1_u8);
		assert_eq!(
			a.left_contraction((&b, 3)),
			out
		);

		// e1 + e12 + e3 + e134 + e234
		let a = SeaBasisBladeBits(24602_u16);
		// e34
		let b = SeaBasisBladeBits(1_u16 << 12);
		// e4
		let out = SeaBasisBladeBits(256_u16);
		assert_eq!(
			a.left_contraction((&b, 3)),
			out
		);

		// 1 + e13 + e4 + e123
		let a = SeaBasisBladeBits(417_u16);
		// e1 + e12 + e134 + e234
		let b = SeaBasisBladeBits(24586_u16);
		// e1 + e12 + e4 + e134 + e234 
		let out = SeaBasisBladeBits(24842_u16);
		assert_eq!(
			a.left_contraction((&b, 3)),
			out
		);

		// e3 + e14 + e124
		let a = SeaBasisBladeBits(2576_u16);
		// e24 + e134
		let b = SeaBasisBladeBits(9216_u16);
		// e14
		let out = SeaBasisBladeBits(512_u16);
		assert_eq!(
			a.left_contraction((&b, 3)),
			out
		);
	}

	#[test]
	fn right_contraction() {
		// 1 + e1 + e2 + e13 + e123
		let a = SeaBasisBladeBits(167_u8);
		// e12
		let b = SeaBasisBladeBits(8_u8);
		// e3
		let out = SeaBasisBladeBits(16_u8);
		assert_eq!(
			a.right_contraction((&b, 3)),
			out
		);

		// 1 + e13 + e123
		let a = SeaBasisBladeBits(161_u8);
		// e1 + e12 + e13 + e23 
		let b = SeaBasisBladeBits(106_u8);
		// 1 + e1 + e2 + e3 + e23
		let out = SeaBasisBladeBits(87_u8);
		assert_eq!(
			a.right_contraction((&b, 3)),
			out
		);

		// e1 + e2
		let a = SeaBasisBladeBits(6_u8);
		// e2 + e3
		let b = SeaBasisBladeBits(20_u8);
		// 1
		let out = SeaBasisBladeBits(1_u8);
		assert_eq!(
			a.right_contraction((&b, 3)),
			out
		);

		// e1 + e12 + e3 + e134 + e234
		let a = SeaBasisBladeBits(24602_u16);
		// e34
		let b = SeaBasisBladeBits(1_u16 << 12);
		// 0
		let out = SeaBasisBladeBits(0_u16);
		assert_eq!(
			a.right_contraction((&b, 3)),
			out
		);

		// 1 + e13 + e4 + e123
		let a = SeaBasisBladeBits(417_u16);
		// e1 + e12 + e134 + e234
		let b = SeaBasisBladeBits(24586_u16);
		// e3 + e23
		let out = SeaBasisBladeBits(80_u16);
		assert_eq!(
			a.right_contraction((&b, 3)),
			out
		);

		// e3 + e14 + e124
		let a = SeaBasisBladeBits(2576_u16);
		// e24 + e134
		let b = SeaBasisBladeBits(9216_u16);
		// 0
		let out = SeaBasisBladeBits(0_u16);
		assert_eq!(
			a.right_contraction((&b, 3)),
			out
		);
	}

	#[test]
	fn inner_product() {
		// 1 + e1 + e2 + e13 + e123
		let a = SeaBasisBladeBits(167_u8);
		// e12
		let b = SeaBasisBladeBits(8_u8);
		// e1 + e2 + e12 + e3
		let out = SeaBasisBladeBits(30_u8);
		
		assert_eq!(
			a.inner_product((&b, 3)),
			out
		);

		// 1 + e13 + e123
		let a = SeaBasisBladeBits(161_u8);
		// e1 + e12 + e13 + e23 
		let b = SeaBasisBladeBits(106_u8);
		// 1 + e2 + e12 + e3 + e13 + e23
		let out = SeaBasisBladeBits(127_u8);
		assert_eq!(
			a.inner_product((&b, 3)),
			out
		);

		// e1 + e2
		let a = SeaBasisBladeBits(6_u8);
		// e2 + e3
		let b = SeaBasisBladeBits(20_u8);
		// 1
		let out = SeaBasisBladeBits(1_u8);
		assert_eq!(
			a.inner_product((&b, 3)),
			out
		);

		// e1 + e12 + e3 + e134 + e234
		let a = SeaBasisBladeBits(24602_u16);
		// e34
		let b = SeaBasisBladeBits(1_u16 << 12);
		// e4
		let out = SeaBasisBladeBits(256_u16);
		assert_eq!(
			a.inner_product((&b, 3)),
			out
		);

		// 1 + e13 + e4 + e123
		let a = SeaBasisBladeBits(417_u16);
		// e1 + e12 + e134 + e234
		let b = SeaBasisBladeBits(24586_u16);
		// e1 + e12 + e3  + e23 + e4 +  e134 + e234
		let out = SeaBasisBladeBits(24922_u16);
		assert_eq!(
			a.inner_product((&b, 3)),
			out
		);

		// e3 + e14 + e124
		let a = SeaBasisBladeBits(2576_u16);
		// e24 + e134
		let b = SeaBasisBladeBits(9216_u16);
		// e14
		let out = SeaBasisBladeBits(512_u16);
		assert_eq!(
			a.inner_product((&b, 3)),
			out
		);

		// 3d
		// usual 1.228µs
		// no right 1.162µs (105% faster)
		// 4d
		// usual 2.589µs
		// no right 1.948µs (132% faster)
		// 5d
		// usual 5.299µs
		// no right 3.488µs (151% faster)
		// 6d
		// usual 11.643µs
		// no right 6.624µs (175% faster)
		// let a = SeaBasisBladeBits((0_u64).overflowing_sub(1).0);
		// let b = SeaBasisBladeBits((0_u64).overflowing_sub(1).0);

		// let now = Instant::now();
		// let num_iters = 10_000;
		// let mut t = a.clone();
		// for _ in 0..num_iters {
		// 	t = t.inner_product((&b, 3));
		// }
		// let elapsed_time = now.elapsed();
		// println!("time time time {:?}", elapsed_time/num_iters);
	}

	#[test]
	fn outer_product() {
		// 1 + e1 + e2 + e13 + e123
		let a = SeaBasisBladeBits(167_u8);
		// e12
		let b = SeaBasisBladeBits(8_u8);
		// e12
		let out = SeaBasisBladeBits(8_u8);
		assert_eq!(
			a.outer_product(&b),
			out
		);

		// 1 + e13 + e123
		let a = SeaBasisBladeBits(161_u8);
		// e1 + e12 + e13 + e23 
		let b = SeaBasisBladeBits(106_u8);
		// e1 + e12 + e13 + e23
		let out = SeaBasisBladeBits(106_u8);
		assert_eq!(
			a.outer_product(&b),
			out
		);

		// e1 + e2
		let a = SeaBasisBladeBits(6_u8);
		// e2 + e3
		let b = SeaBasisBladeBits(20_u8);
		// e12 + e13 + e23
		let out = SeaBasisBladeBits(104_u8);
		assert_eq!(
			a.outer_product(&b),
			out
		);

		// e1 + e12 + e3 + e134 + e234
		let a = SeaBasisBladeBits(24602_u16);
		// e34
		let b = SeaBasisBladeBits(1_u16 << 12);
		// e134 + e1234
		let out = SeaBasisBladeBits(40960_u16);
		assert_eq!(
			a.outer_product(&b),
			out
		);

		// 1 + e13 + e4 + e123
		let a = SeaBasisBladeBits(417_u16);
		// e1 + e12 + e134 + e234
		let b = SeaBasisBladeBits(24586_u16);
		// e1 + e14 + e12 + e124 + e134 + e234
		let out = SeaBasisBladeBits(27146_u16);
		assert_eq!(
			a.outer_product(&b),
			out
		);

		// e3 + e14 + e124
		let a = SeaBasisBladeBits(2576_u16);
		// e24 + e134
		let b = SeaBasisBladeBits(9216_u16);
		// e234
		let out = SeaBasisBladeBits(16384_u16);
		assert_eq!(
			a.outer_product(&b),
			out
		);
	}

	#[test]
	fn geometric_product() {
		// 1 + e1 + e2 + e13 + e123
		let a = SeaBasisBladeBits(167_u8);
		// e12
		let b = SeaBasisBladeBits(8_u8);
		// e1 + e2 - e12 + e3 + e23 
		let out = SeaBasisBladeBits(94_u8);
		assert_eq!(
			a.geometric_product((&b, 3)),
			out
		);

		// 1 + e13 + e123
		let a = SeaBasisBladeBits(161_u8);
		// e1 + e12 + e13 + e23 
		let b = SeaBasisBladeBits(106_u8);
		// 1 + e1 + e2 + e12 + e3 + e13 + e23
		let out = SeaBasisBladeBits(127_u8);
		assert_eq!(
			a.geometric_product((&b, 3)),
			out
		);

		// e1 + e2
		let a = SeaBasisBladeBits(6_u8);
		// e2 + e3
		let b = SeaBasisBladeBits(20_u8);
		// 1 + e12 + e13 + e23
		let out = SeaBasisBladeBits(105_u8);
		assert_eq!(
			a.geometric_product((&b, 3)),
			out
		);

		// e1 + e12 + e3 + e134 + e234
		let a = SeaBasisBladeBits(24602_u16);
		// e34
		let b = SeaBasisBladeBits(1_u16 << 12);
		// e4 + e134 + e1234
		let out = SeaBasisBladeBits(41216_u16);
		assert_eq!(
			a.geometric_product((&b, 3)),
			out
		);

		// 1 + e13 + e4 + e123
		let a = SeaBasisBladeBits(417_u16);
		// e1 + e12 + e134 + e234
		let b = SeaBasisBladeBits(24586_u16);
		// e4 + e1 + e3 + e14 + e24 + e12 + e23 + e124 + e134 + e234
		let out = SeaBasisBladeBits(28506_u16);
		assert_eq!(
			a.geometric_product((&b, 3)),
			out
		);

		// e3 + e14 + e124
		let a = SeaBasisBladeBits(2576_u16);
		// e24 + e134
		let b = SeaBasisBladeBits(9216_u16);
		// e14 + e234
		let out = SeaBasisBladeBits(16896_u16);
		assert_eq!(
			a.geometric_product((&b, 3)),
			out
		);
	}
}
