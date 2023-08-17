use std::ops::{Deref, DerefMut, AddAssign, BitOrAssign, Shl, Add, BitOr, BitAndAssign, BitAnd, Not, ShlAssign, ShrAssign, BitXor};
use std::str::FromStr;

use funty::Unsigned;

use crate::errors::ParseMultivectorError;
use crate::multivector_product_traits::{GeometricProduct, LeftContration, RightContraction, InnerProduct, OuterProduct};

use super::basis_blade_index::SeaBasisBladeIndexLow;


pub trait SeaBasisBladeBitsInner: Unsigned{
	const DIMENSIONS: u8;
	const NUM_BASIS_BLADES: u8 = Self::BITS as u8;
	const EVEN_GRADE_BITS: Self;
	const ODD_GRADE_BITS: Self;
	const REVERSION_NEGATION_BITS: Self;
	const CONJUGATION_NEGATION_BITS: Self;
	const GRADE_BITS: [Self; 7];
	const CONTAINS_BASIS_VEC_BITS: [Self; 7];
	const K: [Self; 7] = Self::CONTAINS_BASIS_VEC_BITS;
	const OUTPUT_REVERSION_NEGATION_BITS: [Self; 7];

	fn from_u8(other: u8) -> Self;

	fn set_lsb_up_to_n(n: u8) -> Self {
		let one = Self::ONE;
		one.checked_shl(n as u32).unwrap_or(Self::ZERO).overflowing_sub(one).0
	}

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
	const OUTPUT_REVERSION_NEGATION_BITS: [Self; 7] = [
		Self::K[1] ^ Self::K[2] ^ Self::K[3] ^ Self::K[4] ^ Self::K[5] ^ Self::K[6],
					 Self::K[2] ^ Self::K[3] ^ Self::K[4] ^ Self::K[5] ^ Self::K[6],
					 			  Self::K[3] ^ Self::K[4] ^ Self::K[5] ^ Self::K[6],
					 			  			   Self::K[4] ^ Self::K[5] ^ Self::K[6],
					 			  			                Self::K[5] ^ Self::K[6],
					 			  			                			 Self::K[6],
					 			  			                			 0,
	];

	fn from_u8(other: u8) -> Self{
		return other as Self;
	}
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
	const OUTPUT_REVERSION_NEGATION_BITS: [Self; 7] = [
		Self::K[1] ^ Self::K[2] ^ Self::K[3] ^ Self::K[4] ^ Self::K[5] ^ Self::K[6],
					 Self::K[2] ^ Self::K[3] ^ Self::K[4] ^ Self::K[5] ^ Self::K[6],
					 			  Self::K[3] ^ Self::K[4] ^ Self::K[5] ^ Self::K[6],
					 			  			   Self::K[4] ^ Self::K[5] ^ Self::K[6],
					 			  			                Self::K[5] ^ Self::K[6],
					 			  			                			 Self::K[6],
					 			  			                			 0,
	];

	fn from_u8(other: u8) -> Self{
		return other as Self;
	}
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
	const OUTPUT_REVERSION_NEGATION_BITS: [Self; 7] = [
		Self::K[1] ^ Self::K[2] ^ Self::K[3] ^ Self::K[4] ^ Self::K[5] ^ Self::K[6],
					 Self::K[2] ^ Self::K[3] ^ Self::K[4] ^ Self::K[5] ^ Self::K[6],
					 			  Self::K[3] ^ Self::K[4] ^ Self::K[5] ^ Self::K[6],
					 			  			   Self::K[4] ^ Self::K[5] ^ Self::K[6],
					 			  			                Self::K[5] ^ Self::K[6],
					 			  			                			 Self::K[6],
					 			  			                			 0,
	];

	fn from_u8(other: u8) -> Self{
		return other as Self;
	}
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
	const OUTPUT_REVERSION_NEGATION_BITS: [Self; 7] = [
		Self::K[1] ^ Self::K[2] ^ Self::K[3] ^ Self::K[4] ^ Self::K[5] ^ Self::K[6],
						Self::K[2] ^ Self::K[3] ^ Self::K[4] ^ Self::K[5] ^ Self::K[6],
									Self::K[3] ^ Self::K[4] ^ Self::K[5] ^ Self::K[6],
													Self::K[4] ^ Self::K[5] ^ Self::K[6],
																Self::K[5] ^ Self::K[6],
																			Self::K[6],
																			0,
	];

	fn from_u8(other: u8) -> Self{
		return other as Self;
	}
}


pub struct BasisBladeBitsIter<I:SeaBasisBladeBitsInner>(I);
impl<I:SeaBasisBladeBitsInner> Iterator for BasisBladeBitsIter<I> {
	type Item = SeaBasisBladeIndexLow;

	fn next(&mut self) -> Option<Self::Item> {
		if self.0 == I::ZERO {
			return None;
		}
		let out = self.0.trailing_zeros();
		self.0 ^= I::ONE << out;
		return Some(SeaBasisBladeIndexLow(out as u8));
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


#[derive(Debug, Clone, PartialEq)]
pub struct SeaBasisBladeBits<I: SeaBasisBladeBitsInner>(pub I);
impl<I: SeaBasisBladeBitsInner> SeaBasisBladeBits<I> {
	pub const DIMENSIONS: u8 = I::DIMENSIONS;
	pub const NUM_BASIS_BLADES: u8 = I::NUM_BASIS_BLADES;
	pub const ZERO: Self = SeaBasisBladeBits(I::ZERO);
	pub const ONE: Self = SeaBasisBladeBits(I::ONE);
	pub const EVEN_GRADES: Self = SeaBasisBladeBits(I::EVEN_GRADE_BITS);
	pub const ODD_GRADES: Self = SeaBasisBladeBits(I::ODD_GRADE_BITS);
	pub const REVERSION_NEGATIONS: Self = SeaBasisBladeBits(I::REVERSION_NEGATION_BITS);
	pub const CONJUGATION_NEGATIONS: Self = SeaBasisBladeBits(I::CONJUGATION_NEGATION_BITS);
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

	pub fn cannonical_neg_mask(index: &SeaBasisBladeIndexLow) -> Self {
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
	
	pub fn is_zero(&self) -> bool {
		*self == Self::ZERO
	}

	pub fn has(&self, index: &SeaBasisBladeIndexLow) -> bool{
		!(self & index).is_zero()		
	}

	pub fn count_basis_blades_less_than(&self, idx:&SeaBasisBladeIndexLow) -> u32{
		let mask = I::set_lsb_up_to_n(**idx);
		(**self & &mask).count_ones()
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
	type IntoIter = BasisBladeBitsIter<I>;
	type Item = SeaBasisBladeIndexLow;
	fn into_iter(self) -> Self::IntoIter {
		BasisBladeBitsIter(self.0.clone())
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

impl<I: SeaBasisBladeBitsInner> ShrAssign<u8> for SeaBasisBladeBits<I>{
	fn shr_assign(&mut self, rhs: u8) {
		**self >>= rhs;
	}
}

impl<I: SeaBasisBladeBitsInner> LeftContration<(&SeaBasisBladeIndexLow, u8)> for SeaBasisBladeBits<I> {
	type Output = Self;
	fn left_contraction(&self, (index, num_metric_non_zeros): (&SeaBasisBladeIndexLow, u8)) -> Self::Output {
		if index.contains_basis_vectors_greater_than(num_metric_non_zeros) {
			return SeaBasisBladeBits::<I>::ZERO;
		}
		let mut curr = self.0.clone() >> index.0;
		for i in index {
			let bv_idx = i+1;
			curr &= !I::CONTAINS_BASIS_VEC_BITS[bv_idx as usize];
		}
		SeaBasisBladeBits(curr)
	}
}
impl<I: SeaBasisBladeBitsInner> LeftContration<(&Self, u8)> for SeaBasisBladeBits<I> {
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
		let mut curr = self.0.clone().reverse_bits() >> (I::NUM_BASIS_BLADES - 1 - index.0);
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
		let mut curr = self.0.clone() << index.0;
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
		let mut curr = self.0.clone();
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
