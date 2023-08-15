use std::ops::{Deref, DerefMut, AddAssign, BitOrAssign, Shl, Add, BitOr, BitAndAssign, BitAnd};
use std::str::FromStr;

use funty::Unsigned;

use crate::errors::ParseMultivectorError;

use super::basis_blade_index::SeaBasisBladeIndexLow;


pub trait SeaBasisBladeBitsInner: Unsigned{
	const DIMENSIONS: u8;
	const NUM_BASIS_BLADES: u8;
	fn from_u8(other: u8) -> Self;
}
impl SeaBasisBladeBitsInner for u8 {
	const DIMENSIONS: u8 = 3;
	const NUM_BASIS_BLADES: u8 = 8;
	fn from_u8(other: u8) -> Self{
		return other as Self;
	}
}
impl SeaBasisBladeBitsInner for u16 {
	const DIMENSIONS: u8 = 4;
	const NUM_BASIS_BLADES: u8 = 16;
	fn from_u8(other: u8) -> Self{
		return other as Self;
	}
}
impl SeaBasisBladeBitsInner for u32 {
	const DIMENSIONS: u8 = 5;
	const NUM_BASIS_BLADES: u8 = 32;
	fn from_u8(other: u8) -> Self{
		return other as Self;
	}
}
impl SeaBasisBladeBitsInner for u64 {
	const DIMENSIONS: u8 = 6;
	const NUM_BASIS_BLADES: u8 = 64;
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
	pub const ZERO: Self = SeaBasisBladeBits(I::ZERO);
	pub const ONE: Self = SeaBasisBladeBits(I::ONE);
	
	pub fn has(&self, index: &SeaBasisBladeIndexLow) -> bool{
		!(self & index).is_zero()		
	}

	pub fn is_zero(&self) -> bool {
		*self == Self::ZERO
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
		let z = I::ONE << *value.deref();
		SeaBasisBladeBits(z)
	}
}

impl<I: SeaBasisBladeBitsInner> BitOrAssign<&SeaBasisBladeIndexLow> for SeaBasisBladeBits<I> {
	fn bitor_assign(&mut self, rhs: &SeaBasisBladeIndexLow) {
		*self |= &Self::from(rhs);
	}
}
impl<I: SeaBasisBladeBitsInner> BitOrAssign<&SeaBasisBladeBits<I>> for SeaBasisBladeBits<I> {
	fn bitor_assign(&mut self, rhs: &SeaBasisBladeBits<I>) {
		*self.deref_mut() |= rhs.deref();
	}
}
impl<I: SeaBasisBladeBitsInner> BitOr<&SeaBasisBladeBits<I>> for &SeaBasisBladeBits<I> {
	type Output = SeaBasisBladeBits<I>;
	fn bitor(self, rhs: &SeaBasisBladeBits<I>) -> Self::Output {
		SeaBasisBladeBits(*self.deref() | rhs.deref())
	}
}

impl<I: SeaBasisBladeBitsInner> BitAndAssign<&SeaBasisBladeBits<I>> for SeaBasisBladeBits<I> {
	fn bitand_assign(&mut self, rhs: &SeaBasisBladeBits<I>) {
		*self.deref_mut() &= rhs.deref();
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
		SeaBasisBladeBits(*self.deref() & rhs.deref())
	}
}