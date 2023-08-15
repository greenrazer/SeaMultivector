use std::ops::{Deref, DerefMut, AddAssign, BitOrAssign, Shl};
use std::str::FromStr;

use num_traits::{Unsigned, Zero, One};

use crate::errors::ParseMultivectorError;

use super::basis_blade_index::SeaBasisBladeIndexLow;

pub trait SeaBasisBladeBitsInner: Unsigned + BitOrAssign + for<'a> Shl<&'a u8, Output = Self>{
	const ZERO: Self;
	const ONE: Self;
	fn from_u8(other: u8) -> Self;
}
impl SeaBasisBladeBitsInner for u8 {
	const ZERO: Self = 0;
	const ONE: Self = 1;
	fn from_u8(other: u8) -> Self{
		return other as Self;
	}
}
impl SeaBasisBladeBitsInner for u16 {
	const ZERO: Self = 0;
	const ONE: Self = 1;
	fn from_u8(other: u8) -> Self{
		return other as Self;
	}
}
impl SeaBasisBladeBitsInner for u32 {
	const ZERO: Self = 0;
	const ONE: Self = 1;
	fn from_u8(other: u8) -> Self{
		return other as Self;
	}
}
impl SeaBasisBladeBitsInner for u64 {
	const ZERO: Self = 0;
	const ONE: Self = 1;
	fn from_u8(other: u8) -> Self{
		return other as Self;
	}
}

#[derive(Debug, PartialEq)]
pub struct SeaBasisBladeBitsWrapper<I: SeaBasisBladeBitsInner>(pub I);
impl<I: SeaBasisBladeBitsInner> Deref for SeaBasisBladeBitsWrapper<I> {
	type Target = I;
	fn deref(&self) -> &Self::Target {
		&self.0
	}
}
impl<I: SeaBasisBladeBitsInner> DerefMut for SeaBasisBladeBitsWrapper<I> {
	fn deref_mut(&mut self) -> &mut Self::Target {
		&mut self.0
	}
}
impl<I: SeaBasisBladeBitsInner> AddAssign<&SeaBasisBladeIndexLow> for SeaBasisBladeBitsWrapper<I> {
	fn add_assign(&mut self, rhs: &SeaBasisBladeIndexLow) {
		*self.deref_mut() |= I::ONE << rhs.deref();
	}
}

pub trait SeaBasisBladeBits: Deref + DerefMut + for<'a> AddAssign<&'a SeaBasisBladeIndexLow> {
	const DIMENSIONS: u8;
	const NUM_BASIS_BLADES: u8 = 1 << Self::DIMENSIONS;
	const ZERO: Self;
	const ONE: Self;
}
impl SeaBasisBladeBits for SeaBasisBladeBitsWrapper<u8>{
	const DIMENSIONS: u8 = 3;
	const ZERO: Self = SeaBasisBladeBitsWrapper(0);
	const ONE: Self = SeaBasisBladeBitsWrapper(1);
}
impl SeaBasisBladeBits for SeaBasisBladeBitsWrapper<u16>{
	const DIMENSIONS: u8 = 4;
	const ZERO: Self = SeaBasisBladeBitsWrapper(0);
	const ONE: Self = SeaBasisBladeBitsWrapper(1);
}
impl SeaBasisBladeBits for SeaBasisBladeBitsWrapper<u32>{
	const DIMENSIONS: u8 = 5;
	const ZERO: Self = SeaBasisBladeBitsWrapper(0);
	const ONE: Self = SeaBasisBladeBitsWrapper(1);
}
impl SeaBasisBladeBits for SeaBasisBladeBitsWrapper<u64>{
	const DIMENSIONS: u8 = 6;
	const ZERO: Self = SeaBasisBladeBitsWrapper(0);
	const ONE: Self = SeaBasisBladeBitsWrapper(1);
}