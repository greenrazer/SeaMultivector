use std::str::FromStr;
use std::ops::{Deref, DerefMut, BitXor, BitOr, Not, BitAnd, Add, Sub};

use crate::errors::ParseMultivectorError;

pub struct BasisBladeIndexLowIter(u8);
impl Iterator for BasisBladeIndexLowIter {
	type Item = u8;

	fn next(&mut self) -> Option<Self::Item> {
		if self.0 == 0 {
			return None;
		}
		let out = self.0.trailing_zeros();
		self.0 ^= 1 << out;
		Some(out as u8)
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

#[derive(Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct SeaBasisBladeIndexLow(pub u8);
impl SeaBasisBladeIndexLow {
	pub const SCALAR: Self = SeaBasisBladeIndexLow(0);
	pub const ALL: Self = SeaBasisBladeIndexLow(255);

	pub fn basis_vectors_less_than_or_equal(n: u8) -> Self {
		let inner = 1_u8.checked_shl(n as u32).unwrap_or(0).overflowing_sub(1).0;
		SeaBasisBladeIndexLow(inner)
	}

	pub fn basis_vectors_greater_than(n: u8) -> Self {
		let all = Self::ALL;
		let bvlt = Self::basis_vectors_less_than_or_equal(n);
		&all ^ &bvlt
	}

	pub fn contains_basis_vectors_greater_than(&self, n: u8) -> bool {
		let bvgt = Self::basis_vectors_greater_than(n);
		let a = self & &bvgt;
		a > Self::SCALAR
	}

	pub fn has(&self, index: u8) -> bool {
		(**self & (1 << (index-1))) > 0
	}
}

impl Deref for SeaBasisBladeIndexLow {
	type Target = u8;
	fn deref(&self) -> &Self::Target {
		&self.0
	}
}
impl DerefMut for SeaBasisBladeIndexLow {
	fn deref_mut(&mut self) -> &mut Self::Target {
		&mut self.0
	}
}

impl IntoIterator for &SeaBasisBladeIndexLow {
	type IntoIter = BasisBladeIndexLowIter;
	type Item = u8;

	fn into_iter(self) -> Self::IntoIter {
		BasisBladeIndexLowIter(self.0)
	}
}

impl FromStr for SeaBasisBladeIndexLow {
	type Err = ParseMultivectorError;
		
	fn from_str(s: &str) -> Result<Self, Self::Err> {
		if s.is_empty() {
			return Err(ParseMultivectorError::InvalidBasisBladeIndex);
		}
        let mut out: u8 = 0;
		let mut last: u8 = 0;

        for c in s.chars() {
			let bv = match c.to_digit(10) {
				Some(0) => return Err(ParseMultivectorError::InvalidBasisVectorZero),
				Some(8) => return Err(ParseMultivectorError::InvalidBasisVectorOutOfRange),
				Some(9) => return Err(ParseMultivectorError::InvalidBasisVectorOutOfRange),
				Some(d) => d as u8,
				None => return Err(ParseMultivectorError::InvalidBasisVector)
			};
			if bv == last {
				return Err(ParseMultivectorError::DuplicateBasisVector);
			}
			if bv < last {
				return Err(ParseMultivectorError::InvalidBasisVectorOrdering);
			}
			last = bv;

			out |= 1 << (bv-1);
        }
		Ok(SeaBasisBladeIndexLow(out))
	}
}

impl Not for &SeaBasisBladeIndexLow {
	type Output = SeaBasisBladeIndexLow;
	fn not(self) -> Self::Output {
		SeaBasisBladeIndexLow(!**self)
	}
}
impl Not for SeaBasisBladeIndexLow {
	type Output = SeaBasisBladeIndexLow;
	fn not(self) -> Self::Output {
		SeaBasisBladeIndexLow(!*self)
	}
}

impl BitXor<Self> for &SeaBasisBladeIndexLow {
	type Output = SeaBasisBladeIndexLow;

	fn bitxor(self, rhs: Self) -> Self::Output {
		SeaBasisBladeIndexLow(**self ^ **rhs)
	}
}

impl BitOr<Self> for &SeaBasisBladeIndexLow {
	type Output = SeaBasisBladeIndexLow;

	fn bitor(self, rhs: Self) -> Self::Output {
		SeaBasisBladeIndexLow(**self | **rhs)
	}
}

impl BitAnd<Self> for &SeaBasisBladeIndexLow {
	type Output = SeaBasisBladeIndexLow;

	fn bitand(self, rhs: Self) -> Self::Output {
		SeaBasisBladeIndexLow(**self & **rhs)
	}
}

impl Add<Self> for &SeaBasisBladeIndexLow {
	type Output = SeaBasisBladeIndexLow;

	fn add(self, rhs: Self) -> Self::Output {
		SeaBasisBladeIndexLow(**self + **rhs)
	}
}

impl Sub<Self> for &SeaBasisBladeIndexLow {
	type Output = SeaBasisBladeIndexLow;

	fn sub(self, rhs: Self) -> Self::Output {
		SeaBasisBladeIndexLow(**self - **rhs)
	}
}
