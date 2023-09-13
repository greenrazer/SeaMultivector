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
	pub const E1: Self = SeaBasisBladeIndexLow(1);
	pub const E2: Self = SeaBasisBladeIndexLow(2);
	pub const E12: Self = SeaBasisBladeIndexLow(3);
	pub const E3: Self = SeaBasisBladeIndexLow(4);
	pub const E13: Self = SeaBasisBladeIndexLow(5);
	pub const E23: Self = SeaBasisBladeIndexLow(6);
	pub const E123: Self = SeaBasisBladeIndexLow(7);
	pub const E4: Self = SeaBasisBladeIndexLow(8);
	pub const E14: Self = SeaBasisBladeIndexLow(9);
	pub const E24: Self = SeaBasisBladeIndexLow(10);
	pub const E124: Self = SeaBasisBladeIndexLow(11);
	pub const E34: Self = SeaBasisBladeIndexLow(12);
	pub const E134: Self = SeaBasisBladeIndexLow(13);
	pub const E234: Self = SeaBasisBladeIndexLow(14);
	pub const E1234: Self = SeaBasisBladeIndexLow(15);
	pub const E5: Self = SeaBasisBladeIndexLow(16);
	pub const E15: Self = SeaBasisBladeIndexLow(17);
	pub const E25: Self = SeaBasisBladeIndexLow(18);
	pub const E125: Self = SeaBasisBladeIndexLow(19);
	pub const E35: Self = SeaBasisBladeIndexLow(20);
	pub const E135: Self = SeaBasisBladeIndexLow(21);
	pub const E235: Self = SeaBasisBladeIndexLow(22);
	pub const E1235: Self = SeaBasisBladeIndexLow(23);
	pub const E45: Self = SeaBasisBladeIndexLow(24);
	pub const E145: Self = SeaBasisBladeIndexLow(25);
	pub const E245: Self = SeaBasisBladeIndexLow(26);
	pub const E1245: Self = SeaBasisBladeIndexLow(27);
	pub const E345: Self = SeaBasisBladeIndexLow(28);
	pub const E1345: Self = SeaBasisBladeIndexLow(29);
	pub const E2345: Self = SeaBasisBladeIndexLow(30);
	pub const E12345: Self = SeaBasisBladeIndexLow(31);
	pub const E6: Self = SeaBasisBladeIndexLow(32);
	pub const E7: Self = SeaBasisBladeIndexLow(64);
	pub const E8: Self = SeaBasisBladeIndexLow(128);
	pub const ALL: Self = SeaBasisBladeIndexLow(255);
}
impl SeaBasisBladeIndexLow {
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
