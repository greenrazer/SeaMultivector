use std::str::FromStr;
use std::ops::{Deref, DerefMut};

use crate::errors::ParseMultivectorError;


#[derive(Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct SeaBasisBladeIndexLow(pub u8);

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