use std::str::FromStr;

use crate::errors::ParseMultivectorError;
use crate::coefficient::SeaCoefficient;

use super::basis_blade_index::SeaBasisBladeIndexLow;


#[derive(Debug, PartialEq)]
pub struct SeaBasisBladeLow<T: SeaCoefficient>{
	pub index: SeaBasisBladeIndexLow,
	pub coefficient: T
}

impl<T: SeaCoefficient> FromStr for SeaBasisBladeLow<T> {
	type Err = ParseMultivectorError;

	fn from_str(s: &str) -> Result<Self, Self::Err> {
		let (raw_coefficient, raw_index) = match s.split_once('e') {
			Some(r) => r,
			None => {
				let parsed_coeff = match T::from_str(s) {
					Ok(a) => a,
					Err(_) => return Err(ParseMultivectorError::InvalidBasisBladeCoefficient)
				};
	
				return Ok(SeaBasisBladeLow{
					index: SeaBasisBladeIndexLow(0),
					coefficient: parsed_coeff
				});
			}
		};

		let parsed_coefficient= match raw_coefficient {
			"" | "+" => T::one(),
			"-" => -T::one(),
			r => match T::from_str(r) {
				Ok(a) => a,
				Err(_) => return Err(ParseMultivectorError::InvalidBasisBladeCoefficient)
			}
		};

		let parsed_index = match SeaBasisBladeIndexLow::from_str(raw_index) {
			Ok(a) => a,
			Err(e) => return Err(e)
		};
		Ok( Self {
			index: parsed_index,
			coefficient: parsed_coefficient
		})
	}
}

#[derive(Debug, PartialEq)]
pub struct SeaBasisBladeLowRef<'a, T: SeaCoefficient>{
	pub index: SeaBasisBladeIndexLow,
	pub coefficient: &'a T
}

#[derive(Debug, PartialEq)]
pub struct SeaBasisBladeLowRefMut<'a, T: SeaCoefficient>{
	pub index: SeaBasisBladeIndexLow,
	pub coefficient: &'a mut T
}

#[cfg(test)]
mod tests {
    use super::*;

	#[test]
	fn from_str() {
		assert_eq!(
			SeaBasisBladeLow::<f32>::from_str("3e12").unwrap(),
			SeaBasisBladeLow{index: SeaBasisBladeIndexLow(3), coefficient:3_f32}
		);
		assert_eq!(
			SeaBasisBladeLow::<f32>::from_str("3e5").unwrap(),
			SeaBasisBladeLow{index: SeaBasisBladeIndexLow(16), coefficient:3_f32}
		);
		assert_eq!(SeaBasisBladeLow::<f32>::from_str("e").unwrap_err(), ParseMultivectorError::InvalidBasisBladeIndex);
		assert_eq!(SeaBasisBladeLow::<f32>::from_str("ae3").unwrap_err(), ParseMultivectorError::InvalidBasisBladeCoefficient);
		assert_eq!(SeaBasisBladeLow::<f32>::from_str("ea").unwrap_err(), ParseMultivectorError::InvalidBasisVector);
		assert_eq!(SeaBasisBladeLow::<f32>::from_str("e0").unwrap_err(), ParseMultivectorError::InvalidBasisVectorZero);
		assert_eq!(SeaBasisBladeLow::<f32>::from_str("e21").unwrap_err(), ParseMultivectorError::InvalidBasisVectorOrdering);
		assert_eq!(SeaBasisBladeLow::<f32>::from_str("e8").unwrap_err(), ParseMultivectorError::InvalidBasisVectorOutOfRange);
		assert_eq!(SeaBasisBladeLow::<f32>::from_str("e11").unwrap_err(), ParseMultivectorError::DuplicateBasisVector);
	}
}
