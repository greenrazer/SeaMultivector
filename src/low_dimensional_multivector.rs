use std::str::FromStr;
use std::ops::{Deref, DerefMut};

use lazy_static::lazy_static;
use regex::Regex;

use crate::errors::ParseMultivectorError;
use crate::metric::SeaMetric;
use crate::coefficient::SeaCoefficient;

#[derive(Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct SeaBasisBladeIndexLow(u8);
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

#[derive(Debug, PartialEq)]
pub struct SeaBasisBladeLow<T: SeaCoefficient>{
	index: SeaBasisBladeIndexLow,
	coefficient: T
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
		return Ok( Self {
			index: parsed_index,
			coefficient: parsed_coefficient
		});
	}
}

mod basis_blade_bits {
	use std::ops::{Deref, DerefMut, AddAssign, BitOrAssign, Shl};
	use std::str::FromStr;
	
	use num_traits::{Unsigned, Zero, One};

	use crate::errors::ParseMultivectorError;

	use super::SeaBasisBladeIndexLow;

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
}
use basis_blade_bits::SeaBasisBladeBits;
pub use basis_blade_bits::SeaBasisBladeBitsWrapper;

#[derive(Debug, PartialEq)]
pub struct MultivectorLowSparseCoefficentsNoMetric<T:SeaCoefficient, I:SeaBasisBladeBits>{
	basis_blade_bits: I,
	coefficients: Vec<T>
}
impl<T:SeaCoefficient, I:SeaBasisBladeBits> FromStr for MultivectorLowSparseCoefficentsNoMetric<T, I> {
	type Err = ParseMultivectorError;
	
	fn from_str(s: &str) -> Result<Self, Self::Err> {
		lazy_static! {
			// captures (+) (123.4e)(1234), (+)  (e)(1234), ()()(1234.45), (-)(10.9e)(1234), ()(10.3e)(123)
			static ref BBS_REGEX: Regex = Regex::new(r"([+-])[ ]*((?:[0-9]*(?:\.[0-9]+)?e)?[\(0-9\)]+(?:\.[0-9]+)?)").unwrap();
		}

		let mut basis_blades = Vec::<SeaBasisBladeLow<T>>::new();

		for (_, [sign, core]) in BBS_REGEX.captures_iter(&format!("+{}", s)).map(|c| c.extract()) {
			let raw_bb_s = format!("{}{}", sign, core);
			let bb = match SeaBasisBladeLow::<T>::from_str(&raw_bb_s) {
				Ok(x) => x,
				Err(e) => return Err(e)
			};
			if *bb.index >= I::NUM_BASIS_BLADES {
				return Err(ParseMultivectorError::InvalidBasisBladeForAlgebra);
			}
			if basis_blades.iter().any(|i| i.index == bb.index) {
				return Err(ParseMultivectorError::DuplicateBasisBlade);
			}
			basis_blades.push(bb);
		}

		basis_blades.sort_by_key(|i| *i.index);

		let mut out_bbbs = I::ZERO;

		for bb in &basis_blades {
			out_bbbs += &bb.index;
		}

		Ok(Self {
			basis_blade_bits: out_bbbs,
			coefficients: basis_blades.iter().map(|i| i.coefficient).collect()
		})

	}
}


#[derive(Debug, PartialEq)]
pub struct MultivectorLowSparseCoefficents<T:SeaCoefficient, I:SeaBasisBladeBits>{
	metric: SeaMetric,
	multivector:MultivectorLowSparseCoefficentsNoMetric<T, I>
}
impl<T:SeaCoefficient, I:SeaBasisBladeBits> FromStr for MultivectorLowSparseCoefficents<T, I> {
	type Err = ParseMultivectorError;
	
	fn from_str(s: &str) -> Result<Self, Self::Err> {
		let (metric_str, mv_str) = match s.split_once(">") {
			Some(r) => r,
			None => return Err(ParseMultivectorError::MissingMetric)
		};
		let metric = match  SeaMetric::from_str(metric_str){
			Ok(v) => v,
			Err(e) => return Err(ParseMultivectorError::InvalidMetric(e))
		};
		if metric.dimensions != I::DIMENSIONS {
			return Err(ParseMultivectorError::InvalidMetricDimensions);
		}
		let mv = MultivectorLowSparseCoefficentsNoMetric::from_str(mv_str)?;

		Ok(Self {
			metric: metric,
			multivector: mv 
		})
	}
}

#[cfg(test)]
mod tests {
    use super::*;

	#[test]
	fn basis_blade_from_str() {
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

    #[test]
    fn multivector_low_sparse_coefficients_no_metric_from_str() {
		let a = MultivectorLowSparseCoefficentsNoMetric::<f32, SeaBasisBladeBitsWrapper<u8>>::from_str("1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123").unwrap();
		assert_eq!(
			a,
			MultivectorLowSparseCoefficentsNoMetric { 
				basis_blade_bits: SeaBasisBladeBitsWrapper(255), 
				coefficients: vec![1_f32, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0] 
			}
		);

		let a = MultivectorLowSparseCoefficentsNoMetric::<f32, SeaBasisBladeBitsWrapper<u16>>::from_str("1 - e1 + 2e2 - 3e12 + 4e3 - 5e13 + 6e23 - 7e123 + 8e4 - 9e14 + 10e24 - 11e124 + 12e34 - 13e134 + 14e234 - 15e1234").unwrap();
		assert_eq!(
			a,
			MultivectorLowSparseCoefficentsNoMetric { 
				basis_blade_bits: SeaBasisBladeBitsWrapper(65535), 
				coefficients: vec![1_f32, -1.0, 2.0, -3.0, 4.0, -5.0, 6.0, -7.0, 8.0, -9.0, 10.0, -11.0, 12.0, -13.0, 14.0, -15.0] 
			}
		);

		assert_eq!(MultivectorLowSparseCoefficentsNoMetric::<f32, SeaBasisBladeBitsWrapper<u8>>::from_str("e4").unwrap_err(), ParseMultivectorError::InvalidBasisBladeForAlgebra);
		assert_eq!(MultivectorLowSparseCoefficentsNoMetric::<f32, SeaBasisBladeBitsWrapper<u8>>::from_str("e21").unwrap_err(), ParseMultivectorError::InvalidBasisVectorOrdering);
		assert_eq!(MultivectorLowSparseCoefficentsNoMetric::<f32, SeaBasisBladeBitsWrapper<u8>>::from_str("e11").unwrap_err(), ParseMultivectorError::DuplicateBasisVector);
		assert_eq!(MultivectorLowSparseCoefficentsNoMetric::<f32, SeaBasisBladeBitsWrapper<u8>>::from_str("e1 + e1").unwrap_err(), ParseMultivectorError::DuplicateBasisBlade);
	}

    #[test]
    fn multivector_low_sparse_coefficients_from_str() {
		let a = MultivectorLowSparseCoefficents::<f32, SeaBasisBladeBitsWrapper<u8>>::from_str("3>1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123").unwrap();
		assert_eq!(
			a,
			MultivectorLowSparseCoefficents{
				metric: SeaMetric::new(3,0,0),
				multivector: MultivectorLowSparseCoefficentsNoMetric { 
					basis_blade_bits: SeaBasisBladeBitsWrapper(255), 
					coefficients: vec![1_f32, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0] 
				}
			}
		);

		let a = MultivectorLowSparseCoefficents::<f32, SeaBasisBladeBitsWrapper<u16>>::from_str("2,2>1 - e1 + 2e2 - 3e12 + 4e3 - 5e13 + 6e23 - 7e123 + 8e4 - 9e14 + 10e24 - 11e124 + 12e34 - 13e134 + 14e234 - 15e1234").unwrap();
		assert_eq!(
			a,
			MultivectorLowSparseCoefficents{
				metric: SeaMetric::new(2,2,0),
				multivector: MultivectorLowSparseCoefficentsNoMetric { 
					basis_blade_bits: SeaBasisBladeBitsWrapper(65535), 
					coefficients: vec![1_f32, -1.0, 2.0, -3.0, 4.0, -5.0, 6.0, -7.0, 8.0, -9.0, 10.0, -11.0, 12.0, -13.0, 14.0, -15.0] 
				}
			}
		);

		let a = MultivectorLowSparseCoefficents::<f32, SeaBasisBladeBitsWrapper<u16>>::from_str("0,0,4>1 - e1 + 2e2 - 3e12 + 4e3 - 5e13 + 6e23 - 7e123 + 8e4 - 9e14 + 10e24 - 11e124 + 12e34 - 13e134 + 14e234 - 15e1234").unwrap();
		assert_eq!(
			a,
			MultivectorLowSparseCoefficents{
				metric: SeaMetric::new(0,0,4),
				multivector: MultivectorLowSparseCoefficentsNoMetric { 
					basis_blade_bits: SeaBasisBladeBitsWrapper(65535), 
					coefficients: vec![1_f32, -1.0, 2.0, -3.0, 4.0, -5.0, 6.0, -7.0, 8.0, -9.0, 10.0, -11.0, 12.0, -13.0, 14.0, -15.0] 
				}
			}
		);

		assert_eq!(MultivectorLowSparseCoefficents::<f32, SeaBasisBladeBitsWrapper<u8>>::from_str("e3").unwrap_err(), ParseMultivectorError::MissingMetric);
		assert_eq!(MultivectorLowSparseCoefficents::<f32, SeaBasisBladeBitsWrapper<u8>>::from_str("a>e3").unwrap_err(), ParseMultivectorError::InvalidMetric(crate::errors::ParseMetricError::InvalidNumPositive));
		assert_eq!(MultivectorLowSparseCoefficents::<f32, SeaBasisBladeBitsWrapper<u8>>::from_str("3,4,5>e1").unwrap_err(), ParseMultivectorError::InvalidMetricDimensions);
	}
}
