use std::iter::repeat;
use std::ops::AddAssign;
use std::str::FromStr;

use lazy_static::lazy_static;
use regex::Regex;

use crate::coefficient::SeaCoefficient;
use crate::errors::ParseMultivectorError;

use super::basis_blade_bits::SeaBasisBladeBits;
use super::basis_blade_bits::SeaBasisBladeBitsWrapper;
use super::basis_blade::SeaBasisBladeLow;

#[derive(Debug, PartialEq)]
pub struct SeaMultivectorLowSparseCoefficentsNoMetric<T:SeaCoefficient, I:SeaBasisBladeBits>{
	pub basis_blade_bits: I,
	pub coefficients: Vec<T>
}
impl<T:SeaCoefficient, I:SeaBasisBladeBits> FromStr for SeaMultivectorLowSparseCoefficentsNoMetric<T, I> {
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

// impl<T:SeaCoefficient, I:SeaBasisBladeBits> AddAssign<&Self> for SeaMultivectorLowSparseCoefficentsNoMetric<T, I> {
// 	fn add_assign(&mut self, rhs: &Self) {
// 		let bb_out = &self.basis_blade_bits + &rhs.basis_blade_bits;

// 		let size_diff = bb_out.count_ones() as usize - self.coefficients.len();
// 		self.coefficients.splice(0..0, repeat(T::zero()).take(size_diff));

// 		let mut self_idx = 0;
// 		let mut rhs_idx = 0;
// 		for (i, idx) in bb_out.into_iter().enumerate() {
// 			if self.basis_blade_bits.has(&idx) {
// 				self.coefficients.swap(i, size_diff+self_idx);
// 				self_idx += 1;
// 			}
			
// 			if rhs.basis_blade_bits.has(&idx) {
// 				self.coefficients[i] += rhs.coefficients[rhs_idx];
// 				rhs_idx += 1;
// 			}
// 		}

// 		self.basis_blade_bits = bb_out;
// 	}
// }

#[cfg(test)]
mod tests{
	use super::*;

	#[test]
	fn multivector_low_sparse_coefficients_no_metric_from_str() {
		let a = SeaMultivectorLowSparseCoefficentsNoMetric::<f32, SeaBasisBladeBitsWrapper<u8>>::from_str("1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123").unwrap();
		assert_eq!(
			a,
			SeaMultivectorLowSparseCoefficentsNoMetric { 
				basis_blade_bits: SeaBasisBladeBitsWrapper(255), 
				coefficients: vec![1_f32, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0] 
			}
		);

		let a = SeaMultivectorLowSparseCoefficentsNoMetric::<f32, SeaBasisBladeBitsWrapper<u16>>::from_str("1 - e1 + 2e2 - 3e12 + 4e3 - 5e13 + 6e23 - 7e123 + 8e4 - 9e14 + 10e24 - 11e124 + 12e34 - 13e134 + 14e234 - 15e1234").unwrap();
		assert_eq!(
			a,
			SeaMultivectorLowSparseCoefficentsNoMetric { 
				basis_blade_bits: SeaBasisBladeBitsWrapper(65535), 
				coefficients: vec![1_f32, -1.0, 2.0, -3.0, 4.0, -5.0, 6.0, -7.0, 8.0, -9.0, 10.0, -11.0, 12.0, -13.0, 14.0, -15.0] 
			}
		);

		assert_eq!(SeaMultivectorLowSparseCoefficentsNoMetric::<f32, SeaBasisBladeBitsWrapper<u8>>::from_str("e4").unwrap_err(), ParseMultivectorError::InvalidBasisBladeForAlgebra);
		assert_eq!(SeaMultivectorLowSparseCoefficentsNoMetric::<f32, SeaBasisBladeBitsWrapper<u8>>::from_str("e21").unwrap_err(), ParseMultivectorError::InvalidBasisVectorOrdering);
		assert_eq!(SeaMultivectorLowSparseCoefficentsNoMetric::<f32, SeaBasisBladeBitsWrapper<u8>>::from_str("e11").unwrap_err(), ParseMultivectorError::DuplicateBasisVector);
		assert_eq!(SeaMultivectorLowSparseCoefficentsNoMetric::<f32, SeaBasisBladeBitsWrapper<u8>>::from_str("e1 + e1").unwrap_err(), ParseMultivectorError::DuplicateBasisBlade);
	}
}