use std::str::FromStr;
use std::ops::{Deref, DerefMut};

use crate::errors::ParseMultivectorError;
use crate::metric::SeaMetric;
use crate::coefficient::SeaCoefficient;

use super::no_metric::SeaMultivectorLowSparseCoefficentsNoMetric; 
use super::basis_blade_bits::SeaBasisBladeBits;
use super::basis_blade_bits::SeaBasisBladeBitsWrapper;


#[derive(Debug, PartialEq)]
pub struct SeaMultivectorLowSparseCoefficents<T:SeaCoefficient, I:SeaBasisBladeBits>{
	metric: SeaMetric,
	multivector:SeaMultivectorLowSparseCoefficentsNoMetric<T, I>
}
impl<T:SeaCoefficient, I:SeaBasisBladeBits> FromStr for SeaMultivectorLowSparseCoefficents<T, I> {
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
		let mv = SeaMultivectorLowSparseCoefficentsNoMetric::from_str(mv_str)?;

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
    fn multivector_low_sparse_coefficients_from_str() {
		let a = SeaMultivectorLowSparseCoefficents::<f32, SeaBasisBladeBitsWrapper<u8>>::from_str("3>1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123").unwrap();
		assert_eq!(
			a,
			SeaMultivectorLowSparseCoefficents{
				metric: SeaMetric::new(3,0,0),
				multivector: SeaMultivectorLowSparseCoefficentsNoMetric { 
					basis_blade_bits: SeaBasisBladeBitsWrapper(255), 
					coefficients: vec![1_f32, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0] 
				}
			}
		);

		let a = SeaMultivectorLowSparseCoefficents::<f32, SeaBasisBladeBitsWrapper<u16>>::from_str("2,2>1 - e1 + 2e2 - 3e12 + 4e3 - 5e13 + 6e23 - 7e123 + 8e4 - 9e14 + 10e24 - 11e124 + 12e34 - 13e134 + 14e234 - 15e1234").unwrap();
		assert_eq!(
			a,
			SeaMultivectorLowSparseCoefficents{
				metric: SeaMetric::new(2,2,0),
				multivector: SeaMultivectorLowSparseCoefficentsNoMetric { 
					basis_blade_bits: SeaBasisBladeBitsWrapper(65535), 
					coefficients: vec![1_f32, -1.0, 2.0, -3.0, 4.0, -5.0, 6.0, -7.0, 8.0, -9.0, 10.0, -11.0, 12.0, -13.0, 14.0, -15.0] 
				}
			}
		);

		let a = SeaMultivectorLowSparseCoefficents::<f32, SeaBasisBladeBitsWrapper<u16>>::from_str("0,0,4>1 - e1 + 2e2 - 3e12 + 4e3 - 5e13 + 6e23 - 7e123 + 8e4 - 9e14 + 10e24 - 11e124 + 12e34 - 13e134 + 14e234 - 15e1234").unwrap();
		assert_eq!(
			a,
			SeaMultivectorLowSparseCoefficents{
				metric: SeaMetric::new(0,0,4),
				multivector: SeaMultivectorLowSparseCoefficentsNoMetric { 
					basis_blade_bits: SeaBasisBladeBitsWrapper(65535), 
					coefficients: vec![1_f32, -1.0, 2.0, -3.0, 4.0, -5.0, 6.0, -7.0, 8.0, -9.0, 10.0, -11.0, 12.0, -13.0, 14.0, -15.0] 
				}
			}
		);

		assert_eq!(SeaMultivectorLowSparseCoefficents::<f32, SeaBasisBladeBitsWrapper<u8>>::from_str("e3").unwrap_err(), ParseMultivectorError::MissingMetric);
		assert_eq!(SeaMultivectorLowSparseCoefficents::<f32, SeaBasisBladeBitsWrapper<u8>>::from_str("a>e3").unwrap_err(), ParseMultivectorError::InvalidMetric(crate::errors::ParseMetricError::InvalidNumPositive));
		assert_eq!(SeaMultivectorLowSparseCoefficents::<f32, SeaBasisBladeBitsWrapper<u8>>::from_str("3,4,5>e1").unwrap_err(), ParseMultivectorError::InvalidMetricDimensions);
	}
}
