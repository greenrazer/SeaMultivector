use std::str::FromStr;
use std::ops::{Deref, DerefMut};

use crate::errors::ParseMultivectorError;
use crate::metric::SeaMetric;
use crate::coefficient::SeaCoefficient;

use super::without_metric::SeaMultivectorLowSparseCoefficents; 
use super::basis_blade_bits::SeaBasisBladeBitsInner;
use super::basis_blade_bits::SeaBasisBladeBits;


#[derive(Debug, PartialEq)]
pub struct SeaMultivectorLowSparseCoefficentsWithMetric<T:SeaCoefficient, I:SeaBasisBladeBitsInner>{
	metric: SeaMetric,
	multivector:SeaMultivectorLowSparseCoefficents<T, I>
}

impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> Deref for SeaMultivectorLowSparseCoefficentsWithMetric<T, I> {
	type Target = SeaMultivectorLowSparseCoefficents<T, I>;
	fn deref(&self) -> &Self::Target {
		&self.multivector
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> DerefMut for SeaMultivectorLowSparseCoefficentsWithMetric<T, I> {
	fn deref_mut(&mut self) -> &mut Self::Target {
		&mut self.multivector
	}
}

impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> FromStr for SeaMultivectorLowSparseCoefficentsWithMetric<T, I> {
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
		let mv = SeaMultivectorLowSparseCoefficents::from_str(mv_str)?;

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
    fn from_str() {
		let a = SeaMultivectorLowSparseCoefficentsWithMetric::<f32, u8>::from_str("3>1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123").unwrap();
		assert_eq!(
			a,
			SeaMultivectorLowSparseCoefficentsWithMetric{
				metric: SeaMetric::new(3,0,0),
				multivector: SeaMultivectorLowSparseCoefficents { 
					basis_blade_bits: SeaBasisBladeBits(255), 
					coefficients: vec![1_f32, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0] 
				}
			}
		);

		let a = SeaMultivectorLowSparseCoefficentsWithMetric::<f32, u16>::from_str("2,2>1 - e1 + 2e2 - 3e12 + 4e3 - 5e13 + 6e23 - 7e123 + 8e4 - 9e14 + 10e24 - 11e124 + 12e34 - 13e134 + 14e234 - 15e1234").unwrap();
		assert_eq!(
			a,
			SeaMultivectorLowSparseCoefficentsWithMetric{
				metric: SeaMetric::new(2,2,0),
				multivector: SeaMultivectorLowSparseCoefficents { 
					basis_blade_bits: SeaBasisBladeBits(65535), 
					coefficients: vec![1_f32, -1.0, 2.0, -3.0, 4.0, -5.0, 6.0, -7.0, 8.0, -9.0, 10.0, -11.0, 12.0, -13.0, 14.0, -15.0] 
				}
			}
		);

		let a = SeaMultivectorLowSparseCoefficentsWithMetric::<f32, u16>::from_str("0,0,4>1 - e1 + 2e2 - 3e12 + 4e3 - 5e13 + 6e23 - 7e123 + 8e4 - 9e14 + 10e24 - 11e124 + 12e34 - 13e134 + 14e234 - 15e1234").unwrap();
		assert_eq!(
			a,
			SeaMultivectorLowSparseCoefficentsWithMetric{
				metric: SeaMetric::new(0,0,4),
				multivector: SeaMultivectorLowSparseCoefficents { 
					basis_blade_bits: SeaBasisBladeBits(65535), 
					coefficients: vec![1_f32, -1.0, 2.0, -3.0, 4.0, -5.0, 6.0, -7.0, 8.0, -9.0, 10.0, -11.0, 12.0, -13.0, 14.0, -15.0] 
				}
			}
		);

		assert_eq!(SeaMultivectorLowSparseCoefficentsWithMetric::<f32, u8>::from_str("e3").unwrap_err(), ParseMultivectorError::MissingMetric);
		assert_eq!(SeaMultivectorLowSparseCoefficentsWithMetric::<f32, u8>::from_str("a>e3").unwrap_err(), ParseMultivectorError::InvalidMetric(crate::errors::ParseMetricError::InvalidNumPositive));
		assert_eq!(SeaMultivectorLowSparseCoefficentsWithMetric::<f32, u8>::from_str("3,4,5>e1").unwrap_err(), ParseMultivectorError::InvalidMetricDimensions);
	}
}
