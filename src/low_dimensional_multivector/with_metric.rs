use std::str::FromStr;
use std::ops::{Deref, DerefMut, Shl, Shr, BitOr, BitXor, Mul, Not};

use crate::errors::ParseMultivectorError;
use crate::metric::SeaMetric;
use crate::coefficient::SeaCoefficient;
use crate::multivector_product_traits::{LeftContraction, RightContraction, InnerProduct, OuterProduct, GeometricProduct};
use crate::multivector_traits::{Multivector, MultivectorMut};

use super::without_metric::SeaMultivectorLowSparseCoefficents; 
use super::basis_blade_bits::SeaBasisBladeBitsInner;
use super::basis_blade_bits::SeaBasisBladeBits;


#[derive(Debug, PartialEq)]
pub struct SeaMultivectorLowSparseCoefficentsWithMetric<T:SeaCoefficient, I:SeaBasisBladeBitsInner>{
	metric: SeaMetric,
	multivector:SeaMultivectorLowSparseCoefficents<T, I>
}
impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> SeaMultivectorLowSparseCoefficentsWithMetric<T, I> {
	pub fn new(p:u8, n:u8, z:u8) -> Self {
		Self {
			metric: SeaMetric::new(p,n,z),
			multivector: SeaMultivectorLowSparseCoefficents::<T, I>::new()
		}
	}

	pub fn dual(&self) -> Self {
		Self {
			metric: self.metric.clone(),
			multivector: self.multivector.dual(&self.metric)
		}
	}

	pub fn undual(&self) -> Self {
		Self {
			metric: self.metric.clone(),
			multivector: self.multivector.undual(&self.metric)
		}
	}

	pub fn dualization(&mut self) {
		self.multivector.dualization(&self.metric);
	}

	pub fn undualization(&mut self) {
		self.multivector.undualization(&self.metric);
	}
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
		let (metric_str, mv_str) = match s.split_once('>') {
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
			metric,
			multivector: mv 
		})
	}
}

impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> Multivector for SeaMultivectorLowSparseCoefficentsWithMetric<T, I> {
	fn reversion_conjugated(&self) -> Self {
		Self {
			metric: self.metric.clone(),
			multivector: self.multivector.reversion_conjugated()
		}
	}
	
	fn involution_conjugated(&self) -> Self {
		Self {
			metric: self.metric.clone(),
			multivector: self.multivector.involution_conjugated()
		}
	}
	
	fn clifford_conjugated(&self) -> Self {
		Self {
			metric: self.metric.clone(),
			multivector: self.multivector.clifford_conjugated()
		}
	}

	fn poincare_dual(&self) -> Self {
		Self {
			metric: self.metric.clone(),
			multivector: self.multivector.poincare_dual()
		}
	}

	fn right_complimented(&self) -> Self {
		Self {
			metric: self.metric.clone(),
			multivector: self.multivector.right_complimented()
		}
	}

	fn left_complimented(&self) -> Self {
		Self {
			metric: self.metric.clone(),
			multivector: self.multivector.left_complimented()
		}
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> MultivectorMut for SeaMultivectorLowSparseCoefficentsWithMetric<T, I> {
	fn reversion_conjugation(&mut self) {
		self.multivector.reversion_conjugation()
	}

	fn involution_conjugation(&mut self) {
		self.multivector.involution_conjugation()
	}

	fn clifford_conjugation(&mut self) {
		self.multivector.clifford_conjugation()
	}

	fn poincare_dualization(&mut self) {
		self.multivector.poincare_dualization()
	}

	fn right_compliment(&mut self) {
		self.multivector.right_compliment()
	}

	fn left_compliment(&mut self) {
		self.multivector.left_compliment()
	}
}

impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> Not for &SeaMultivectorLowSparseCoefficentsWithMetric<T, I> {
	type Output = SeaMultivectorLowSparseCoefficentsWithMetric<T, I>;

	fn not(self) -> Self::Output {
		self.dual()
	}
}

impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> LeftContraction<&Self> for SeaMultivectorLowSparseCoefficentsWithMetric<T, I> {
	type Output = Self;
	fn left_contraction(&self, rhs: &Self) -> Self::Output {
		assert_eq!(self.metric, rhs.metric);
		Self {
			metric: self.metric.clone(),
			multivector: self.multivector.left_contraction((&rhs.multivector, &self.metric))
		}
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> Shl<Self> for &SeaMultivectorLowSparseCoefficentsWithMetric<T, I> {
	type Output = SeaMultivectorLowSparseCoefficentsWithMetric<T, I>;
	fn shl(self, rhs: Self) -> Self::Output {
		self.left_contraction(rhs)
	}
}

impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> RightContraction<&Self> for SeaMultivectorLowSparseCoefficentsWithMetric<T, I> {
	type Output = Self;
	fn right_contraction(&self, rhs: &Self) -> Self::Output {
		assert_eq!(self.metric, rhs.metric);
		Self {
			metric: self.metric.clone(),
			multivector: self.multivector.right_contraction((&rhs.multivector, &self.metric))
		}
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> Shr<Self> for &SeaMultivectorLowSparseCoefficentsWithMetric<T, I> {
	type Output = SeaMultivectorLowSparseCoefficentsWithMetric<T, I>;
	fn shr(self, rhs: Self) -> Self::Output {
		self.right_contraction(rhs)
	}
}

impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> InnerProduct<&Self> for SeaMultivectorLowSparseCoefficentsWithMetric<T, I> {
	type Output = Self;
	fn inner_product(&self, rhs: &Self) -> Self::Output {
		assert_eq!(self.metric, rhs.metric);
		Self {
			metric: self.metric.clone(),
			multivector: self.multivector.inner_product((&rhs.multivector, &self.metric))
		}
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> BitOr<Self> for &SeaMultivectorLowSparseCoefficentsWithMetric<T, I> {
	type Output = SeaMultivectorLowSparseCoefficentsWithMetric<T, I>;
	fn bitor(self, rhs: Self) -> Self::Output {
		self.inner_product(rhs)
	}
}

impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> OuterProduct<&Self> for SeaMultivectorLowSparseCoefficentsWithMetric<T, I> {
	type Output = Self;
	fn outer_product(&self, rhs: &Self) -> Self::Output {
		assert_eq!(self.metric.dimensions, rhs.metric.dimensions);
		Self {
			metric: self.metric.clone(),
			multivector: self.multivector.outer_product(&rhs.multivector)
		}
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> BitXor<Self> for &SeaMultivectorLowSparseCoefficentsWithMetric<T, I> {
	type Output = SeaMultivectorLowSparseCoefficentsWithMetric<T, I>;
	fn bitxor(self, rhs: Self) -> Self::Output {
		self.outer_product(rhs)
	}
}

impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> GeometricProduct<&Self> for SeaMultivectorLowSparseCoefficentsWithMetric<T, I> {
	type Output = Self;
	fn geometric_product(&self, rhs: &Self) -> Self::Output {
		assert_eq!(self.metric, rhs.metric);
		Self {
			metric: self.metric.clone(),
			multivector: self.multivector.geometric_product((&rhs.multivector, &self.metric))
		}
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> Mul<Self> for &SeaMultivectorLowSparseCoefficentsWithMetric<T, I> {
	type Output = SeaMultivectorLowSparseCoefficentsWithMetric<T, I>;
	fn mul(self, rhs: Self) -> Self::Output {
		self.geometric_product(rhs)
	}
}

#[cfg(test)]
mod tests {
	use crate::tests::catch_unwind_silent;

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

	#[test]
	fn dual() {
		let a = SeaMultivectorLowSparseCoefficentsWithMetric::<f32, u32>::from_str("2,3>1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e34 + 13e134 + 14e234 + 15e1234 + e5 + e15 + 2e25 + 3e125 + 4e35 + 5e135 + 6e235 + 7e1235 + 8e45 + 9e145 + 10e245 + 11e1245 + 12e345 + 13e1345 + 14e2345 + 15e12345").unwrap();
		let out_dual = SeaMultivectorLowSparseCoefficentsWithMetric::<f32, u32>::from_str("2,3>15 +14e1 -13e2 -11e3 +7e4 -15e5 -12e12 -10e13 +6e14 -14e15 +9e23 -5e24 +13e25 -3e34 +11e35 -7e45 +8e123 -4e124 +12e125 -2e134 +10e135 -6e145 +e234 -9e235 +5e245 +3e345 +e1234 -8e1235 +4e1245 +2e1345 -e2345 -e12345").unwrap();
		assert_eq!(
			!&a,
			out_dual
		);
	}

	mod left_contraction {
		use super::*;

		#[test]
		fn works() {
			let a = SeaMultivectorLowSparseCoefficentsWithMetric::<f32, u64>::from_str("3,3> + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345 +e6 + 12e126 + 8e46").unwrap();
			let b = SeaMultivectorLowSparseCoefficentsWithMetric::<f32, u64>::from_str("3,3>  + 3e13 + 1e4 + 2e124  +4e345 +83e1235 + 70e1245 + 4e1345 + 8e2345 + 12e123456").unwrap();

			let out = SeaMultivectorLowSparseCoefficentsWithMetric::<f32, u64>::from_str("3,3>115+56e1+78e2+3e3-5e4+189e5+180e6-16e12+3e13-4e14+202e15+168e16+2e24-215e25-156e26-333e35-262e45+2e124-228e125-144e126-198e135-156e145+19e235+38e245-120e345-132e356-84e456+179e1235+70e1245+4e1345-120e1356-72e1456+8e2345+108e2356+60e2456-36e3456+12e12345+96e12356+48e12456-24e13456+12e23456+12e123456").unwrap();
			assert_eq!(
				a.left_contraction(&b),
				out
			);
		}
	
		#[test]
		fn panics_on_mismatched_metrics() {
			let a = SeaMultivectorLowSparseCoefficentsWithMetric::<f32, u64>::from_str("1,2,3> + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345 +e6 + 12e126 + 8e46").unwrap();
			let b = SeaMultivectorLowSparseCoefficentsWithMetric::<f32, u64>::from_str("3,2,1>  + 3e13 + 1e4 + 2e124  +4e345 +83e1235 + 70e1245 + 4e1345 + 8e2345 + 12e123456").unwrap();
	
			let result = catch_unwind_silent(|| a.left_contraction(&b));
			assert!(result.is_err());
		}
	}

	mod right_contraction {
		use super::*;

		#[test]
		fn works() {
			let a = SeaMultivectorLowSparseCoefficentsWithMetric::<f32, u64>::from_str("3,3> + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345 +e6 + 12e126 + 8e46").unwrap();
			let b = SeaMultivectorLowSparseCoefficentsWithMetric::<f32, u64>::from_str("3,3>  + 3e13 + 1e4 + 2e124  +4e345 +83e1235 + 70e1245 + 4e1345 + 8e2345 + 12e123456").unwrap();

			let out = SeaMultivectorLowSparseCoefficentsWithMetric::<f32, u64>::from_str("3,3>115+59e1-105e2+1050e3+1245e4+8e6-71e12-18e35-39e45+13e135+14e235+45e245+15e1235").unwrap();
			assert_eq!(
				a.right_contraction(&b),
				out
			);
		}
	
		#[test]
		fn panics_on_mismatched_metrics() {
			let a = SeaMultivectorLowSparseCoefficentsWithMetric::<f32, u64>::from_str("1,2,3> + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345 +e6 + 12e126 + 8e46").unwrap();
			let b = SeaMultivectorLowSparseCoefficentsWithMetric::<f32, u64>::from_str("3,2,1>  + 3e13 + 1e4 + 2e124  +4e345 +83e1235 + 70e1245 + 4e1345 + 8e2345 + 12e123456").unwrap();
	
			let result = catch_unwind_silent(|| a.right_contraction(&b));
			assert!(result.is_err());
		}
	}

	mod inner_product {
		use super::*;

		#[test]
		fn works() {
			let a = SeaMultivectorLowSparseCoefficentsWithMetric::<f32, u64>::from_str("3,3> + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345 +e6 + 12e126 + 8e46").unwrap();
			let b = SeaMultivectorLowSparseCoefficentsWithMetric::<f32, u64>::from_str("3,3>  + 3e13 + 1e4 + 2e124  +4e345 +83e1235 + 70e1245 + 4e1345 + 8e2345 + 12e123456").unwrap();

			let out = SeaMultivectorLowSparseCoefficentsWithMetric::<f32, u64>::from_str("3,3>115+115e1 -27e2 +1053e3 +1240e4 +189e5 +188e6 -87e12 +3e13 -4e14 +202e15 +168e16 +2e24 -215e25 -156e26 -351e35 -301e45 +2e124 -228e125 -144e126 -185e135 -156e145 +33e235 +83e245 -120e345 -132e356 -84e456 +194e1235 +70e1245 +4e1345 -120e1356 -72e1456 +8e2345 +108e2356 +60e2456 -36e3456 +12e12345 +96e12356 +48e12456 -24e13456 +12e23456 +12e123456​").unwrap();
			assert_eq!(
				a.inner_product(&b),
				out
			);
		}
	
		#[test]
		fn panics_on_mismatched_metrics() {
			let a = SeaMultivectorLowSparseCoefficentsWithMetric::<f32, u64>::from_str("1,2,3> + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345 +e6 + 12e126 + 8e46").unwrap();
			let b = SeaMultivectorLowSparseCoefficentsWithMetric::<f32, u64>::from_str("3,2,1>  + 3e13 + 1e4 + 2e124  +4e345 +83e1235 + 70e1245 + 4e1345 + 8e2345 + 12e123456").unwrap();
	
			let result = catch_unwind_silent(|| a.inner_product(&b));
			assert!(result.is_err());
		}
	}

	mod outer_product {
		use super::*;

		#[test]
		fn works() {
			let a = SeaMultivectorLowSparseCoefficentsWithMetric::<f32, u64>::from_str("3,3> + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345 +e6 + 12e126 + 8e46").unwrap();
			let b = SeaMultivectorLowSparseCoefficentsWithMetric::<f32, u64>::from_str("3,3>  + 3e13 + 1e4 + 2e124  +4e345 +83e1235 + 70e1245 + 4e1345 + 8e2345 + 12e123456").unwrap();

			let out = SeaMultivectorLowSparseCoefficentsWithMetric::<f32, u64>::from_str("3,3>e4 +3e13 +e14 +2e24 +4e34 -e46 -6e123 +5e124 +29e134 +3e136 +6e234 +4e345 -15e1234 +83e1235 +70e1245 -14e1246 +8e1345 +24e1346 +16e2345 -4e3456 -372e12345 +83e12356 +70e12456 +4e13456 +8e23456 -700e123456").unwrap();
			assert_eq!(
				a.outer_product(&b),
				out
			);
		}
	}

	mod geometric_product {
		use super::*;

		#[test]
		fn works() {
			let a = SeaMultivectorLowSparseCoefficentsWithMetric::<f32, u64>::from_str("3,3> + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345 +e6 + 12e126 + 8e46").unwrap();
			let b = SeaMultivectorLowSparseCoefficentsWithMetric::<f32, u64>::from_str("3,3>  + 3e13 + 1e4 + 2e124  +4e345 +83e1235 + 70e1245 + 4e1345 + 8e2345 + 12e123456").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficentsWithMetric::<f32, u64>::from_str("3,3>115+115e1-27e2+1053e3+1240e4+189e5+188e6-21e12+983e13+1159e14+202e15+168e16-919e23-1075e24-215e25-156e26+17e34-351e35-301e45+23e46-846e123-991e124-228e125-128e126+17e134-265e135+3e136-228e145-17e234+143e235-36e236+135e245+303e345-1160e356-924e456-15e1234+230e1235+124e1245-14e1246-560e1256+442e1345+24e1346-152e1356-72e1456-393e2345+44e2356+60e2456-40e3456-360e12345+179e12356+118e12456+76e13456-28e23456-700e123456").unwrap();
			assert_eq!(
				a.geometric_product(&b),
				out
			);
		}
	
		#[test]
		fn panics_on_mismatched_metrics() {
			let a = SeaMultivectorLowSparseCoefficentsWithMetric::<f32, u64>::from_str("1,2,3> + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345 +e6 + 12e126 + 8e46").unwrap();
			let b = SeaMultivectorLowSparseCoefficentsWithMetric::<f32, u64>::from_str("3,2,1>  + 3e13 + 1e4 + 2e124  +4e345 +83e1235 + 70e1245 + 4e1345 + 8e2345 + 12e123456").unwrap();
	
			let result = catch_unwind_silent(|| a.geometric_product(&b));
			assert!(result.is_err());
		}
	}
}
