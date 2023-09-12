use std::str::FromStr;
use std::ops::{Deref, DerefMut, Shl, Shr, BitOr, BitXor, Mul, Not, BitAnd, Index, AddAssign, Add, Neg, Sub, Div, SubAssign, MulAssign, DivAssign};

use crate::errors::ParseMultivectorError;
use crate::grade::SeaGrade;
use crate::metric::SeaMetric;
use crate::coefficient::SeaCoefficient;
use crate::multivector_product_traits::{LeftContraction, RightContraction, InnerProduct, OuterProduct, GeometricProduct, RegressiveProduct, ScalarProduct};
use crate::multivector_traits::{Multivector, MultivectorMut};

use super::basis_blade::SeaBasisBladeLow;
use super::basis_blade_index::SeaBasisBladeIndexLow;
use super::without_metric::{SeaMultivectorLowSparseCoefficents, SparseMultivectorIter, SparseMultivectorGradeIter}; 
use super::basis_blade_set::SeaBasisBladeSet;

mod linear_operations;
mod non_geometric_linear_operations;
mod scalar_operations;

#[derive(Debug, PartialEq, Clone)]
pub struct SeaMultivectorLowSparseCoefficentsWithMetric<T:SeaCoefficient, I:SeaBasisBladeSet>{
	pub metric: SeaMetric,
	pub multivector:SeaMultivectorLowSparseCoefficents<T, I>
}
impl<T:SeaCoefficient, I:SeaBasisBladeSet> SeaMultivectorLowSparseCoefficentsWithMetric<T, I> {
	pub fn new(p:u8, n:u8, z:u8) -> Self {
		Self {
			metric: SeaMetric::new(p,n,z),
			multivector: SeaMultivectorLowSparseCoefficents::<T, I>::new()
		}
	}

	pub fn iter(&self) -> SparseMultivectorIter<'_, T, I> {
		self.multivector.iter()
	}

	pub fn grade_iter(&self, grade: &SeaGrade) -> SparseMultivectorGradeIter<'_, T, I> {
		self.multivector.grade_iter(grade)
	}

	pub fn remove_zeros(&mut self) {
		self.multivector.remove_zeros();
	}

	pub fn remove_approx_zeros(&mut self, epsilon: T) {
		self.multivector.remove_approx_zeros(epsilon);
	}

	pub fn round_to_decimals(&mut self, decimals:i32) {
		self.multivector.round_to_decimals(decimals);
	}
}

impl<T:SeaCoefficient, I:SeaBasisBladeSet> FromStr for SeaMultivectorLowSparseCoefficentsWithMetric<T, I> {
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

impl<T:SeaCoefficient, I:SeaBasisBladeSet> Multivector<T> for SeaMultivectorLowSparseCoefficentsWithMetric<T, I> {
	type MultivectorIndexType = SeaBasisBladeIndexLow;

	fn basis_blade_coefficient(&self, index: &Self::MultivectorIndexType) -> T {
		self.multivector.basis_blade_coefficient(index)
	}

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

	fn dual(&self) -> Self {
		Self {
			metric: self.metric.clone(),
			multivector: self.multivector.dual(&self.metric)
		}
	}

	fn undual(&self) -> Self {
		Self {
			metric: self.metric.clone(),
			multivector: self.multivector.undual(&self.metric)
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

	fn scalar_norm(&self) -> T {
		self.multivector.scalar_norm(&self.metric)
	}

	fn scalar_normalized(&self) -> Self {
		Self { 
			metric: self.metric.clone(), 
			multivector: self.multivector.scalar_normalized(&self.metric) 
		}
	}

	fn exp(&self) -> Self {
		Self {
			metric: self.metric.clone(),
			multivector: self.multivector.exp(&self.metric)
		}
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeSet> MultivectorMut<T> for SeaMultivectorLowSparseCoefficentsWithMetric<T, I> {
	fn reversion_conjugation(&mut self) {
		self.multivector.reversion_conjugation()
	}

	fn involution_conjugation(&mut self) {
		self.multivector.involution_conjugation()
	}

	fn clifford_conjugation(&mut self) {
		self.multivector.clifford_conjugation()
	}

	
	fn dualization(&mut self) {
		self.multivector.dualization(&self.metric);
	}
	
	fn undualization(&mut self) {
		self.multivector.undualization(&self.metric);
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

	fn scalar_normalize(&mut self) {
		self.multivector.scalar_normalize(&self.metric);
	}
}

impl<T:SeaCoefficient, I:SeaBasisBladeSet> Not for SeaMultivectorLowSparseCoefficentsWithMetric<T, I> {
	type Output = SeaMultivectorLowSparseCoefficentsWithMetric<T, I>;

	fn not(self) -> Self::Output {
		self.dual()
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeSet> Not for &SeaMultivectorLowSparseCoefficentsWithMetric<T, I> {
	type Output = SeaMultivectorLowSparseCoefficentsWithMetric<T, I>;

	fn not(self) -> Self::Output {
		self.dual()
	}
}

impl<T:SeaCoefficient, I:SeaBasisBladeSet> Neg for &SeaMultivectorLowSparseCoefficentsWithMetric<T, I> {
	type Output = SeaMultivectorLowSparseCoefficentsWithMetric<T, I>;

	fn neg(self) -> Self::Output {
		SeaMultivectorLowSparseCoefficentsWithMetric {
			metric: self.metric.clone(),
			multivector: -&self.multivector
		}
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeSet> Neg for SeaMultivectorLowSparseCoefficentsWithMetric<T, I> {
	type Output = SeaMultivectorLowSparseCoefficentsWithMetric<T, I>;

	fn neg(self) -> Self::Output {
		-&self
	}
}

impl<T:SeaCoefficient, I:SeaBasisBladeSet> Shl<Self> for SeaMultivectorLowSparseCoefficentsWithMetric<T, I> {
	type Output = SeaMultivectorLowSparseCoefficentsWithMetric<T, I>;
	fn shl(self, rhs: Self) -> Self::Output {
		self.left_contraction(rhs)
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeSet> Shl<&Self> for SeaMultivectorLowSparseCoefficentsWithMetric<T, I> {
	type Output = SeaMultivectorLowSparseCoefficentsWithMetric<T, I>;
	fn shl(self, rhs: &Self) -> Self::Output {
		self.left_contraction(rhs)
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeSet> Shl<SeaMultivectorLowSparseCoefficentsWithMetric<T, I>> for &SeaMultivectorLowSparseCoefficentsWithMetric<T, I> {
	type Output = SeaMultivectorLowSparseCoefficentsWithMetric<T, I>;
	fn shl(self, rhs: SeaMultivectorLowSparseCoefficentsWithMetric<T, I>) -> Self::Output {
		self.left_contraction(rhs)
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeSet> Shl<Self> for &SeaMultivectorLowSparseCoefficentsWithMetric<T, I> {
	type Output = SeaMultivectorLowSparseCoefficentsWithMetric<T, I>;
	fn shl(self, rhs: Self) -> Self::Output {
		self.left_contraction(rhs)
	}
}

impl<T:SeaCoefficient, I:SeaBasisBladeSet> Shr<Self> for SeaMultivectorLowSparseCoefficentsWithMetric<T, I> {
	type Output = SeaMultivectorLowSparseCoefficentsWithMetric<T, I>;
	fn shr(self, rhs: Self) -> Self::Output {
		self.right_contraction(rhs)
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeSet> Shr<&Self> for SeaMultivectorLowSparseCoefficentsWithMetric<T, I> {
	type Output = SeaMultivectorLowSparseCoefficentsWithMetric<T, I>;
	fn shr(self, rhs: &Self) -> Self::Output {
		self.right_contraction(rhs)
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeSet> Shr<SeaMultivectorLowSparseCoefficentsWithMetric<T, I>> for &SeaMultivectorLowSparseCoefficentsWithMetric<T, I> {
	type Output = SeaMultivectorLowSparseCoefficentsWithMetric<T, I>;
	fn shr(self, rhs: SeaMultivectorLowSparseCoefficentsWithMetric<T, I>) -> Self::Output {
		self.right_contraction(rhs)
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeSet> Shr<Self> for &SeaMultivectorLowSparseCoefficentsWithMetric<T, I> {
	type Output = SeaMultivectorLowSparseCoefficentsWithMetric<T, I>;
	fn shr(self, rhs: Self) -> Self::Output {
		self.right_contraction(rhs)
	}
}

impl<T:SeaCoefficient, I:SeaBasisBladeSet> BitOr<Self> for SeaMultivectorLowSparseCoefficentsWithMetric<T, I> {
	type Output = SeaMultivectorLowSparseCoefficentsWithMetric<T, I>;
	fn bitor(self, rhs: Self) -> Self::Output {
		self.inner_product(rhs)
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeSet> BitOr<&Self> for SeaMultivectorLowSparseCoefficentsWithMetric<T, I> {
	type Output = SeaMultivectorLowSparseCoefficentsWithMetric<T, I>;
	fn bitor(self, rhs: &Self) -> Self::Output {
		self.inner_product(rhs)
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeSet> BitOr<SeaMultivectorLowSparseCoefficentsWithMetric<T, I>> for &SeaMultivectorLowSparseCoefficentsWithMetric<T, I> {
	type Output = SeaMultivectorLowSparseCoefficentsWithMetric<T, I>;
	fn bitor(self, rhs: SeaMultivectorLowSparseCoefficentsWithMetric<T, I>) -> Self::Output {
		self.inner_product(rhs)
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeSet> BitOr<Self> for &SeaMultivectorLowSparseCoefficentsWithMetric<T, I> {
	type Output = SeaMultivectorLowSparseCoefficentsWithMetric<T, I>;
	fn bitor(self, rhs: Self) -> Self::Output {
		self.inner_product(rhs)
	}
}

impl<T:SeaCoefficient, I:SeaBasisBladeSet> BitXor<Self> for SeaMultivectorLowSparseCoefficentsWithMetric<T, I> {
	type Output = SeaMultivectorLowSparseCoefficentsWithMetric<T, I>;
	fn bitxor(self, rhs: Self) -> Self::Output {
		self.outer_product(rhs)
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeSet> BitXor<&Self> for SeaMultivectorLowSparseCoefficentsWithMetric<T, I> {
	type Output = SeaMultivectorLowSparseCoefficentsWithMetric<T, I>;
	fn bitxor(self, rhs: &Self) -> Self::Output {
		self.outer_product(rhs)
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeSet> BitXor<SeaMultivectorLowSparseCoefficentsWithMetric<T, I>> for &SeaMultivectorLowSparseCoefficentsWithMetric<T, I> {
	type Output = SeaMultivectorLowSparseCoefficentsWithMetric<T, I>;
	fn bitxor(self, rhs: SeaMultivectorLowSparseCoefficentsWithMetric<T, I>) -> Self::Output {
		self.outer_product(rhs)
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeSet> BitXor<Self> for &SeaMultivectorLowSparseCoefficentsWithMetric<T, I> {
	type Output = SeaMultivectorLowSparseCoefficentsWithMetric<T, I>;
	fn bitxor(self, rhs: Self) -> Self::Output {
		self.outer_product(rhs)
	}
}

impl<T:SeaCoefficient, I:SeaBasisBladeSet> BitAnd<Self> for SeaMultivectorLowSparseCoefficentsWithMetric<T, I> {
	type Output = SeaMultivectorLowSparseCoefficentsWithMetric<T, I>;
	
	fn bitand(self, rhs: Self) -> Self::Output {
		self.regressive_product(rhs)
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeSet> BitAnd<&Self> for SeaMultivectorLowSparseCoefficentsWithMetric<T, I> {
	type Output = SeaMultivectorLowSparseCoefficentsWithMetric<T, I>;
	
	fn bitand(self, rhs: &Self) -> Self::Output {
		self.regressive_product(rhs)
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeSet> BitAnd<SeaMultivectorLowSparseCoefficentsWithMetric<T, I>> for &SeaMultivectorLowSparseCoefficentsWithMetric<T, I> {
	type Output = SeaMultivectorLowSparseCoefficentsWithMetric<T, I>;
	
	fn bitand(self, rhs: SeaMultivectorLowSparseCoefficentsWithMetric<T, I>) -> Self::Output {
		self.regressive_product(rhs)
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeSet> BitAnd<Self> for &SeaMultivectorLowSparseCoefficentsWithMetric<T, I> {
	type Output = SeaMultivectorLowSparseCoefficentsWithMetric<T, I>;
	
	fn bitand(self, rhs: Self) -> Self::Output {
		self.regressive_product(rhs)
	}
}

impl<T:SeaCoefficient, I:SeaBasisBladeSet> Mul<Self> for SeaMultivectorLowSparseCoefficentsWithMetric<T, I> {
	type Output = SeaMultivectorLowSparseCoefficentsWithMetric<T, I>;
	fn mul(self, rhs: Self) -> Self::Output {
		self.geometric_product(rhs)
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeSet> Mul<&Self> for SeaMultivectorLowSparseCoefficentsWithMetric<T, I> {
	type Output = SeaMultivectorLowSparseCoefficentsWithMetric<T, I>;
	fn mul(self, rhs: &Self) -> Self::Output {
		self.geometric_product(rhs)
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeSet> Mul<SeaMultivectorLowSparseCoefficentsWithMetric<T, I>> for &SeaMultivectorLowSparseCoefficentsWithMetric<T, I> {
	type Output = SeaMultivectorLowSparseCoefficentsWithMetric<T, I>;
	fn mul(self, rhs: SeaMultivectorLowSparseCoefficentsWithMetric<T, I>) -> Self::Output {
		self.geometric_product(rhs)
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeSet> Mul<Self> for &SeaMultivectorLowSparseCoefficentsWithMetric<T, I> {
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
					basis_blade_bits: 255, 
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
					basis_blade_bits: 65535, 
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
					basis_blade_bits: 65535, 
					coefficients: vec![1_f32, -1.0, 2.0, -3.0, 4.0, -5.0, 6.0, -7.0, 8.0, -9.0, 10.0, -11.0, 12.0, -13.0, 14.0, -15.0] 
				}
			}
		);

		assert_eq!(SeaMultivectorLowSparseCoefficentsWithMetric::<f32, u8>::from_str("e3").unwrap_err(), ParseMultivectorError::MissingMetric);
		assert_eq!(SeaMultivectorLowSparseCoefficentsWithMetric::<f32, u8>::from_str("a>e3").unwrap_err(), ParseMultivectorError::InvalidMetric(crate::errors::ParseMetricError::InvalidNumPositive));
		assert_eq!(SeaMultivectorLowSparseCoefficentsWithMetric::<f32, u8>::from_str("3,4,5>e1").unwrap_err(), ParseMultivectorError::InvalidMetricDimensions);
	}

	#[test]
	fn iter() {
		let a = SeaMultivectorLowSparseCoefficentsWithMetric::<f32, u16>::from_str("4>1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e34 + 13e134 + 14e234 + 15e1234").unwrap();
		let mut iter = a.iter();
		assert_eq!(iter.next(), Some(SeaBasisBladeLow::from_str("1").unwrap()));
		assert_eq!(iter.next(), Some(SeaBasisBladeLow::from_str("1e1").unwrap()));
		assert_eq!(iter.next(), Some(SeaBasisBladeLow::from_str("2e2").unwrap()));
		assert_eq!(iter.next(), Some(SeaBasisBladeLow::from_str("3e12").unwrap()));
		assert_eq!(iter.next(), Some(SeaBasisBladeLow::from_str("4e3").unwrap()));
		assert_eq!(iter.next(), Some(SeaBasisBladeLow::from_str("5e13").unwrap()));
		assert_eq!(iter.next(), Some(SeaBasisBladeLow::from_str("6e23").unwrap()));
		assert_eq!(iter.next(), Some(SeaBasisBladeLow::from_str("7e123").unwrap()));
		assert_eq!(iter.next(), Some(SeaBasisBladeLow::from_str("8e4").unwrap()));
		assert_eq!(iter.next(), Some(SeaBasisBladeLow::from_str("9e14").unwrap()));
		assert_eq!(iter.next(), Some(SeaBasisBladeLow::from_str("10e24").unwrap()));
		assert_eq!(iter.next(), Some(SeaBasisBladeLow::from_str("11e124").unwrap()));
		assert_eq!(iter.next(), Some(SeaBasisBladeLow::from_str("12e34").unwrap()));
		assert_eq!(iter.next(), Some(SeaBasisBladeLow::from_str("13e134").unwrap()));
		assert_eq!(iter.next(), Some(SeaBasisBladeLow::from_str("14e234").unwrap()));
		assert_eq!(iter.next(), Some(SeaBasisBladeLow::from_str("15e1234").unwrap()));
		assert_eq!(iter.next(), None);
	}

	#[test]
	fn grade_iter() {
		let a = SeaMultivectorLowSparseCoefficentsWithMetric::<f32, u16>::from_str("4>1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e34 + 13e134 + 14e234 + 15e1234").unwrap();
		
		let mut iter = a.grade_iter(&SeaGrade(0));
		assert_eq!(iter.next(), Some(SeaBasisBladeLow::from_str("1").unwrap()));
		assert_eq!(iter.next(), None);

		let mut iter = a.grade_iter(&SeaGrade(1));
		assert_eq!(iter.next(), Some(SeaBasisBladeLow::from_str("1e1").unwrap()));
		assert_eq!(iter.next(), Some(SeaBasisBladeLow::from_str("2e2").unwrap()));
		assert_eq!(iter.next(), Some(SeaBasisBladeLow::from_str("4e3").unwrap()));
		assert_eq!(iter.next(), Some(SeaBasisBladeLow::from_str("8e4").unwrap()));
		assert_eq!(iter.next(), None);

		let mut iter = a.grade_iter(&SeaGrade(2));
		assert_eq!(iter.next(), Some(SeaBasisBladeLow::from_str("3e12").unwrap()));
		assert_eq!(iter.next(), Some(SeaBasisBladeLow::from_str("5e13").unwrap()));
		assert_eq!(iter.next(), Some(SeaBasisBladeLow::from_str("6e23").unwrap()));
		assert_eq!(iter.next(), Some(SeaBasisBladeLow::from_str("9e14").unwrap()));
		assert_eq!(iter.next(), Some(SeaBasisBladeLow::from_str("10e24").unwrap()));
		assert_eq!(iter.next(), Some(SeaBasisBladeLow::from_str("12e34").unwrap()));
		assert_eq!(iter.next(), None);

		let mut iter = a.grade_iter(&SeaGrade(3));
		assert_eq!(iter.next(), Some(SeaBasisBladeLow::from_str("7e123").unwrap()));
		assert_eq!(iter.next(), Some(SeaBasisBladeLow::from_str("11e124").unwrap()));
		assert_eq!(iter.next(), Some(SeaBasisBladeLow::from_str("13e134").unwrap()));
		assert_eq!(iter.next(), Some(SeaBasisBladeLow::from_str("14e234").unwrap()));
		assert_eq!(iter.next(), None);

		let mut iter = a.grade_iter(&SeaGrade(4));
		assert_eq!(iter.next(), Some(SeaBasisBladeLow::from_str("15e1234").unwrap()));
		assert_eq!(iter.next(), None);
	}

	#[test]
	fn dual() {
		let a = SeaMultivectorLowSparseCoefficentsWithMetric::<f32, u32>::from_str("2,3>1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e34 + 13e134 + 14e234 + 15e1234 + e5 + e15 + 2e25 + 3e125 + 4e35 + 5e135 + 6e235 + 7e1235 + 8e45 + 9e145 + 10e245 + 11e1245 + 12e345 + 13e1345 + 14e2345 + 15e12345").unwrap();
		let out_dual = SeaMultivectorLowSparseCoefficentsWithMetric::<f32, u32>::from_str("2,3>15 +14e1 -13e2 -11e3 +7e4 -15e5 -12e12 -10e13 +6e14 -14e15 +9e23 -5e24 +13e25 -3e34 +11e35 -7e45 +8e123 -4e124 +12e125 -2e134 +10e135 -6e145 +e234 -9e235 +5e245 +3e345 +e1234 -8e1235 +4e1245 +2e1345 -e2345 -e12345").unwrap();
		assert_eq!(
			!&a,
			out_dual
		);
		assert_eq!(
			!a,
			out_dual
		);
	}

	#[test]
	fn exp() {
		let a = SeaMultivectorLowSparseCoefficentsWithMetric::<f32, u8>::from_str("2,1>e12 + 2e13 + 3e23").unwrap();
		let out = SeaMultivectorLowSparseCoefficentsWithMetric::<f32, u8>::from_str("2,1>15.989522934 +4.6067419052e12 +9.2134838104e13 +13.8202257156e23").unwrap();
		for (a_bb, o_bb) in a.exp().iter().zip(out.iter()) {
			assert_eq!(a_bb.index, o_bb.index);
			let abs_diff = (a_bb.coefficient - o_bb.coefficient).abs();
			assert!(abs_diff < 1e-6, "left != right(epsilon={})\n{:?}\n{:?}", 1e-6, a.exp(), out);
		}

		let a = SeaMultivectorLowSparseCoefficentsWithMetric::<f32, u32>::from_str("4,1>-e15 +4e25+e45 -e35").unwrap();
		let out = SeaMultivectorLowSparseCoefficentsWithMetric::<f32, u32>::from_str("4,1>39.0919036865 -8.9653635025e15 +35.86145401e25 -8.9653635025e35 +8.9653635025e45").unwrap();
		for (a_bb, o_bb) in a.exp().iter().zip(out.iter()) {
			assert_eq!(a_bb.index, o_bb.index);
			let abs_diff = (a_bb.coefficient - o_bb.coefficient).abs();
			assert!(abs_diff < 1e-5, "left != right(epsilon={})\n{:?}\n{:?}", 1e-5, a.exp(), out);
		}

		let a = SeaMultivectorLowSparseCoefficentsWithMetric::<f32, u16>::from_str("3,0,1>e12 + 2e13 + 3e23").unwrap();
		let out = SeaMultivectorLowSparseCoefficentsWithMetric::<f32, u16>::from_str("3,0,1>-0.8252990842 -0.1509213299e12 -0.3018426597e13 -0.4527639747e23").unwrap();
		for (a_bb, o_bb) in a.exp().iter().zip(out.iter()) {
			assert_eq!(a_bb.index, o_bb.index);
			let abs_diff = (a_bb.coefficient - o_bb.coefficient).abs();
			assert!(abs_diff < 1e-6, "left != right(epsilon={})\n{:?}\n{:?}", 1e-6, a.exp(), out);
		}

		let a = SeaMultivectorLowSparseCoefficentsWithMetric::<f32, u16>::from_str("3,0,1>4e14 + 5e24 + 6e34").unwrap();
		let out = SeaMultivectorLowSparseCoefficentsWithMetric::<f32, u16>::from_str("3,0,1>1 +4e14 +5e24 +6e34").unwrap();
		for (a_bb, o_bb) in a.exp().iter().zip(out.iter()) {
			assert_eq!(a_bb.index, o_bb.index);
			let abs_diff = (a_bb.coefficient - o_bb.coefficient).abs();
			assert!(abs_diff < 1e-6, "left != right(epsilon={})\n{:?}\n{:?}", 1e-6, a.exp(), out);
		}
	}

	mod norm {
		use super::*;

		#[test]
		fn scalar_norm() {
			let a = SeaMultivectorLowSparseCoefficentsWithMetric::<f32, u16>::from_str("3,0,1>1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e34 + 13e134 + 14e234 + 15e1234").unwrap();
			assert_eq!(a.scalar_norm(), 11.8743418201);

			let a = SeaMultivectorLowSparseCoefficentsWithMetric::<f32, u16>::from_str("3,1,0>1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e34 + 13e134 + 14e234 + 15e1234").unwrap();
			assert_eq!(a.scalar_norm(), 30.9677251344);
		}

		#[test]
		fn normalize() {
			let a = SeaMultivectorLowSparseCoefficentsWithMetric::<f32, u16>::from_str("3,0,1>1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e34 + 13e134 + 14e234 + 15e1234").unwrap();
			let out = SeaMultivectorLowSparseCoefficentsWithMetric::<f32, u16>::from_str("3,0,1>0.084215194 +0.67372155e4 +0.084215194e1 +0.168430388e2 +0.3368607759e3 +0.7579367e14 +0.84215194e24 +1.0105823e34 +0.252645582e12 +0.4210759699e13 +0.5052911639e23 +0.9263671e124 +1.0947975e134 +1.1790127e234 +0.5895063e123 +1.2632279e1234").unwrap();
			for (a_bb, o_bb) in a.scalar_normalized().iter().zip(out.iter()) {
				assert_eq!(a_bb.index, o_bb.index);
				let abs_diff = (a_bb.coefficient - o_bb.coefficient).abs();
				assert!(abs_diff < 1e-7);
			}

			let mut a = SeaMultivectorLowSparseCoefficentsWithMetric::<f32, u16>::from_str("3,1,0>1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e34 + 13e134 + 14e234 + 15e1234").unwrap();
			let out = SeaMultivectorLowSparseCoefficentsWithMetric::<f32, u16>::from_str("3,1,0>0.0322916843 +0.0322916843e1 +0.0645833686e2 +0.1291667372e3 +0.2583334744e4 +0.0968750492e12 +0.1614584178e13 +0.290625155e14 +0.1937500983e23 +0.3229168355e24 +0.3875001967e34 +0.22604178e123 +0.3552085161e124 +0.41979188e134 +0.4520835876e234 +0.48437524e1234").unwrap();
			a.scalar_normalize();

			for (a_bb, o_bb) in a.iter().zip(out.iter()) {
				assert_eq!(a_bb.index, o_bb.index);
				let abs_diff = (a_bb.coefficient - o_bb.coefficient).abs();
				assert!(abs_diff < 1e-7);
			}
		}
	}

}
