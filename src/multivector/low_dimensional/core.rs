use std::cmp::Ordering;
use std::iter::repeat;
use std::ops::{Add, AddAssign, BitXor, Neg, Index, Mul, MulAssign, Sub, Div, DivAssign, SubAssign};
use std::str::FromStr;

use lazy_static::lazy_static;
use num_traits::{Zero, ToPrimitive};
use regex::Regex;

use crate::coefficient::{SeaCoefficient, self};
use crate::errors::ParseMultivectorError;

use crate::multivector::metric::SeaMetric;
use crate::multivector::multivector_traits::{Multivector, MultivectorMut};
use crate::multivector::multivector_product_traits::{LeftContraction, RightContraction, InnerProduct, OuterProduct, GeometricProduct, RegressiveProduct, ScalarProduct};
use crate::multivector::grade::SeaGrade;
use crate::multivector::basis_blade::SeaBasisBladeLow;
use crate::multivector::basis_blade_index::SeaBasisBladeIndexLow;

use super::basis_blade_set::{SeaBasisBladeSet, SeaBasisBladeSetIter};

mod scalar_operations;
mod basis_blade_operations;
mod linear_operations;
mod non_geometric_linear_operations;
mod non_linear_operations;

pub struct SparseMultivectorIter<'a, T:SeaCoefficient, I:SeaBasisBladeSet>{
	multivector: &'a SeaLowDimensionalSparseMultivector<T,I>,
	curr_iter: I::IntoIterType
}
impl<'a, T:SeaCoefficient, I:SeaBasisBladeSet> SparseMultivectorIter<'a, T, I> {
	pub fn new(multivector: &'a SeaLowDimensionalSparseMultivector<T,I>) -> Self{
		SparseMultivectorIter {
			multivector,
			curr_iter: multivector.basis_blade_set.iter()
		}
	}
}
impl<'a, T:SeaCoefficient, I:SeaBasisBladeSet> Iterator for SparseMultivectorIter<'a, T, I> {
	type Item = SeaBasisBladeLow<T>;
	fn next(&mut self) -> Option<Self::Item> {
		let next_index = self.curr_iter.next()?;
		let coeff = self.multivector.basis_blade_coefficient(&next_index);
		Some(SeaBasisBladeLow {
			index: next_index,
			coefficient: coeff
		})
	}

	fn size_hint(&self) -> (usize, Option<usize>) {
		self.curr_iter.size_hint()
	}

	fn count(self) -> usize
		where
			Self: Sized, {
		self.curr_iter.count()
	} 
}

pub struct SparseMultivectorGradeIter<'a, T:SeaCoefficient, I:SeaBasisBladeSet>{
	multivector: &'a SeaLowDimensionalSparseMultivector<T,I>,
	curr_iter: I::IntoIterType
}
impl<'a, T:SeaCoefficient, I:SeaBasisBladeSet> SparseMultivectorGradeIter<'a, T, I> {
	pub fn new(multivector: &'a SeaLowDimensionalSparseMultivector<T,I>, grade: &SeaGrade) -> Self{
		let iter_bits =  multivector.basis_blade_set.intersection(&I::grade(grade));
		SparseMultivectorGradeIter {
			multivector,
			curr_iter: iter_bits.iter()
		}
	}
}
impl<'a, T:SeaCoefficient, I:SeaBasisBladeSet> Iterator for SparseMultivectorGradeIter<'a, T, I> {
	type Item = SeaBasisBladeLow<T>;
	fn next(&mut self) -> Option<Self::Item> {
		let next_index = self.curr_iter.next()?;
		let coeff = self.multivector.basis_blade_coefficient(&next_index);
		Some(SeaBasisBladeLow {
			index: next_index,
			coefficient: coeff
		})
	}

	fn size_hint(&self) -> (usize, Option<usize>) {
		self.curr_iter.size_hint()
	}

	fn count(self) -> usize
		where
			Self: Sized, {
		self.curr_iter.count()
	} 
}


#[derive(Debug, Clone, PartialEq)]
pub struct SeaLowDimensionalSparseMultivector<T:SeaCoefficient, I:SeaBasisBladeSet>{
	pub basis_blade_set: I,
	pub coefficients: Vec<T>
}
// New
impl<T:SeaCoefficient, I:SeaBasisBladeSet> SeaLowDimensionalSparseMultivector<T, I> {
	pub fn new(basis_blade_set: I, coefficients: Vec<T>) -> Self {
		Self { 
			basis_blade_set, 
			coefficients
		}
	}
}

// size queries
impl<T:SeaCoefficient, I:SeaBasisBladeSet> SeaLowDimensionalSparseMultivector<T, I> {
	pub fn is_zero(&self) -> bool {
		self.basis_blade_set.is_empty()
	}
}

// simplify
impl<T:SeaCoefficient, I:SeaBasisBladeSet> SeaLowDimensionalSparseMultivector<T, I> {	
	pub fn remove_zeros(&mut self) {
		for (i, idx) in self.basis_blade_set.iter().enumerate(){
			if self.coefficients[i].is_zero() {
				self.basis_blade_set.toggle(&idx);		
			}
		}
		self.coefficients.retain(|x| !x.is_zero())
	}

	pub fn remove_approx_zeros(&mut self, epsilon: T) {
		for (i, idx) in self.basis_blade_set.iter().enumerate(){
			if self.coefficients[i].abs() < epsilon {
				self.basis_blade_set.toggle(&idx);
			}
		}

		self.coefficients.retain(|x| x.abs() >= epsilon)
	}

	pub fn round_to_decimals(&mut self, decimals: i32) {
		let sf = T::from(10.0).unwrap().powi(decimals);
		for c in &mut self.coefficients {
			*c = (*c * sf).round() / sf;
		}
	}
}

// Multivector Index
impl<T:SeaCoefficient, I:SeaBasisBladeSet> SeaLowDimensionalSparseMultivector<T, I> {	
	pub fn basis_blade_coefficient(&self, index: &SeaBasisBladeIndexLow) -> T {
		if !self.basis_blade_set.has(index) {
			return T::zero();
		}
		let coeff_idx = self.basis_blade_set.count_basis_blades_less_than(index);
		self.coefficients[coeff_idx as usize]
	}
}

// Multivector Product Helpers
impl<T:SeaCoefficient, I:SeaBasisBladeSet> SeaLowDimensionalSparseMultivector<T, I> {	
	fn metric_product_constants(metric: &SeaMetric) -> (u8, SeaBasisBladeIndexLow, SeaBasisBladeIndexLow) {
		let num_metric_non_zeros = metric.positive + metric.negitive;
		let metric_non_zeros = SeaBasisBladeIndexLow::basis_vectors_less_than_or_equal(num_metric_non_zeros);
		let metric_positive = SeaBasisBladeIndexLow::basis_vectors_less_than_or_equal(metric.positive);
		let metric_negs = &metric_non_zeros ^ &metric_positive;
		(num_metric_non_zeros, metric_non_zeros, metric_negs)
	}
	
	/// dot two multivectors with their coefficients together
	/// based on a product mask and a neg mask telling what to negate
	/// allong with the current output coefficent index
	fn special_dot_product(&self, rhs: &Self, dot_mask: &I, neg_mask: &I, output_index: &SeaBasisBladeIndexLow) -> T {
		let mut out = T::zero();
		for i in dot_mask.iter() {
			let a_idx = self.basis_blade_set.count_basis_blades_less_than(&i) as usize;
			let b_look_at_idx_offset = &i ^ output_index;
			let b_idx = rhs.basis_blade_set.count_basis_blades_less_than(&b_look_at_idx_offset) as usize;
			let c_a = self.coefficients[a_idx];
			let c_b = rhs.coefficients[b_idx];
			let mut c = c_a * c_b;
			
			if neg_mask.has(&i){
				c = -c; 
			}

			out += c;
		}

		out
	}
	
	#[allow(clippy::uninit_vec)]
	unsafe fn create_coefficients_array(output_basis_blade_set: &I) -> Vec<T> {
		let output_num_coeffs = output_basis_blade_set.count_basis_blades() as usize;
		let mut output_coefficients = Vec::<T>::with_capacity(output_num_coeffs);
		unsafe { output_coefficients.set_len(output_num_coeffs)}
		output_coefficients
	}

	fn is_scalar(&self) -> bool{
		return if self.basis_blade_set.is_scalar() {
			true
		}
		else if self.basis_blade_set.has(&SeaBasisBladeIndexLow::SCALAR){
			self.coefficients[1..].iter().all(|x| x == &T::zero())
		}
		else {
			false
		}
	}
}

// Multivector Iter
impl<T:SeaCoefficient, I:SeaBasisBladeSet> SeaLowDimensionalSparseMultivector<T, I> {
	pub fn iter(&self) -> SparseMultivectorIter<'_, T, I> {
		SparseMultivectorIter::new(self)
	}

	pub fn grade_iter<'a>(&'a self, grade: &SeaGrade) -> SparseMultivectorGradeIter<'a, T, I> {
		SparseMultivectorGradeIter::new(self, grade)
	}
}

impl<T:SeaCoefficient, I:SeaBasisBladeSet> Default for SeaLowDimensionalSparseMultivector<T, I> {
	fn default() -> Self {
		Self { 
			basis_blade_set: SeaBasisBladeSet::EMPTY, 
			coefficients: vec![]
		}
	}
}

impl<T:SeaCoefficient, I:SeaBasisBladeSet> FromStr for SeaLowDimensionalSparseMultivector<T, I> {
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
			if *bb.index as usize >= I::NUM_BASIS_BLADES {
				return Err(ParseMultivectorError::InvalidBasisBladeForAlgebra);
			}
			if basis_blades.iter().any(|i| i.index == bb.index) {
				return Err(ParseMultivectorError::DuplicateBasisBlade);
			}
			basis_blades.push(bb);
		}

		basis_blades.sort_by_key(|i| *i.index);

		let mut out_bbbs = I::EMPTY;

		for bb in &basis_blades {
			out_bbbs.union_mut(&I::from(&bb.index));
		}

		Ok(Self {
			basis_blade_set: out_bbbs,
			coefficients: basis_blades.iter().map(|i| i.coefficient).collect()
		})

	}
}

impl<T:SeaCoefficient, I:SeaBasisBladeSet> BitXor<Self> for &SeaLowDimensionalSparseMultivector<T, I> {
	type Output = SeaLowDimensionalSparseMultivector<T, I>;
	fn bitxor(self, rhs: Self) -> Self::Output {
		self.outer_product(rhs)
	}
}

impl<T:SeaCoefficient, I:SeaBasisBladeSet> Neg for &SeaLowDimensionalSparseMultivector<T, I> {
	type Output = SeaLowDimensionalSparseMultivector<T, I>;
	fn neg(self) -> Self::Output {
		self.negated()
	}
}

#[cfg(test)]
mod tests{
	use std::time::Instant;

	use crate::tests::catch_unwind_silent;

	use super::*;

	#[test]
	fn from_str() {
		let a = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123").unwrap();
		assert_eq!(
			a,
			SeaLowDimensionalSparseMultivector { 
				basis_blade_set: 255, 
				coefficients: vec![1_f32, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0] 
			}
		);

		let a = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("1 - e1 + 2e2 - 3e12 + 4e3 - 5e13 + 6e23 - 7e123 + 8e4 - 9e14 + 10e24 - 11e124 + 12e34 - 13e134 + 14e234 - 15e1234").unwrap();
		assert_eq!(
			a,
			SeaLowDimensionalSparseMultivector { 
				basis_blade_set: 65535, 
				coefficients: vec![1_f32, -1.0, 2.0, -3.0, 4.0, -5.0, 6.0, -7.0, 8.0, -9.0, 10.0, -11.0, 12.0, -13.0, 14.0, -15.0] 
			}
		);

		assert_eq!(SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("e4").unwrap_err(), ParseMultivectorError::InvalidBasisBladeForAlgebra);
		assert_eq!(SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("e21").unwrap_err(), ParseMultivectorError::InvalidBasisVectorOrdering);
		assert_eq!(SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("e11").unwrap_err(), ParseMultivectorError::DuplicateBasisVector);
		assert_eq!(SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("e1 + e1").unwrap_err(), ParseMultivectorError::DuplicateBasisBlade);
	}

	
	#[test]
	fn multivector() {
		let a = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123").unwrap();
		let out = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("1 + e1 + 2e2 - 3e12 + 4e3 - 5e13 - 6e23 - 7e123").unwrap();
		assert_eq!(
			a.reversion_conjugated(),
			out
		);

		let a = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123").unwrap();
		let out = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("1 - e1 - 2e2 + 3e12 - 4e3 + 5e13 + 6e23 - 7e123").unwrap();
		assert_eq!(
			a.involution_conjugated(),
			out
		);

		let a = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123").unwrap();
		let out = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("1 - e1 - 2e2 - 3e12 - 4e3 - 5e13 - 6e23 + 7e123").unwrap();
		assert_eq!(
			a.clifford_conjugated(),
			out
		);

		let a = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123").unwrap();
		let out = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("7 + 6e1 + 5e2 + 4e12 + 3e3 + 2e13 + 1e23 + 1e123").unwrap();
		assert_eq!(
			a.poincare_dual(),
			out
		);

		let a = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123").unwrap();
		let out = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("7 + 6e1 - 5e2 + 4e12 + 3e3 - 2e13 + 1e23 + 1e123").unwrap();
		assert_eq!(
			a.right_complimented(),
			out
		);

		let a = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e34 + 13e134 + 14e234 + 15e1234").unwrap();
		let out = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("15 + 14e1 + 13e2 + 12e12 + 11e3 + 10e13 + 9e23 + 8e123 + 7e4 + 6e14 + 5e24 + 4e124 + 3e34 + 2e134 + 1e234 + 1e1234").unwrap();
		assert_eq!(
			a.poincare_dual(),
			out
		);

		let a = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e34 + 13e134 + 14e234 + 15e1234").unwrap();
		let out = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("15 + 14e1 - 13e2 + 12e12 + 11e3 - 10e13 + 9e23 + 8e123 - 7e4 + 6e14 - 5e24 - 4e124 + 3e34 + 2e134 - 1e234 + 1e1234").unwrap();
		assert_eq!(
			a.right_complimented(),
			out
		);

		let a = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e34 + 13e134 + 14e234 + 15e1234").unwrap();
		let out = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("15 - 14e1 + 13e2 + 12e12 - 11e3 - 10e13 + 9e23 - 8e123 + 7e4 + 6e14 - 5e24 + 4e124 + 3e34 - 2e134 + 1e234 + 1e1234").unwrap();
		assert_eq!(
			a.left_complimented(),
			out
		);
	}

	#[test]
	fn multivector_mut(){
		let mut a = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123").unwrap();
		let out = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("1 + e1 + 2e2 - 3e12 + 4e3 - 5e13 - 6e23 - 7e123").unwrap();
		a.reversion_conjugation();
		assert_eq!(
			a,
			out
		);

		let mut a = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123").unwrap();
		let out = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("1 - e1 - 2e2 + 3e12 - 4e3 + 5e13 + 6e23 - 7e123").unwrap();
		a.involution_conjugation();
		assert_eq!(
			a,
			out
		);

		let mut a = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123").unwrap();
		let out = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("1 - e1 - 2e2 - 3e12 - 4e3 - 5e13 - 6e23 + 7e123").unwrap();
		a.clifford_conjugation();
		assert_eq!(
			a,
			out
		);

		let mut a = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123").unwrap();
		let out = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("7 + 6e1 + 5e2 + 4e12 + 3e3 + 2e13 + 1e23 + 1e123").unwrap();
		a.poincare_dualization();
		assert_eq!(
			a,
			out
		);

		let mut a = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123").unwrap();
		let out = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("7 + 6e1 - 5e2 + 4e12 + 3e3 - 2e13 + 1e23 + 1e123").unwrap();
		a.right_compliment();
		assert_eq!(
			a,
			out
		);

		let mut a = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123").unwrap();
		let out = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("7 + 6e1 - 5e2 + 4e12 + 3e3 - 2e13 + 1e23 + 1e123").unwrap();
		a.left_compliment();
		assert_eq!(
			a,
			out
		);

		let mut a = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e34 + 13e134 + 14e234 + 15e1234").unwrap();
		let out = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("15 + 14e1 - 13e2 + 12e12 + 11e3 - 10e13 + 9e23 + 8e123 - 7e4 + 6e14 - 5e24 - 4e124 + 3e34 + 2e134 - 1e234 + 1e1234").unwrap();
		a.right_compliment();
		assert_eq!(
			a,
			out
		);

		let mut a = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e34 + 13e134 + 14e234 + 15e1234").unwrap();
		let out = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("15 - 14e1 + 13e2 + 12e12 - 11e3 - 10e13 + 9e23 - 8e123 + 7e4 + 6e14 - 5e24 + 4e124 + 3e34 - 2e134 + 1e234 + 1e1234").unwrap();
		a.left_compliment();
		assert_eq!(
			a,
			out
		);
	}
	
	#[test]
	fn iter() {
		let a = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e34 + 13e134 + 14e234 + 15e1234").unwrap();
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

		let a = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str(" e1 + 3e12 + 4e3 + 13e134 + 15e1234").unwrap();
		let mut iter = a.iter();
		assert_eq!(iter.next(), Some(SeaBasisBladeLow::from_str("1e1").unwrap()));
		assert_eq!(iter.next(), Some(SeaBasisBladeLow::from_str("3e12").unwrap()));
		assert_eq!(iter.next(), Some(SeaBasisBladeLow::from_str("4e3").unwrap()));
		assert_eq!(iter.next(), Some(SeaBasisBladeLow::from_str("13e134").unwrap()));
		assert_eq!(iter.next(), Some(SeaBasisBladeLow::from_str("15e1234").unwrap()));
		assert_eq!(iter.next(), None);
	}

	#[test]
	fn grade_iter() {
		let a = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e34 + 13e134 + 14e234 + 15e1234").unwrap();
		
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
	fn remove_zeros() {
		let mut a = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("1 + e1 + 0e2 + 0e12 + 0e3 + 0e13 + 0e23 + 0e123 + 0e4 + 0e14 + 0e24 + 0e124 + 0e34 + 0e134 + 0e234 + 0e1234").unwrap();
		let out = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("1 + e1").unwrap();
		a.remove_zeros();
		assert_eq!(
			a,
			out
		);

		let mut a = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("0 + 0e1 + 0e2 + 0e12 + 0e3 + 1e13 + 0e23 + 0e123 + 0e4 + 0e14 + 0e24 + 0e124 + 10e34 + 0e134 + 0e234 + 0e1234").unwrap();
		let out = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("1e13 + 10e34").unwrap();
		a.remove_zeros();
		assert_eq!(
			a,
			out
		)
	}

	#[test]
	fn remove_approx_zeros() {
		let mut a = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("1 + e1 + 0.0000000000000000001e2 + 0.0000000000000000001e12 + 0.0000000000000000001e3 + 0.0000000000000000001e13 + 0.0000000000000000001e23 + 0.0000000000000000001e123 + 0.0000000000000000001e4 + 0.0000000000000000001e14 + 0.0000000000000000001e24 + 0.0000000000000000001e124 + 0.0000000000000000001e34 + 0.0000000000000000001e134 + 0.0000000000000000001e234 + 0.0000000000000000001e1234").unwrap();
		let out = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("1 + e1").unwrap();
		a.remove_approx_zeros(1e-15);
		assert_eq!(
			a,
			out
		);

		let mut a = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("0.0000000000000000001 + 0.0000000000000000001e1 + 0.0000000000000000001e2 + 0.0000000000000000001e12 + 0.0000000000000000001e3 + 1e13 + 0.0000000000000000001e23 + 0.0000000000000000001e123 + 0.0000000000000000001e4 + 0.0000000000000000001e14 + 0.0000000000000000001e24 + 0.0000000000000000001e124 + 10.0000000000000000001e34 + 0.0000000000000000001e134 + 0.0000000000000000001e234 + 0.0000000000000000001e1234").unwrap();
		let out = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("1e13 + 10e34").unwrap();
		a.remove_approx_zeros(1e-15);
		assert_eq!(
			a,
			out
		)
	}
}