use std::cmp::Ordering;
use std::iter::repeat;
use std::ops::{Add, AddAssign, BitXor, Neg, Index, Mul, MulAssign, Sub, Div, DivAssign, SubAssign};
use std::str::FromStr;

use lazy_static::lazy_static;
use num_traits::{Zero, ToPrimitive};
use regex::Regex;

use crate::coefficient::SeaCoefficient;
use crate::errors::ParseMultivectorError;
use crate::metric::SeaMetric;
use crate::multivector_traits::{Multivector, MultivectorMut};
use crate::multivector_product_traits::{LeftContraction, RightContraction, InnerProduct, OuterProduct, GeometricProduct, RegressiveProduct, ScalarProduct};
use crate::grade::SeaGrade;

use super::basis_blade_bits::{SeaBasisBladeBitsInner, SeaBasisBladeBitsIter};
use super::basis_blade_bits::SeaBasisBladeBits;
use super::basis_blade::SeaBasisBladeLow;
use super::basis_blade_index::SeaBasisBladeIndexLow;

mod scalar_operations;
mod basis_blade_operations;
mod linear_operations;
mod non_geometric_linear_operations;
mod non_linear_operations;

pub struct SparseMultivectorIter<'a, T:SeaCoefficient, I:SeaBasisBladeBitsInner>{
	multivector: &'a SeaMultivectorLowSparseCoefficents<T,I>,
	curr_iter: SeaBasisBladeBitsIter<I>
}
impl<'a, T:SeaCoefficient, I:SeaBasisBladeBitsInner> SparseMultivectorIter<'a, T, I> {
	pub fn new(multivector: &'a SeaMultivectorLowSparseCoefficents<T,I>) -> Self{
		SparseMultivectorIter {
			multivector,
			curr_iter: multivector.basis_blade_bits.into_iter()
		}
	}
}
impl<'a, T:SeaCoefficient, I:SeaBasisBladeBitsInner> Iterator for SparseMultivectorIter<'a, T, I> {
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

pub struct SparseMultivectorGradeIter<'a, T:SeaCoefficient, I:SeaBasisBladeBitsInner>{
	multivector: &'a SeaMultivectorLowSparseCoefficents<T,I>,
	curr_iter: SeaBasisBladeBitsIter<I>
}
impl<'a, T:SeaCoefficient, I:SeaBasisBladeBitsInner> SparseMultivectorGradeIter<'a, T, I> {
	pub fn new(multivector: &'a SeaMultivectorLowSparseCoefficents<T,I>, grade: &SeaGrade) -> Self{
		let iter_bits =  &multivector.basis_blade_bits & &SeaBasisBladeBits::<I>::GRADES[**grade as usize];
		SparseMultivectorGradeIter {
			multivector,
			curr_iter: iter_bits.into_iter()
		}
	}
}
impl<'a, T:SeaCoefficient, I:SeaBasisBladeBitsInner> Iterator for SparseMultivectorGradeIter<'a, T, I> {
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
pub struct SeaMultivectorLowSparseCoefficents<T:SeaCoefficient, I:SeaBasisBladeBitsInner>{
	pub basis_blade_bits: SeaBasisBladeBits<I>,
	pub coefficients: Vec<T>
}
// New
impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> SeaMultivectorLowSparseCoefficents<T, I> {
	pub fn new() -> Self {
		Self { 
			basis_blade_bits: SeaBasisBladeBits::ZERO, 
			coefficients: vec![]
		}
	}
}

// simplify
impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> SeaMultivectorLowSparseCoefficents<T, I> {	
	pub fn simplify(&mut self) {
		for (i, idx) in self.basis_blade_bits.into_iter().enumerate(){
			if self.coefficients[i].is_zero() {
				self.basis_blade_bits.toggle(&idx);
			}
		}

		self.coefficients.retain(|x| !x.is_zero());
	}
}

// Multivector Index
impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> SeaMultivectorLowSparseCoefficents<T, I> {	
	pub fn basis_blade_coefficient(&self, index: &SeaBasisBladeIndexLow) -> T {
		if !self.basis_blade_bits.has(index) {
			return T::zero();
		}
		let coeff_idx = self.basis_blade_bits.count_basis_blades_less_than(index);
		self.coefficients[coeff_idx as usize]
	}
}

// Multivector Product Helpers
impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> SeaMultivectorLowSparseCoefficents<T, I> {	
	fn metric_product_constants(metric: &SeaMetric) -> (u8, SeaBasisBladeIndexLow, SeaBasisBladeIndexLow) {
		let num_metric_non_zeros = metric.positive + metric.negitive;
		let metric_non_zeros = SeaBasisBladeIndexLow::basis_vectors_less_than_or_equal(num_metric_non_zeros);
		let metric_positive = SeaBasisBladeIndexLow::basis_vectors_less_than_or_equal(metric.positive);
		let metric_negs = &metric_non_zeros ^ &metric_positive;
		(num_metric_non_zeros, metric_non_zeros, metric_negs)
	}

	fn update_zero_cache(zero_cache: &mut (SeaBasisBladeIndexLow, SeaBasisBladeBits<I>), output_index: &SeaBasisBladeIndexLow, metric_non_zeros: &SeaBasisBladeIndexLow) {
		let metric_zeros_lookup = output_index | metric_non_zeros;
		if zero_cache.0 != metric_zeros_lookup {
			let row = SeaBasisBladeBits::<I>::zero_mask_for_output_index(&metric_zeros_lookup);
			zero_cache.0 = metric_zeros_lookup;
			zero_cache.1 = row;
		}
	}
	
	fn total_neg_mask(output_index: &SeaBasisBladeIndexLow, metric_negs:&SeaBasisBladeIndexLow) -> SeaBasisBladeBits<I> {
		let cannonical_neg_mask = SeaBasisBladeBits::<I>::cannonical_reordering_neg_mask(output_index);
		let metric_neg_mask = SeaBasisBladeBits::<I>::metric_neg_mask(output_index, metric_negs);
		&cannonical_neg_mask ^ &metric_neg_mask
	}

	/// dot two multivectors with their coefficients together
	/// based on a product mask and a neg mask telling what to negate
	/// allong with the current output coefficent index
	fn special_dot_product(&self, rhs: &Self, dot_mask: &SeaBasisBladeBits<I>, neg_mask: &SeaBasisBladeBits<I>, output_index: &SeaBasisBladeIndexLow) -> T {
		let mut out = T::zero();
		for i in dot_mask {
			let a_idx = self.basis_blade_bits.count_basis_blades_less_than(&i) as usize;
			let b_look_at_idx_offset = &i ^ output_index;
			let b_idx = rhs.basis_blade_bits.count_basis_blades_less_than(&b_look_at_idx_offset) as usize;
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
	unsafe fn create_coefficients_array(output_basis_blade_bits: &SeaBasisBladeBits<I>) -> Vec<T> {
		let output_num_coeffs = output_basis_blade_bits.count_ones() as usize;
		let mut output_coefficients = Vec::<T>::with_capacity(output_num_coeffs);
		unsafe { output_coefficients.set_len(output_num_coeffs)}
		output_coefficients
	}

	fn geometric_product_swap(&mut self, swap_by:&SeaBasisBladeIndexLow){
		for i in swap_by{
			let bv_idx = i+1;
			
			let lookup_idx = bv_idx as usize;
			let bitmask = &SeaBasisBladeBits::<I>::CONTAINS_BASIS_VECS[lookup_idx];
			let mut msb = &self.basis_blade_bits & bitmask;
			let mut lsb = &self.basis_blade_bits & &!bitmask;
			
			let shift_amnt = 1_u8 << (bv_idx-1);
			let window = 1_u8 << bv_idx;
			let mask = SeaBasisBladeBits::<I>::basis_blades_less_than(&SeaBasisBladeIndexLow(window));
			
			let mut num_seen = 0;
			for i in (0..I::NUM_BASIS_BLADES).step_by(window as usize) {
				let t_mask = &mask << i;
				let msb_ones = (&msb & &t_mask).count_ones() as usize;
				let lsb_ones = (&lsb & &t_mask).count_ones() as usize;
				if msb_ones > 0 && lsb_ones > 0 {
					let (left, right) = self.coefficients.split_at_mut(num_seen + lsb_ones);
					left[num_seen..].reverse();
					right[..msb_ones].reverse();
					self.coefficients[num_seen..num_seen+lsb_ones+msb_ones].reverse();
				}
				num_seen += msb_ones + lsb_ones;
			}
		
			msb >>= shift_amnt;
			lsb <<= shift_amnt;
			self.basis_blade_bits = &msb | &lsb;
		}
	}

	fn is_scalar(&self) -> bool{
		return if self.basis_blade_bits == SeaBasisBladeBits::<I>::ONE{
			true
		}
		else if self.basis_blade_bits.has(&SeaBasisBladeIndexLow::ZERO){
			self.coefficients[1..].iter().all(|x| x == &T::zero())
		}
		else {
			false
		}
	}
}

// Multivector Iter
impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> SeaMultivectorLowSparseCoefficents<T, I> {
	pub fn iter(&self) -> SparseMultivectorIter<'_, T, I> {
		SparseMultivectorIter::new(self)
	}

	pub fn grade_iter<'a>(&'a self, grade: &SeaGrade) -> SparseMultivectorGradeIter<'a, T, I> {
		SparseMultivectorGradeIter::new(self, grade)
	}
}

impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> Default for SeaMultivectorLowSparseCoefficents<T, I> {
	fn default() -> Self {
		Self::new()
	}
}

impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> FromStr for SeaMultivectorLowSparseCoefficents<T, I> {
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

		let mut out_bbbs = SeaBasisBladeBits::<I>::ZERO;

		for bb in &basis_blades {
			out_bbbs |= &bb.index;
		}

		Ok(Self {
			basis_blade_bits: out_bbbs,
			coefficients: basis_blades.iter().map(|i| i.coefficient).collect()
		})

	}
}

impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> BitXor<Self> for &SeaMultivectorLowSparseCoefficents<T, I> {
	type Output = SeaMultivectorLowSparseCoefficents<T, I>;
	fn bitxor(self, rhs: Self) -> Self::Output {
		self.outer_product(rhs)
	}
}

impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> Neg for &SeaMultivectorLowSparseCoefficents<T, I> {
	type Output = SeaMultivectorLowSparseCoefficents<T, I>;
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
		let a = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123").unwrap();
		assert_eq!(
			a,
			SeaMultivectorLowSparseCoefficents { 
				basis_blade_bits: SeaBasisBladeBits(255), 
				coefficients: vec![1_f32, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0] 
			}
		);

		let a = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("1 - e1 + 2e2 - 3e12 + 4e3 - 5e13 + 6e23 - 7e123 + 8e4 - 9e14 + 10e24 - 11e124 + 12e34 - 13e134 + 14e234 - 15e1234").unwrap();
		assert_eq!(
			a,
			SeaMultivectorLowSparseCoefficents { 
				basis_blade_bits: SeaBasisBladeBits(65535), 
				coefficients: vec![1_f32, -1.0, 2.0, -3.0, 4.0, -5.0, 6.0, -7.0, 8.0, -9.0, 10.0, -11.0, 12.0, -13.0, 14.0, -15.0] 
			}
		);

		assert_eq!(SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("e4").unwrap_err(), ParseMultivectorError::InvalidBasisBladeForAlgebra);
		assert_eq!(SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("e21").unwrap_err(), ParseMultivectorError::InvalidBasisVectorOrdering);
		assert_eq!(SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("e11").unwrap_err(), ParseMultivectorError::DuplicateBasisVector);
		assert_eq!(SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("e1 + e1").unwrap_err(), ParseMultivectorError::DuplicateBasisBlade);
	}

	
	#[test]
	fn multivector() {
		let a = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123").unwrap();
		let out = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("1 + e1 + 2e2 - 3e12 + 4e3 - 5e13 - 6e23 - 7e123").unwrap();
		assert_eq!(
			a.reversion_conjugated(),
			out
		);

		let a = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123").unwrap();
		let out = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("1 - e1 - 2e2 + 3e12 - 4e3 + 5e13 + 6e23 - 7e123").unwrap();
		assert_eq!(
			a.involution_conjugated(),
			out
		);

		let a = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123").unwrap();
		let out = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("1 - e1 - 2e2 - 3e12 - 4e3 - 5e13 - 6e23 + 7e123").unwrap();
		assert_eq!(
			a.clifford_conjugated(),
			out
		);

		let a = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123").unwrap();
		let out = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("7 + 6e1 + 5e2 + 4e12 + 3e3 + 2e13 + 1e23 + 1e123").unwrap();
		assert_eq!(
			a.poincare_dual(),
			out
		);

		let a = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123").unwrap();
		let out = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("7 + 6e1 - 5e2 + 4e12 + 3e3 - 2e13 + 1e23 + 1e123").unwrap();
		assert_eq!(
			a.right_complimented(),
			out
		);

		let a = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e34 + 13e134 + 14e234 + 15e1234").unwrap();
		let out = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("15 + 14e1 + 13e2 + 12e12 + 11e3 + 10e13 + 9e23 + 8e123 + 7e4 + 6e14 + 5e24 + 4e124 + 3e34 + 2e134 + 1e234 + 1e1234").unwrap();
		assert_eq!(
			a.poincare_dual(),
			out
		);

		let a = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e34 + 13e134 + 14e234 + 15e1234").unwrap();
		let out = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("15 + 14e1 - 13e2 + 12e12 + 11e3 - 10e13 + 9e23 + 8e123 - 7e4 + 6e14 - 5e24 - 4e124 + 3e34 + 2e134 - 1e234 + 1e1234").unwrap();
		assert_eq!(
			a.right_complimented(),
			out
		);

		let a = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e34 + 13e134 + 14e234 + 15e1234").unwrap();
		let out = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("15 - 14e1 + 13e2 + 12e12 - 11e3 - 10e13 + 9e23 - 8e123 + 7e4 + 6e14 - 5e24 + 4e124 + 3e34 - 2e134 + 1e234 + 1e1234").unwrap();
		assert_eq!(
			a.left_complimented(),
			out
		);
	}

	#[test]
	fn multivector_mut(){
		let mut a = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123").unwrap();
		let out = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("1 + e1 + 2e2 - 3e12 + 4e3 - 5e13 - 6e23 - 7e123").unwrap();
		a.reversion_conjugation();
		assert_eq!(
			a,
			out
		);

		let mut a = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123").unwrap();
		let out = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("1 - e1 - 2e2 + 3e12 - 4e3 + 5e13 + 6e23 - 7e123").unwrap();
		a.involution_conjugation();
		assert_eq!(
			a,
			out
		);

		let mut a = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123").unwrap();
		let out = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("1 - e1 - 2e2 - 3e12 - 4e3 - 5e13 - 6e23 + 7e123").unwrap();
		a.clifford_conjugation();
		assert_eq!(
			a,
			out
		);

		let mut a = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123").unwrap();
		let out = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("7 + 6e1 + 5e2 + 4e12 + 3e3 + 2e13 + 1e23 + 1e123").unwrap();
		a.poincare_dualization();
		assert_eq!(
			a,
			out
		);

		let mut a = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123").unwrap();
		let out = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("7 + 6e1 - 5e2 + 4e12 + 3e3 - 2e13 + 1e23 + 1e123").unwrap();
		a.right_compliment();
		assert_eq!(
			a,
			out
		);

		let mut a = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123").unwrap();
		let out = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("7 + 6e1 - 5e2 + 4e12 + 3e3 - 2e13 + 1e23 + 1e123").unwrap();
		a.left_compliment();
		assert_eq!(
			a,
			out
		);

		let mut a = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e34 + 13e134 + 14e234 + 15e1234").unwrap();
		let out = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("15 + 14e1 - 13e2 + 12e12 + 11e3 - 10e13 + 9e23 + 8e123 - 7e4 + 6e14 - 5e24 - 4e124 + 3e34 + 2e134 - 1e234 + 1e1234").unwrap();
		a.right_compliment();
		assert_eq!(
			a,
			out
		);

		let mut a = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e34 + 13e134 + 14e234 + 15e1234").unwrap();
		let out = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("15 - 14e1 + 13e2 + 12e12 - 11e3 - 10e13 + 9e23 - 8e123 + 7e4 + 6e14 - 5e24 + 4e124 + 3e34 - 2e134 + 1e234 + 1e1234").unwrap();
		a.left_compliment();
		assert_eq!(
			a,
			out
		);
	}
	
	#[test]
	fn iter() {
		let a = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e34 + 13e134 + 14e234 + 15e1234").unwrap();
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

		let a = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str(" e1 + 3e12 + 4e3 + 13e134 + 15e1234").unwrap();
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
		let a = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e34 + 13e134 + 14e234 + 15e1234").unwrap();
		
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
	fn simplify() {
		let mut a = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("1 + e1 + 0e2 + 0e12 + 0e3 + 0e13 + 0e23 + 0e123 + 0e4 + 0e14 + 0e24 + 0e124 + 0e34 + 0e134 + 0e234 + 0e1234").unwrap();
		let out = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("1 + e1").unwrap();
		a.simplify();
		assert_eq!(
			a,
			out
		);

		let mut a = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("0 + 0e1 + 0e2 + 0e12 + 0e3 + 1e13 + 0e23 + 0e123 + 0e4 + 0e14 + 0e24 + 0e124 + 10e34 + 0e134 + 0e234 + 0e1234").unwrap();
		let out = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("1e13 + 10e34").unwrap();
		a.simplify();
		assert_eq!(
			a,
			out
		)
	}
}