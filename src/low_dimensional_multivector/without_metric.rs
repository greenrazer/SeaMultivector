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
		let coeff = self.multivector[&next_index];
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
		let coeff = self.multivector[&next_index];
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

	fn negated_basis_blades(&self, blade_bits: &SeaBasisBladeBits<I>) -> Self{
		let mut new_coeffs = self.coefficients.clone();
		let negation_bits = &self.basis_blade_bits & blade_bits;
		for idx in &negation_bits {
			let i = self.basis_blade_bits.count_basis_blades_less_than(&idx); 
			new_coeffs[i as usize] = -new_coeffs[i as usize];
		}

		Self { 
			basis_blade_bits: self.basis_blade_bits.clone(), 
			coefficients: new_coeffs
		}
	}

	fn negate_basis_blades(&mut self, blade_bits: &SeaBasisBladeBits<I>){
		let negation_bits = &self.basis_blade_bits & blade_bits;
		for idx in &negation_bits {
			let i = self.basis_blade_bits.count_basis_blades_less_than(&idx); 
			self.coefficients[i as usize] = -self.coefficients[i as usize];
		}
	}

	fn negate(&mut self) {
		for i in 0..self.coefficients.len() {
			self.coefficients[i] = -self.coefficients[i]
		}
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

// Dual
impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> SeaMultivectorLowSparseCoefficents<T, I> {
	pub fn dual(&self, metric: &SeaMetric) -> Self {
		let mut output = self.poincare_dual();
		output.dualization_helper(metric, false);
		output
	}

	pub fn undual(&self, metric: &SeaMetric) -> Self {
		let mut output = self.poincare_dual();
		output.dualization_helper(metric, true);
		output
	}

	pub fn dualization(&mut self, metric: &SeaMetric) {
		self.poincare_dualization();
		self.dualization_helper(metric, false)
	}

	pub fn undualization(&mut self, metric: &SeaMetric) {
		self.poincare_dualization();
		self.dualization_helper(metric, true)
	}

	fn dualization_helper(&mut self, metric: &SeaMetric, inverse: bool) {
		if metric.has_zero() {
			let negations = if inverse {
				SeaBasisBladeBits::<I>::LEFT_COMPLIMENT_NEGATIONS
			}
			else {
				SeaBasisBladeBits::<I>::RIGHT_COMPLIMENT_NEGATIONS
			};
			self.negate_basis_blades(&negations);
			return;
		}
		let mut negations = SeaBasisBladeBits::<I>::DUAL_NEGATIONS;

		if !inverse && metric.pseudoscalar_inverse_negitive() {
			negations = !negations
		}

		if metric.has_negitive() {
			let (_, _, metric_negs) = Self::metric_product_constants(metric);
			let metric_neg_mask = SeaBasisBladeBits::<I>::dual_metric_neg_mask(&metric_negs);
			negations ^= &metric_neg_mask;
		}

		self.negate_basis_blades(&negations);
	}
}

// Scalar Norm
impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> SeaMultivectorLowSparseCoefficents<T, I> {
	pub fn scalar_norm(&self, metric: &SeaMetric) -> T {
		let (_, metric_non_zeros, metric_negs) = Self::metric_product_constants(metric);

		
		let mut total_neg_mask = Self::total_neg_mask(&SeaBasisBladeIndexLow::ZERO, &metric_negs);
		total_neg_mask ^= &SeaBasisBladeBits::<I>::REVERSION_NEGATIONS;
		
		let zeros = SeaBasisBladeBits::<I>::zero_mask_for_output_index(&metric_non_zeros);
		let dot_mask = &self.basis_blade_bits & &zeros;

		let mut out = T::zero();
		for i in &dot_mask {
			let idx = self.basis_blade_bits.count_basis_blades_less_than(&i) as usize;
			let mut c = self.coefficients[idx];
			c *= c;
			out += if total_neg_mask.has(&i) {-c} else {c};
		}

		out.abs().sqrt()
	}

	pub fn scalar_normalized(&self, metric: &SeaMetric) -> Self {
		self / self.scalar_norm(metric)
	}

	pub fn scalar_normalize(&mut self, metric: &SeaMetric) {
		*self /= self.scalar_norm(metric);
	}
}

// Multivector
impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> SeaMultivectorLowSparseCoefficents<T, I>{
	pub fn reversion_conjugated(&self) -> Self {
		self.negated_basis_blades(&SeaBasisBladeBits::<I>::REVERSION_NEGATIONS)
	}

	pub fn involution_conjugated(&self) -> Self {
		self.negated_basis_blades(&SeaBasisBladeBits::<I>::ODD_GRADES)
	}

	pub fn clifford_conjugated(&self) -> Self {
		self.negated_basis_blades(&SeaBasisBladeBits::<I>::CONJUGATION_NEGATIONS)
	}
	
	pub fn poincare_dual(&self) -> Self {
		Self { 
			basis_blade_bits: self.basis_blade_bits.reverse_bits(), 
			coefficients: self.coefficients.iter().copied().rev().collect()
		}
	}

	pub fn right_complimented(&self) -> Self {
		let mut out = self.poincare_dual();
		out.negate_basis_blades(&SeaBasisBladeBits::<I>::RIGHT_COMPLIMENT_NEGATIONS);
		out
	}

	pub fn left_complimented(&self) -> Self {
		let mut out = self.poincare_dual();
		out.negate_basis_blades(&SeaBasisBladeBits::<I>::LEFT_COMPLIMENT_NEGATIONS);
		out
	}
}
// Multivector Mut
impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> SeaMultivectorLowSparseCoefficents<T, I>{
	pub fn reversion_conjugation(&mut self) {
		self.negate_basis_blades(&SeaBasisBladeBits::<I>::REVERSION_NEGATIONS)
	}

	pub fn involution_conjugation(&mut self) {
		self.negate_basis_blades(&SeaBasisBladeBits::<I>::ODD_GRADES)
	}

	pub fn clifford_conjugation(&mut self) {
		self.negate_basis_blades(&SeaBasisBladeBits::<I>::CONJUGATION_NEGATIONS)
	}

	pub fn poincare_dualization(&mut self){
		self.basis_blade_bits = self.basis_blade_bits.reverse_bits();
		self.coefficients.reverse();
	}

	pub fn right_compliment(&mut self) {
		self.poincare_dualization();
		self.negate_basis_blades(&SeaBasisBladeBits::<I>::RIGHT_COMPLIMENT_NEGATIONS);
	}

	pub fn left_compliment(&mut self) {
		self.poincare_dualization();
		self.negate_basis_blades(&SeaBasisBladeBits::<I>::LEFT_COMPLIMENT_NEGATIONS);
	}
}

// Left Contraction
impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> SeaMultivectorLowSparseCoefficents<T, I> {
	pub fn left_contraction(&self, rhs: &Self, metric: &SeaMetric) -> Self {
		let (num_metric_non_zeros, metric_non_zeros, metric_negs) = Self::metric_product_constants(metric);
		
		let output_basis_blade_bits = {
			self.basis_blade_bits.left_contraction(
				(&rhs.basis_blade_bits, num_metric_non_zeros)
			)
		};

		let mut output_coefficients = unsafe {
			Self::create_coefficients_array(&output_basis_blade_bits)
		};

		let mut zero_cache: (SeaBasisBladeIndexLow, SeaBasisBladeBits<I>) = (SeaBasisBladeIndexLow::ZERO, SeaBasisBladeBits::<I>::ONE);
		for (out_idx, o) in output_basis_blade_bits.into_iter().enumerate() {
			Self::update_zero_cache(&mut zero_cache, &o, &metric_non_zeros);

			let total_neg_mask = Self::total_neg_mask(&o, &metric_negs);

			let mut curr_rhs_bbb = &rhs.basis_blade_bits >> &o;
			curr_rhs_bbb.remove_basis_blades_with_any(&o);

			let mut dot_mask = &self.basis_blade_bits & &curr_rhs_bbb;
			dot_mask &= &zero_cache.1;
			
			output_coefficients[out_idx] = Self::special_dot_product(
				self,
				rhs,
				&dot_mask,
				&total_neg_mask,
				&o
			);
		}

		Self {
			basis_blade_bits: output_basis_blade_bits,
			coefficients: output_coefficients
		}
	}
}

// Right Contraction
impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> SeaMultivectorLowSparseCoefficents<T, I> {
	pub fn right_contraction(&self, rhs: &Self, metric: &SeaMetric) -> Self {
		let (num_metric_non_zeros, metric_non_zeros, metric_negs) = Self::metric_product_constants(metric);

		let output_basis_blade_bits = {
			self.basis_blade_bits.right_contraction(
				(&rhs.basis_blade_bits, num_metric_non_zeros)
			)
		};

		let output_num_coeffs = output_basis_blade_bits.count_ones() as usize;
		let mut output_coefficients = vec![T::zero(); output_num_coeffs];

		let mut zero_cache: (SeaBasisBladeIndexLow, SeaBasisBladeBits<I>) = (SeaBasisBladeIndexLow::ZERO, SeaBasisBladeBits::<I>::ONE);
		for (out_idx, o) in output_basis_blade_bits.into_iter().enumerate() {
			Self::update_zero_cache(&mut zero_cache, &o, &metric_non_zeros);

			let total_neg_mask = Self::total_neg_mask(&o, &metric_negs);

			let mut curr_rhs_bbb = &rhs.basis_blade_bits << &o;
			curr_rhs_bbb.remove_basis_blades_without_all(&o);

			let mut dot_mask = &self.basis_blade_bits & &curr_rhs_bbb;
			dot_mask &= &zero_cache.1;

			output_coefficients[out_idx] = Self::special_dot_product(
				self,
				rhs,
				&dot_mask,
				&total_neg_mask,
				&o
			);
		}

		Self {
			basis_blade_bits: output_basis_blade_bits,
			coefficients: output_coefficients
		}
	}
}

// Inner Product
impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> SeaMultivectorLowSparseCoefficents<T, I> {
	pub fn inner_product(&self, rhs: &Self, metric: &SeaMetric) -> Self {
		let (num_metric_non_zeros, metric_non_zeros, metric_negs) = Self::metric_product_constants(metric);

		let output_basis_blade_bits = {
			self.basis_blade_bits.inner_product(
				(&rhs.basis_blade_bits, num_metric_non_zeros)
			)
		};

		let mut output_coefficients = unsafe {
			Self::create_coefficients_array(&output_basis_blade_bits)
		};

		let mut zero_cache: (SeaBasisBladeIndexLow, SeaBasisBladeBits<I>) = (SeaBasisBladeIndexLow::ZERO, SeaBasisBladeBits::<I>::ONE);
		for (out_idx, o) in output_basis_blade_bits.into_iter().enumerate() {
			Self::update_zero_cache(&mut zero_cache, &o, &metric_non_zeros);

			let total_neg_mask = Self::total_neg_mask(&o, &metric_negs);

			let mut left_contraction_rhs_bbb = &rhs.basis_blade_bits >> &o;
			left_contraction_rhs_bbb.remove_basis_blades_with_any(&o);

			let mut dot_mask_left = &self.basis_blade_bits & &left_contraction_rhs_bbb;
			dot_mask_left &= &zero_cache.1;

			output_coefficients[out_idx] = Self::special_dot_product(
				self,
				rhs,
				&dot_mask_left,
				&total_neg_mask,
				&o
			);

			if !o.is_zero() {
				let mut right_contraction_rhs_bbb = &rhs.basis_blade_bits << &o;
				right_contraction_rhs_bbb.remove_basis_blades_without_all(&o);

				let mut dot_mask_right = &self.basis_blade_bits & &right_contraction_rhs_bbb;
				dot_mask_right &= &zero_cache.1;

				output_coefficients[out_idx] += Self::special_dot_product(
					self,
					rhs,
					&dot_mask_right,
					&total_neg_mask,
					&o
				);
			}
		}

		Self {
			basis_blade_bits: output_basis_blade_bits,
			coefficients: output_coefficients
		}
	}
}

// Scalar Product
impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> SeaMultivectorLowSparseCoefficents<T, I> {
	pub fn scalar_product(&self, rhs: &Self, metric: &SeaMetric) -> T {
		let (_, metric_non_zeros, metric_negs) = Self::metric_product_constants(metric);
		
		let total_neg_mask = Self::total_neg_mask(&SeaBasisBladeIndexLow::ZERO, &metric_negs);
		
		let zeros = SeaBasisBladeBits::<I>::zero_mask_for_output_index(&metric_non_zeros);
		let dot_mask = &(&self.basis_blade_bits & &rhs.basis_blade_bits) & &zeros;

		let mut out = T::zero();
		for i in &dot_mask {
			let a_idx = self.basis_blade_bits.count_basis_blades_less_than(&i) as usize;
			let b_idx = rhs.basis_blade_bits.count_basis_blades_less_than(&i) as usize;
			let c_a = self.coefficients[a_idx];
			let c_b = rhs.coefficients[b_idx];
			let c = c_a*c_b;
			
			out += if total_neg_mask.has(&i) {-c} else {c};
		}

		out
	}
}

// Regressive Product
impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> SeaMultivectorLowSparseCoefficents<T, I> {
	pub fn regressive_product(&self, rhs: &Self, metric: &SeaMetric) -> Self {
		let (a_dual, b_dual) = (self.dual(metric), rhs.dual(metric));

		let mut out_dual = &a_dual ^ &b_dual;
		
		out_dual.undualization(metric);
		
		out_dual
	}
}

// Geometric Product
impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> SeaMultivectorLowSparseCoefficents<T, I>{ 
	pub fn geometric_product(&self, rhs: &Self, metric: &SeaMetric) -> Self {
		let (num_metric_non_zeros, metric_non_zeros, metric_negs) = Self::metric_product_constants(metric);

		let output_basis_blade_bits = {
			self.basis_blade_bits.geometric_product(
				(&rhs.basis_blade_bits, num_metric_non_zeros)
			)
		};

		let mut output_coefficients = unsafe {
			Self::create_coefficients_array(&output_basis_blade_bits)
		};

		let mut zero_cache: (SeaBasisBladeIndexLow, SeaBasisBladeBits<I>) = (SeaBasisBladeIndexLow::ZERO, SeaBasisBladeBits::<I>::ONE);
		for (out_idx, o) in output_basis_blade_bits.into_iter().enumerate() {
			Self::update_zero_cache(&mut zero_cache, &o, &metric_non_zeros);

			let total_neg_mask = Self::total_neg_mask(&o, &metric_negs);

			let curr_rhs_bbb = rhs.basis_blade_bits.geometric_product_swap(&o);

			let mut dot_mask = &self.basis_blade_bits & &curr_rhs_bbb;
			dot_mask &= &zero_cache.1;

			output_coefficients[out_idx] = Self::special_dot_product(
				self,
				rhs,
				&dot_mask,
				&total_neg_mask,
				&o
			);
		}

		Self {
			basis_blade_bits: output_basis_blade_bits,
			coefficients: output_coefficients
		}
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

impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> Index<&SeaBasisBladeIndexLow> for SeaMultivectorLowSparseCoefficents<T, I> {
	type Output = T;
	fn index(&self, index: &SeaBasisBladeIndexLow) -> &Self::Output {
		if !self.basis_blade_bits.has(index) {
			panic!("{:?} doesn't exist inside of multivector", index)
		}
		let coeff_idx = self.basis_blade_bits.count_basis_blades_less_than(index);
		&self.coefficients[coeff_idx as usize]
	}
}

impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> AddAssign<&SeaBasisBladeLow<T>> for SeaMultivectorLowSparseCoefficents<T, I> {
	#[allow(clippy::suspicious_op_assign_impl)]
	fn add_assign(&mut self, rhs: &SeaBasisBladeLow<T>) {
		let new_idx = self.basis_blade_bits.count_basis_blades_less_than(&rhs.index);
		if !self.basis_blade_bits.has(&rhs.index) {
			self.coefficients.insert(new_idx as usize, T::zero());
		}
		self.basis_blade_bits |= &rhs.index;
		self.coefficients[new_idx as usize] += rhs.coefficient;
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> Add<&SeaBasisBladeLow<T>> for &SeaMultivectorLowSparseCoefficents<T, I> {
	type Output = SeaMultivectorLowSparseCoefficents<T, I>;
	#[allow(clippy::suspicious_op_assign_impl)]
	fn add(self, rhs: &SeaBasisBladeLow<T>) -> Self::Output{
		let mut out = self.clone();
		out += rhs;
		out
	}
}

impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> AddAssign<&Self> for SeaMultivectorLowSparseCoefficents<T, I> {
	fn add_assign(&mut self, rhs: &Self) {
		let bb_out = &self.basis_blade_bits | &rhs.basis_blade_bits;

		let size_diff = bb_out.count_ones() as usize - self.coefficients.len();
		self.coefficients.splice(0..0, repeat(T::zero()).take(size_diff));

		let mut self_idx = 0;
		let mut rhs_idx = 0;
		for (i, idx) in bb_out.into_iter().enumerate() {
			if self.basis_blade_bits.has(&idx) {
				self.coefficients.swap(i, size_diff+self_idx);
				self_idx += 1;
			}
			
			if rhs.basis_blade_bits.has(&idx) {
				self.coefficients[i] += rhs.coefficients[rhs_idx];
				rhs_idx += 1;
			}
		}

		self.basis_blade_bits = bb_out;
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> Add<Self> for &SeaMultivectorLowSparseCoefficents<T, I> {
    type Output = SeaMultivectorLowSparseCoefficents<T, I>;

    fn add(self, rhs: Self) -> Self::Output {
		let mut out = self.clone();
		out += rhs;
		out
	}
}

impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> SubAssign<&Self> for SeaMultivectorLowSparseCoefficents<T, I> {
	fn sub_assign(&mut self, rhs: &Self) {
		let bb_out = &self.basis_blade_bits | &rhs.basis_blade_bits;

		let size_diff = bb_out.count_ones() as usize - self.coefficients.len();
		self.coefficients.splice(0..0, repeat(T::zero()).take(size_diff));

		let mut self_idx = 0;
		let mut rhs_idx = 0;
		for (i, idx) in bb_out.into_iter().enumerate() {
			if self.basis_blade_bits.has(&idx) {
				self.coefficients.swap(i, size_diff+self_idx);
				self_idx += 1;
			}
			
			if rhs.basis_blade_bits.has(&idx) {
				self.coefficients[i] -= rhs.coefficients[rhs_idx];
				rhs_idx += 1;
			}
		}

		self.basis_blade_bits = bb_out;
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> Sub<Self> for &SeaMultivectorLowSparseCoefficents<T, I> {
    type Output = SeaMultivectorLowSparseCoefficents<T, I>;

    fn sub(self, rhs: Self) -> Self::Output {
		let mut out = self.clone();
		out -= rhs;
		out
	}
}

impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> AddAssign<T> for SeaMultivectorLowSparseCoefficents<T, I> {
    fn add_assign(&mut self, rhs: T) {
		*self += &SeaBasisBladeLow {
			index: SeaBasisBladeIndexLow::ZERO,
			coefficient: rhs
		};
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> Add<T> for &SeaMultivectorLowSparseCoefficents<T, I> {
    type Output = SeaMultivectorLowSparseCoefficents<T, I>;

    fn add(self, rhs: T) -> Self::Output {
		let mut out = self.clone();
		out += &SeaBasisBladeLow {
			index: SeaBasisBladeIndexLow::ZERO,
			coefficient: rhs
		};
		out 
	}
}

impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> SubAssign<T> for SeaMultivectorLowSparseCoefficents<T, I> {
    fn sub_assign(&mut self, rhs: T) {
		*self += &SeaBasisBladeLow {
			index: SeaBasisBladeIndexLow::ZERO,
			coefficient: -rhs
		};
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> Sub<T> for &SeaMultivectorLowSparseCoefficents<T, I> {
    type Output = SeaMultivectorLowSparseCoefficents<T, I>;

    fn sub(self, rhs: T) -> Self::Output {
		let mut out = self.clone();
		out += &SeaBasisBladeLow {
			index: SeaBasisBladeIndexLow::ZERO,
			coefficient: -rhs
		};
		out 
	}
}

impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> MulAssign<T> for SeaMultivectorLowSparseCoefficents<T, I> {
	fn mul_assign(&mut self, rhs: T) {
		for i in self.coefficients.iter_mut() {
            *i *= rhs;
        }
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> Mul<T> for &SeaMultivectorLowSparseCoefficents<T, I> {
	type Output = SeaMultivectorLowSparseCoefficents<T, I>;
	fn mul(self, rhs: T) -> Self::Output {
		SeaMultivectorLowSparseCoefficents {
			basis_blade_bits: self.basis_blade_bits.clone(),
			coefficients: self.coefficients.iter().map(|x| *x*rhs).collect()
		}
	}
}

impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> DivAssign<T> for SeaMultivectorLowSparseCoefficents<T, I> {
	fn div_assign(&mut self, rhs: T) {
		for i in self.coefficients.iter_mut() {
            *i /= rhs;
        }
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> Div<T> for &SeaMultivectorLowSparseCoefficents<T, I> {
	type Output = SeaMultivectorLowSparseCoefficents<T, I>;
	fn div(self, rhs: T) -> Self::Output {
		SeaMultivectorLowSparseCoefficents {
			basis_blade_bits: self.basis_blade_bits.clone(),
			coefficients: self.coefficients.iter().map(|x| *x/rhs).collect()
		}
	}
}

impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> Neg for &SeaMultivectorLowSparseCoefficents<T, I> {
	type Output = SeaMultivectorLowSparseCoefficents<T, I>;

	fn neg(self) -> Self::Output {
		let mut out = self.clone();
		out.negate();
		out
	}
}

impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> OuterProduct<&Self> for SeaMultivectorLowSparseCoefficents<T, I> {
	type Output = Self;
	fn outer_product(&self, rhs: &Self) -> Self::Output {
		let output_basis_blade_bits = {
			self.basis_blade_bits.outer_product(&rhs.basis_blade_bits)
		};

		let mut output_coefficients = unsafe {
			Self::create_coefficients_array(&output_basis_blade_bits)
		};

		let reverse_rhs_bbbs = rhs.basis_blade_bits.reverse_bits();
		let initial_shift = &SeaBasisBladeBits::<I>::NUM_BASIS_BLADES - &SeaBasisBladeIndexLow::ONE;

		for (out_idx, o) in output_basis_blade_bits.into_iter().enumerate() {
			let cannonical_neg_mask = SeaBasisBladeBits::<I>::cannonical_reordering_neg_mask(&o);

			let mut curr_rhs_bbb = &reverse_rhs_bbbs >> &(&initial_shift - &o);
			curr_rhs_bbb &= &SeaBasisBladeBits::<I>::zero_mask_for_output_index(&o);

			let dot_mask = &self.basis_blade_bits & &curr_rhs_bbb;

			output_coefficients[out_idx] = Self::special_dot_product(
				self,
				rhs,
				&dot_mask,
				&cannonical_neg_mask,
				&o
			);
		}

		Self {
			basis_blade_bits: output_basis_blade_bits,
			coefficients: output_coefficients
		}
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> BitXor<Self> for &SeaMultivectorLowSparseCoefficents<T, I> {
	type Output = SeaMultivectorLowSparseCoefficents<T, I>;
	fn bitxor(self, rhs: Self) -> Self::Output {
		self.outer_product(rhs)
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
	fn add() {
		let mut a = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("-1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123").unwrap();
		let b = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("-1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123").unwrap();
		assert_eq!(
			&a + &b,
			SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("-2 + 2e1 + 4e2 + 6e12 + 8e3 + 10e13 + 12e23 + 14e123").unwrap()
		);

		a += &b;
		assert_eq!(
			a,
			SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("-2 + 2e1 + 4e2 + 6e12 + 8e3 + 10e13 + 12e23 + 14e123").unwrap()
		);

		let mut a = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("-1 + e1 + 2e2 + 3e12 + 5e13 + 6e23 + 7e123").unwrap();
		a += &SeaBasisBladeLow::from_str("4e3").unwrap();
		assert_eq!(
			a,
			SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("-1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123").unwrap()
		);

		let mut a = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("-1 + e1 + 2e2 + 3e12 + 5e13").unwrap();
		a += &SeaBasisBladeLow::from_str("4e3").unwrap();
		a += &SeaBasisBladeLow::from_str("7e123").unwrap();
		a += &SeaBasisBladeLow::from_str("6e23").unwrap();
		assert_eq!(
			a,
			SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("-1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123").unwrap()
		);
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
	
	mod scalar_norm {
		use super::*;

		#[test]
		fn works() {
			let m = SeaMetric::new(3, 0, 0);
			let a = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123").unwrap();
			assert_eq!(
				a.scalar_norm(&m),
				11.874342
			);
	
			let m = SeaMetric::new(3, 0, 0);
			let a = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123").unwrap();
			let out = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("0.084215194 + 0.084215194e1 + 0.16843039e2 + 0.25264558e12 + 0.33686078e3 + 0.42107597e13 + 0.50529116e23 + 0.5895064e123").unwrap();
			for (a_bb, o_bb) in a.scalar_normalized(&m).iter().zip(out.iter()) {
				assert_eq!(a_bb.index, o_bb.index);
				let abs_diff = (a_bb.coefficient - o_bb.coefficient).abs();
				assert!(abs_diff < 1e-7);
			}
			
			
			let m = SeaMetric::new(4, 0, 0);
			let a = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e34 + 13e134 + 14e234 + 15e1234").unwrap();
			let out = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("0.028386649 + 0.028386649e1 + 0.056773297e2 + 0.08515994e12 + 0.113546595e3 + 0.14193325e13 + 0.17031989e23 + 0.19870654e123 + 0.22709319e4 + 0.25547984e14 + 0.2838665e24 + 0.31225315e124 + 0.34063977e34 + 0.36902642e134 + 0.39741307e234 + 0.42579973e1234").unwrap();
			assert_eq!(
				a.scalar_normalized(&m),
				out
			);
	
			let m = SeaMetric::new(3, 0, 0);
			let a = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 1e123").unwrap();
			let out = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("0.10369517 +0.103695172e1 +0.20739034e2 +0.41478068e3 +0.31108552e12 +0.51847583e13 +0.62217104e23 +0.10369517e123").unwrap();
			for (a_bb, o_bb) in a.scalar_normalized(&m).iter().zip(out.iter()) {
				assert_eq!(a_bb.index, o_bb.index);
				let abs_diff = (a_bb.coefficient - o_bb.coefficient).abs();
				assert!(abs_diff < 1e-7);
			}

			let m = SeaMetric::new(3, 0, 1);
			let a = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e34 + 13e134 + 14e234 + 15e1234").unwrap();
			let out = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("0.084215194 +0.67372155e4 +0.084215194e1 +0.168430388e2 +0.3368607759e3 +0.7579367e14 +0.84215194e24 +1.0105823e34 +0.252645582e12 +0.4210759699e13 +0.5052911639e23 +0.92636716e124 +1.0947975e134 +1.1790127e234 +0.5895063281e123 +1.2632279e1234").unwrap();
			for (a_bb, o_bb) in a.scalar_normalized(&m).iter().zip(out.iter()) {
				assert_eq!(a_bb.index, o_bb.index);
				let abs_diff = (a_bb.coefficient - o_bb.coefficient).abs();
				assert!(abs_diff < 1e-7);
			}
		}
	}

	mod dual {
		use super::*;

		#[test]
		fn sanity_check() {
			let m = SeaMetric::new(3, 0, 0);
			let a = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("1e13").unwrap();
			let out_dual = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("1e123").unwrap();
			assert_eq!(
				&a ^ &a.dual(&m),
				out_dual
			);

			let m = SeaMetric::new(4, 0, 0);
			let a = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("1e13").unwrap();
			let out_dual = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("-1e1234").unwrap();
			assert_eq!(
				&a ^ &a.dual(&m),
				out_dual
			);

			let m = SeaMetric::new(5, 0, 0);
			let a = SeaMultivectorLowSparseCoefficents::<f32, u32>::from_str("1e13").unwrap();
			let out_dual = SeaMultivectorLowSparseCoefficents::<f32, u32>::from_str("-1e12345").unwrap();
			assert_eq!(
				&a ^ &a.dual(&m),
				out_dual
			);

			let m = SeaMetric::new(6, 0, 0);
			let a = SeaMultivectorLowSparseCoefficents::<f32, u64>::from_str("1e13").unwrap();
			let out_dual = SeaMultivectorLowSparseCoefficents::<f32, u64>::from_str("1e123456").unwrap();
			assert_eq!(
				&a ^ &a.dual(&m),
				out_dual
			);

			let m = SeaMetric::new(4, 0, 1);
			let a = SeaMultivectorLowSparseCoefficents::<f32, u32>::from_str("1e13").unwrap();
			let out_dual = SeaMultivectorLowSparseCoefficents::<f32, u32>::from_str("1e12345").unwrap();
			assert_eq!(
				&a ^ &a.dual(&m),
				out_dual
			);

			let m = SeaMetric::new(4, 0, 1);
			let a = SeaMultivectorLowSparseCoefficents::<f32, u32>::from_str("1e124").unwrap();
			let out_dual = SeaMultivectorLowSparseCoefficents::<f32, u32>::from_str("1e12345").unwrap();
			assert_eq!(
				&a ^ &a.dual(&m),
				out_dual
			);
		}

		#[test]
		fn positive_metric() {
			let m = SeaMetric::new(3, 0, 0);

			let a = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123").unwrap();
			let out_dual = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("7 + 6e1 - 5e2 - 4e12 + 3e3 + 2e13 - 1e23 - 1e123").unwrap();
			assert_eq!(
				a.dual(&m),
				out_dual
			);
			let out_undual = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("-7 - 6e1 + 5e2 + 4e12 - 3e3 - 2e13 + 1e23 + 1e123").unwrap();
			assert_eq!(
				a.undual(&m),
				out_undual
			);
			assert_eq!(
				a.undual(&m).dual(&m),
				a
			);
	
			let m = SeaMetric::new(4, 0, 0);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e34 + 13e134 + 14e234 + 15e1234").unwrap();
			let out_dual = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("15 +14e1 -13e2 +11e3 -7e4 -12e12 +10e13 -6e14 -9e23 +5e24 -3e34 -8e123 +4e124 -2e134 +e234 +e1234").unwrap();
			assert_eq!(
				a.dual(&m),
				out_dual
			);
			let out_undual = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("15 +14e1 -13e2 +11e3 -7e4 -12e12 +10e13 -6e14 -9e23 +5e24 -3e34 -8e123 +4e124 -2e134 +e234 +e1234").unwrap();
			assert_eq!(
				a.undual(&m),
				out_undual
			);
			assert_eq!(
				a.dual(&m).undual(&m),
				a
			);

			let m = SeaMetric::new(5, 0, 0);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u32>::from_str("1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e34 + 13e134 + 14e234 + 15e1234 + e5 + e15 + 2e25 + 3e125 + 4e35 + 5e135 + 6e235 + 7e1235 + 8e45 + 9e145 + 10e245 + 11e1245 + 12e345 + 13e1345 + 14e2345 + 15e12345").unwrap();
			let out_dual = SeaMultivectorLowSparseCoefficents::<f32, u32>::from_str("15 +14e1 -13e2 +11e3 -7e4 +15e5 -12e12 +10e13 -6e14 +14e15 -9e23 +5e24 -13e25 -3e34 +11e35 -7e45 -8e123 +4e124 -12e125 -2e134 +10e135 -6e145 +e234 -9e235 +5e245 -3e345 +e1234 -8e1235 +4e1245 -2e1345 +e2345 +e12345").unwrap();
			assert_eq!(
				a.dual(&m),
				out_dual
			);
			let out_undual = SeaMultivectorLowSparseCoefficents::<f32, u32>::from_str("15 +14e1 -13e2 +11e3 -7e4 +15e5 -12e12 +10e13 -6e14 +14e15 -9e23 +5e24 -13e25 -3e34 +11e35 -7e45 -8e123 +4e124 -12e125 -2e134 +10e135 -6e145 +e234 -9e235 +5e245 -3e345 +e1234 -8e1235 +4e1245 -2e1345 +e2345 +e12345").unwrap();
			assert_eq!(
				a.undual(&m),
				out_undual
			);
			assert_eq!(
				a.dual(&m).undual(&m),
				a
			);

			let m = SeaMetric::new(6, 0, 0);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u64>::from_str("1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e34 + 13e134 + 14e234 + 15e1234 + e5 + e15 + 2e25 + 3e125 + 4e35 + 5e135 + 6e235 + 7e1235 + 8e45 + 9e145 + 10e245 + 11e1245 + 12e345 + 13e1345 + 14e2345 + 15e12345 + e6 + e16 + 2e26 + 3e126 + 4e36 + 5e136 + 6e236 + 7e1236 + 8e46 + 9e146 + 10e246 + 11e1246 + 12e346 + 13e1346 + 14e2346 + 15e12346 + e56 + e156 + 2e256 + 3e1256 + 4e356 + 5e1356 + 6e2356 + 7e12356 + 8e456 + 9e1456 + 10e2456 + 11e12456 + 12e3456 + 13e13456 + 14e23456 + 15e123456").unwrap();
			let out_dual = SeaMultivectorLowSparseCoefficents::<f32, u64>::from_str("15 +14e1 -13e2 +11e3 -7e4 +15e5 -15e6 -12e12 +10e13 -6e14 +14e15 -14e16 -9e23 +5e24 -13e25 +13e26 -3e34 +11e35 -11e36 -7e45 +7e46 -15e56 -8e123 +4e124 -12e125 +12e126 -2e134 +10e135 -10e136 -6e145 +6e146 -14e156 +e234 -9e235 +9e236 +5e245 -5e246 +13e256 -3e345 +3e346 -11e356 +7e456 +e1234 -8e1235 +8e1236 +4e1245 -4e1246 +12e1256 -2e1345 +2e1346 -10e1356 +6e1456 +e2345 -e2346 +9e2356 -5e2456 +3e3456 +e12345 -e12346 +8e12356 -4e12456 +2e13456 -e23456 -e123456").unwrap();
			assert_eq!(
				a.dual(&m),
				out_dual
			);
			let out_undual = SeaMultivectorLowSparseCoefficents::<f32, u64>::from_str("-15 -14e1 +13e2 -11e3 +7e4 -15e5 +15e6 +12e12 -10e13 +6e14 -14e15 +14e16 +9e23 -5e24 +13e25 -13e26 +3e34 -11e35 +11e36 +7e45 -7e46 +15e56 +8e123 -4e124 +12e125 -12e126 +2e134 -10e135 +10e136 +6e145 -6e146 +14e156 -e234 +9e235 -9e236 -5e245 +5e246 -13e256 +3e345 -3e346 +11e356 -7e456 -e1234 +8e1235 -8e1236 -4e1245 +4e1246 -12e1256 +2e1345 -2e1346 +10e1356 -6e1456 -e2345 +e2346 -9e2356 +5e2456 -3e3456 -e12345 +e12346 -8e12356 +4e12456 -2e13456 +e23456 +e123456").unwrap();
			assert_eq!(
				a.undual(&m),
				out_undual
			);
			assert_eq!(
				a.dual(&m).undual(&m),
				a
			);
		}
		
		#[test]
		fn negitive_metric() {
			let m = SeaMetric::new(3, 1, 0);
			let a = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e34 + 13e134 + 14e234 + 15e1234").unwrap();
			let out_dual = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("15 +14e1 -13e2 +11e3 +7e4 -12e12 +10e13 +6e14 -9e23 -5e24 +3e34 -8e123 -4e124 +2e134 -e234 -e1234").unwrap();
			assert_eq!(
				a.dual(&m),
				out_dual
			);
			let out_undual = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str(" -15 -14e1 +13e2 -11e3 -7e4 +12e12 -10e13 -6e14 +9e23 +5e24 -3e34 +8e123 +4e124 -2e134 +e234 +e1234").unwrap();
			assert_eq!(
				a.undual(&m),
				out_undual
			);
			assert_eq!(
				a.undual(&m).dual(&m),
				a
			);
	
			let m = SeaMetric::new(2, 2, 0);
			let a = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e34 + 13e134 + 14e234 + 15e1234").unwrap();
			let out_dual = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("15 +14e1 -13e2 -11e3 +7e4 -12e12 -10e13 +6e14 +9e23 -5e24 -3e34 +8e123 -4e124 -2e134 +e234 +e1234").unwrap();
			//846 ns
			assert_eq!(
				a.dual(&m),
				out_dual
			);
			let out_undual = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("15 +14e1 -13e2 -11e3 +7e4 -12e12 -10e13 +6e14 +9e23 -5e24 -3e34 +8e123 -4e124 -2e134 +e234 +e1234").unwrap();
			//835 ns
			assert_eq!(
				a.undual(&m),
				out_undual
			);
			assert_eq!(
				a.dual(&m).undual(&m),
				a
			);
	
			let m = SeaMetric::new(2, 3, 0);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u32>::from_str("1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e34 + 13e134 + 14e234 + 15e1234 + e5 + e15 + 2e25 + 3e125 + 4e35 + 5e135 + 6e235 + 7e1235 + 8e45 + 9e145 + 10e245 + 11e1245 + 12e345 + 13e1345 + 14e2345 + 15e12345").unwrap();
			let out_dual = SeaMultivectorLowSparseCoefficents::<f32, u32>::from_str("15 +14e1 -13e2 -11e3 +7e4 -15e5 -12e12 -10e13 +6e14 -14e15 +9e23 -5e24 +13e25 -3e34 +11e35 -7e45 +8e123 -4e124 +12e125 -2e134 +10e135 -6e145 +e234 -9e235 +5e245 +3e345 +e1234 -8e1235 +4e1245 +2e1345 -e2345 -e12345").unwrap();
			assert_eq!(
				a.dual(&m),
				out_dual
			);
			let out_undual = SeaMultivectorLowSparseCoefficents::<f32, u32>::from_str(" -15 -14e1 +13e2 +11e3 -7e4 +15e5 +12e12 +10e13 -6e14 +14e15 -9e23 +5e24 -13e25 +3e34 -11e35 +7e45 -8e123 +4e124 -12e125 +2e134 -10e135 +6e145 -e234 +9e235 -5e245 -3e345 -e1234 +8e1235 -4e1245 -2e1345 +e2345 +e12345").unwrap();
			assert_eq!(
				a.undual(&m),
				out_undual
			);
			assert_eq!(
				a.dual(&m).undual(&m),
				a
			);
		}

		#[test]
		fn zero_metric() {
			let m = SeaMetric::new(2, 0, 1);

			let a = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123").unwrap();
			let out = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("7 +3e3 +6e1 -5e2 -2e13 +e23 +4e12 +e123").unwrap();
			assert_eq!(
				a.dual(&m),
				out
			);
			assert_eq!(
				a.dual(&m).undual(&m),
				a
			);
	
			let m = SeaMetric::new(1, 0, 2);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123").unwrap();
			let out = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("7 -5e2 +3e3 +6e1 +e23 +4e12 -2e13 +e123").unwrap();
			assert_eq!(
				a.dual(&m),
				out
			);
			assert_eq!(
				a.dual(&m).undual(&m),
				a
			);
	
			let m = SeaMetric::new(3, 0, 1);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e34 + 13e134 + 14e234 + 15e1234").unwrap();
			let out = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("15 + 14e1 - 13e2 + 12e12 + 11e3 - 10e13 + 9e23 + 8e123 - 7e4 + 6e14 - 5e24 - 4e124 + 3e34 + 2e134 - 1e234 + 1e1234").unwrap();
			assert_eq!(
				a.dual(&m),
				out
			);
			assert_eq!(
				a.dual(&m).undual(&m),
				a
			);
	
			let m = SeaMetric::new(4, 0, 1);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u32>::from_str("1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e34 + 13e134 + 14e234 + 15e1234 + e5 + e15 + 2e25 + 3e125 + 4e35 + 5e135 + 6e235 + 7e1235 + 8e45 + 9e145 + 10e245 + 11e1245 + 12e345 + 13e1345 + 14e2345 + 15e12345").unwrap();
			let out = SeaMultivectorLowSparseCoefficents::<f32, u32>::from_str("15 +15e5 +14e1 -13e2 +11e3 -7e4 -14e15 +13e25 -11e35 +7e45 +12e12 -10e13 +6e14 +9e23 -5e24 +3e34 +12e125 -10e135 +6e145 +9e235 -5e245 +3e345 +8e123 -4e124 +2e134 -e234 -8e1235 +4e1245 -2e1345 +e2345 +e1234 +e12345").unwrap();
			assert_eq!(
				a.dual(&m),
				out
			);
			assert_eq!(
				a.dual(&m).undual(&m),
				a
			);

			let m = SeaMetric::new(3, 0, 2);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u32>::from_str(" + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345").unwrap();
			let out = SeaMultivectorLowSparseCoefficents::<f32, u32>::from_str("15 +14e1 -13e2 +7e45 -11e35 +12e12 +6e145 -5e245 +3e345 -10e135 +9e235 +4e1245 -2e1345 +e2345 -8e1235 +e12345").unwrap();
			assert_eq!(
				a.dual(&m),
				out
			);
			assert_eq!(
				a.dual(&m).undual(&m),
				a
			);

			let a = SeaMultivectorLowSparseCoefficents::<f32, u32>::from_str("  + 3e13 + 1e4 + 2e124  +4e345").unwrap();
			let out = SeaMultivectorLowSparseCoefficents::<f32, u32>::from_str("+4e12 -2e35 -3e245 -e1235").unwrap();
			assert_eq!(
				a.dual(&m),
				out
			);
			assert_eq!(
				a.dual(&m).undual(&m),
				a
			);

			let a = SeaMultivectorLowSparseCoefficents::<f32, u32>::from_str(" + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345").unwrap();
			let out = SeaMultivectorLowSparseCoefficents::<f32, u32>::from_str("15 +14e1 -13e2 +7e45 -11e35 +12e12 +6e145 -5e245 +3e345 -10e135 +9e235 +4e1245 -2e1345 +e2345 -8e1235 +e12345").unwrap();
			assert_eq!(
				a.dual(&m),
				out
			);
			assert_eq!(
				a.dual(&m).undual(&m),
				a
			);

			let a = SeaMultivectorLowSparseCoefficents::<f32, u32>::from_str("  + 3e12 + 1e4 + 2e124  +4e345").unwrap();
			let out = SeaMultivectorLowSparseCoefficents::<f32, u32>::from_str("+4e12 -2e35 -e1235 + 3e345").unwrap();
			assert_eq!(
				a.dual(&m),
				out
			);
			assert_eq!(
				a.dual(&m).undual(&m),
				a
			);
		
		}
	}

	mod geometric_product {
		use super::*;

		#[test]
		fn euclidean() {
			let m = SeaMetric::new(3, 0, 0);
			
			let a = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("1").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("3e12").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("3e12").unwrap();
			//774ns
			assert_eq!(
				a.geometric_product(&b, &m),
				out
			);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("e1 + e2").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("e2 + e3").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("1 + e12 + e13 + e23").unwrap();
			//3.382µs
			assert_eq!(
				a.geometric_product(&b, &m),
				out
			);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("1 + e1 + e3").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("     e1  + 3e12").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("1 + e1 + 3e2 +3e12 -e13 + 3e123").unwrap();
			assert_eq!(
				a.geometric_product(&b, &m),
				out
			);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("1 + 2e1 + 1e12 + 2e3").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("  + 1e2 + 2e13 + 1e23").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("-3e1 - e2 + 2e12 + 4e3 + 3e13 -3e23 + 2e123").unwrap();
			//7.542µs
			assert_eq!(
				a.geometric_product(&b, &m),
				out
			);
	
			let m = SeaMetric::new(4, 0, 0);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str(" + 1 + e1 + 2e12 + e13 + e234 ").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("  + 1e4 + 2e124 + 3e13 +4e34").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("-3 -4e2 + 3e3 -3e4 +e13 +5e14 -5e23 + 2e24 + 4e34 + 7e124 + 5e134 + 2e234 + 8e1234").unwrap();
			//21.11µs
			//10.495µs
			assert_eq!(
				a.geometric_product(&b, &m),
				out
			);
		}
	
		#[test]
		fn projective() {
			let m = SeaMetric::new(3, 0, 1);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("1 + 2e1 + 1e12 + 2e3").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("  + 1e2 + 2e13 + 1e23").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("-3e1 - e2 + 2e12 + 4e3 + 3e13 -3e23 + 2e123").unwrap();
			assert_eq!(
				a.geometric_product(&b, &m),
				out
			);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str(" + 2e3 + 2e14 + 1e124 ").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("  + 1e24 + 2e134").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("-4e14 -2e234").unwrap();
			assert_eq!(
				a.geometric_product(&b, &m),
				out
			);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str(" + 1 + e1 + 2e12 + e13 + e234 ").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("  + 3e13 + 1e4 + 2e124  +4e34").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("-3 +3e3 + 3e13 -6e23 -3e4  + 5e14 +2e24 + 7e124 + 4e34  +5e134 + 2e234 + 8e1234").unwrap();
			//12.67µs
			assert_eq!(
				a.geometric_product(&b, &m),
				out
			);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("  + e1 + 2e2 - e234 + e1234 ").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("  - 3e13 + 1e4 - 2e124  +4e34").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("-3e3 +6e123 + 5e14 - 3e24 +3e124 + 4e134 + 8e234 ").unwrap();
			//8.403µs
			assert_eq!(
				a.geometric_product(&b, &m),
				out
			);
		}
	
		#[test]
		fn hyperbolic() {
			let m = SeaMetric::new(3, 1, 0);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("1 + 2e1 + 1e12 + 2e3").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("  + 1e2 + 2e13 + 1e23").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("-3e1 - e2 + 2e12 + 4e3 + 3e13 -3e23 + 2e123").unwrap();
			assert_eq!(
				a.geometric_product(&b, &m),
				out
			);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str(" + e1 + 2e2 - e234 + e1234  ").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("  - 3e13 + 1e4 - 2e124  +4e34").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("-4e2 -5e3 + 4e12 + 2e13 + 5e14 + e23 -3e24 + 5e123 + 3e124 + 4e134 + 8e234").unwrap();
			assert_eq!(
				a.geometric_product(&b, &m),
				out
			);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str(" + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e34 + 13e134 + 14e234 + 15e1234").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("  + 3e13 + 1e4 + 2e124  +4e34").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("47 +51e1 + 49e2 + 53e3 - 28e4 + 51e12 + 54e13 - 19e14 - 9e23 + 73e24 + 21e34 - e123 + 75e124 + 21e134 -9e234 -3e1234").unwrap();
			//21.547µs
			assert_eq!(
				a.geometric_product(&b, &m),
				out
			);
		}
	
		#[test]
		fn other() {
			let m = SeaMetric::new(3, 2, 0);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u32>::from_str(" + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u32>::from_str("  + 3e13 + 1e4 + 2e124  +4e345").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u32>::from_str("-49 -53e1 -63e2 +3e3 -5e4 -69e12 + 3e13 - 3e14 - 9e23 + 4e24 + 17e34 + 14e35 -23e45 -6e123 + 5e124 +17e134 + 21e135 - 16e145 -17e234 + 80e235 + 69e245 + 4e345 - 15e1234 +83e1235 + 70e1245 + 4e1345 + 8e2345 + 12e12345").unwrap();
			assert_eq!(
				a.geometric_product(&b, &m),
				out
			);
	
			let m = SeaMetric::new(3, 1, 1);
			
			let a = SeaMultivectorLowSparseCoefficents::<f32, u32>::from_str(" + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u32>::from_str("  + 3e13 + 1e4 + 2e124  +4e345").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u32>::from_str("-1 -e1 -7e2 +3e3 - 5e4 +14e35 -23e45 -9e12 +3e13 - 3e14 -9e23 +4e24 +17e34 +21e135 - 16e145 +80e235 +69e245 + 4e345 -6e123 + 5e124 + 17e134 - 17e234 +83e1235 +70e1245 +4e1345 + 8e2345 -15e1234 +12e12345").unwrap();
			// 45.752µs
			// new 15.492µs
			assert_eq!(
				a.geometric_product(&b, &m),
				out
			);
		}
	
		#[test]
		fn dimensions_6() {
			let m = SeaMetric::new(3, 3, 0);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u64>::from_str(" + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345 +e6 + 12e126 + 8e46").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u64>::from_str("  + 3e13 + 1e4 + 2e124  +4e345 +83e1235 + 70e1245 + 4e1345 + 8e2345 + 12e123456").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u64>::from_str("115+115e1-27e2+1053e3+1240e4+189e5+188e6-21e12+983e13+1159e14+202e15+168e16-919e23-1075e24-215e25-156e26+17e34-351e35-301e45+23e46-846e123-991e124-228e125-128e126+17e134-265e135+3e136-228e145-17e234+143e235-36e236+135e245+303e345-1160e356-924e456-15e1234+230e1235+124e1245-14e1246-560e1256+442e1345+24e1346-152e1356-72e1456-393e2345+44e2356+60e2456-40e3456-360e12345+179e12356+118e12456+76e13456-28e23456-700e123456").unwrap();
			// 171.677µs
			// new: 58.328µs
			// with unsafe: 56.79µs
			// 49.433µs
			assert_eq!(
				a.geometric_product(&b, &m),
				out
			);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u64>::from_str("1 + 2e1 + 1e12 + 2e3").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u64>::from_str("  + 1e2 + 2e13 + 1e23").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u64>::from_str("-3e1 - e2 + 2e12 + 4e3 + 3e13 -3e23 + 2e123").unwrap();
			//13.991µs
			//3.334µs
			assert_eq!(
				a.geometric_product(&b, &m),
				out
			);

			// let now = Instant::now();
			// let num_iters = 10_000;
			// let mut t = a.clone();
			// for _ in 0..num_iters {
			// 	t = t.geometric_product(&b, &m);
			// }
			// let elapsed_time = now.elapsed();
			// println!("time time time {:?}", elapsed_time/num_iters);
		}
	
	}

	mod left_contraction {
		use super::*;
		#[test]
		fn euclidean() {
			let m = SeaMetric::new(3, 0, 0);
			
			let a = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("1").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("3e12").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("3e12").unwrap();
			assert_eq!(
				a.left_contraction(&b, &m),
				out
			);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("e1 + e2").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("e2 + e3").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("1").unwrap();
			assert_eq!(
				a.left_contraction(&b, &m),
				out
			);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("1 + e1 + e3").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("     e1  + 3e12").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("1 + e1 + 3e2 +3e12").unwrap();
			assert_eq!(
				a.left_contraction(&b, &m),
				out
			);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("1 + 2e1 + 1e12 + 2e3").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("  + 1e2 + 2e13 + 1e23").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("-4e1 - e2 + 4e3 + 2e13 +1e23").unwrap();
			assert_eq!(
				a.left_contraction(&b, &m),
				out
			);
	
			let m = SeaMetric::new(4, 0, 0);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str(" + 1 + e1 + 2e12 + e13 + e234 ").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("  + 1e4 + 2e124 + 3e13 +4e34").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("-3 + 3e3 -3e4 + 3e13 + 2e24 + 4e34 + 2e124").unwrap();
			assert_eq!(
				a.left_contraction(&b, &m),
				out
			);
		}
	
		#[test]
		fn projective() {
			let m = SeaMetric::new(3, 0, 1);
			let a = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("  + e1 + 2e2 - e234 + e1234 ").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("  - 3e13 + 1e4 - 2e124  +4e34").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("-3e3 +4e14 -2e24").unwrap();
			assert_eq!(
				a.left_contraction(&b, &m),
				out
			);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str(" + 1 + e1 + 2e12 + e13 + e234 ").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("  + 3e13 + 1e4 + 2e124  +4e34").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("-3 -3e4 + 3e3 +2e24 + 4e34 + 3e13 + 2e124").unwrap();
			assert_eq!(
				a.left_contraction(&b, &m),
				out
			);
		}
	
		#[test]
		fn hyperbolic() {
			let m = SeaMetric::new(3, 1, 0);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("1 + 2e1 + 1e12 + 2e3").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("  + 1e2 + 2e13 + 1e23").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("-4e1 -e2 +4e3 +2e13 +e23").unwrap();
			assert_eq!(
				a.left_contraction(&b, &m),
				out
			);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str(" + e1 + 2e2 - e234 + e1234  ").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("  - 3e13 + 1e4 - 2e124  +4e34").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("-3e3 + 4e14 -2e24").unwrap();
			assert_eq!(
				a.left_contraction(&b, &m),
				out
			);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str(" + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e34 + 13e134 + 14e234 + 15e1234").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("  + 3e13 + 1e4 + 2e124  +4e34").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("47 +8e1 -18e2 +35e3 + 11e4 -16e12 + 3e13 -4e14 +2e24 +4e34 + 2e124").unwrap();
			//4.371µs
			assert_eq!(
				a.left_contraction(&b, &m),
				out
			);
		}
	
		#[test]
		fn other() {
			let m = SeaMetric::new(3, 2, 0);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u32>::from_str(" + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u32>::from_str("  + 3e13 + 1e4 + 2e124  +4e345").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u32>::from_str("-49 +8e1 -18e2 +3e3 -5e4 -16e12 +3e13 -4e14 +2e24 +32e35 +16e45 +2e124 + 4e345").unwrap();
			assert_eq!(
				a.left_contraction(&b, &m),
				out
			);
	
			let m = SeaMetric::new(3, 1, 1);
			
			let a = SeaMultivectorLowSparseCoefficents::<f32, u32>::from_str(" + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u32>::from_str("  + 3e13 + 1e4 + 2e124  +4e345").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u32>::from_str("-1 +8e1 -18e2 + 3e3 -5e4 +32e35 +16e45 -16e12 + 3e13 -4e14 +2e24 +4e345 +2e124").unwrap();
			//3.867µs
			assert_eq!(
				a.left_contraction(&b, &m),
				out
			);
		}
	
		#[test]
		fn dimensions_6() {
			let m = SeaMetric::new(3, 3, 0);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u64>::from_str(" + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345 +e6 + 12e126 + 8e46").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u64>::from_str("  + 3e13 + 1e4 + 2e124  +4e345 +83e1235 + 70e1245 + 4e1345 + 8e2345 + 12e123456").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u64>::from_str("115+56e1+78e2+3e3-5e4+189e5+180e6-16e12+3e13-4e14+202e15+168e16+2e24-215e25-156e26-333e35-262e45+2e124-228e125-144e126-198e135-156e145+19e235+38e245-120e345-132e356-84e456+179e1235+70e1245+4e1345-120e1356-72e1456+8e2345+108e2356+60e2456-36e3456+12e12345+96e12356+48e12456-24e13456+12e23456+12e123456").unwrap();
			//25.833µs
			//23.094µs
			assert_eq!(
				a.left_contraction(&b, &m),
				out
			);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u64>::from_str("1 + 2e1 + 1e12 + 2e3").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u64>::from_str("  + 1e2 + 2e13 + 1e23").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u64>::from_str("-4e1 -e2 +4e3 +2e13 +e23").unwrap();
			//2.081µs
			assert_eq!(
				a.left_contraction(&b, &m),
				out
			);

			// let now = Instant::now();
			// let num_iters = 10_000;
			// let mut t = a.clone();
			// for _ in 0..num_iters {
			// 	t = t.left_contraction(&b, &m);
			// }
			// let elapsed_time = now.elapsed();
			// println!("time time time {:?}", elapsed_time/num_iters);
		}
	}

	mod right_contraction {
		use super::*;
		#[test]
		fn euclidean() {
			let m = SeaMetric::new(3, 0, 0);
			
			let a = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("1").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("3e12").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("").unwrap();
			assert_eq!(
				a.right_contraction(&b, &m),
				out
			);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("e1 + e2").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("e2 + e3").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("1").unwrap();
			assert_eq!(
				a.right_contraction(&b, &m),
				out
			);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("1 + e1 + e3").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("     e1  + 3e12").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("1").unwrap();
			assert_eq!(
				a.right_contraction(&b, &m),
				out
			);

			let a = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("e12").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("e2").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("e1").unwrap();
			assert_eq!(
				a.right_contraction(&b, &m),
				out
			);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("1 + 2e1 + 1e12 + 2e3").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("  + 1e2 + 2e13 + 1e23").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("e1").unwrap();
			assert_eq!(
				a.right_contraction(&b, &m),
				out
			);
	
			let m = SeaMetric::new(4, 0, 0);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str(" + 1 + e1 + 2e12 + e13 + e234 ").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("  + 1e4 + 2e124 + 3e13 +4e34").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("-3 -4e2 + e23").unwrap();
			assert_eq!(
				a.right_contraction(&b, &m),
				out
			);
		}
	
		#[test]
		fn projective() {
			let m = SeaMetric::new(3, 0, 1);
			let a = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("  + e1 + 2e2 - e234 + e1234 ").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("  - 3e13 + 1e4 - 2e124  +4e34").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("-3e24").unwrap();
			assert_eq!(
				a.right_contraction(&b, &m),
				out
			);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str(" + 1 + e1 + 2e12 + e13 + e234 ").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("  + 3e13 + 1e4 + 2e124  +4e34").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("-3").unwrap();
			assert_eq!(
				a.right_contraction(&b, &m),
				out
			);
		}
	
		#[test]
		fn hyperbolic() {
			let m = SeaMetric::new(3, 1, 0);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("1 + 2e1 + 1e12 + 2e3").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("  + 1e2 + 2e13 + 1e23").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("e1").unwrap();
			assert_eq!(
				a.right_contraction(&b, &m),
				out
			);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str(" + e1 + 2e2 - e234 + e1234  ").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("  - 3e13 + 1e4 - 2e124  +4e34").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("-4e2 - 2e3 +4e12 + e23 - 3e24 - e123").unwrap();
			assert_eq!(
				a.right_contraction(&b, &m),
				out
			);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str(" + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e34 + 13e134 + 14e234 + 15e1234").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("  + 3e13 + 1e4 + 2e124  +4e34").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("47+43e1+67e2+18e3-39e4+49e12-13e13-14e23+45e24-15e123").unwrap();
			assert_eq!(
				a.right_contraction(&b, &m),
				out
			);
		}
	
		#[test]
		fn other() {
			let m = SeaMetric::new(3, 2, 0);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u32>::from_str(" + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u32>::from_str("  + 3e13 + 1e4 + 2e124  +4e345").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u32>::from_str("-49-61e1-45e2-71e12-18e35-39e45+13e135+14e235+45e245+15e1235").unwrap();
			assert_eq!(
				a.right_contraction(&b, &m),
				out
			);
	
			let m = SeaMetric::new(3, 1, 1);
			
			let a = SeaMultivectorLowSparseCoefficents::<f32, u32>::from_str(" + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u32>::from_str("  + 3e13 + 1e4 + 2e124  +4e345").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u32>::from_str("-1 -9e1 +11e2 -18e35 -39e45 -11e12 +13e135 +14e235 +45e245 +15e1235").unwrap();
			assert_eq!(
				a.right_contraction(&b, &m),
				out
			);
		}
	
		#[test]
		fn dimensions_6() {
			let m = SeaMetric::new(3, 3, 0);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u64>::from_str(" + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345 +e6 + 12e126 + 8e46").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u64>::from_str("  + 3e13 + 1e4 + 2e124  +4e345 +83e1235 + 70e1245 + 4e1345 + 8e2345 + 12e123456").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u64>::from_str("115+59e1-105e2+1050e3+1245e4+8e6-71e12-18e35-39e45+13e135+14e235+45e245+15e1235​").unwrap();
			assert_eq!(
				a.right_contraction(&b, &m),
				out
			);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u64>::from_str("1 + 2e1 + 1e12 + 2e3").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u64>::from_str("  + 1e2 + 2e13 + 1e23").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u64>::from_str("e1").unwrap();
			assert_eq!(
				a.right_contraction(&b, &m),
				out
			);
		}
	}

	mod inner_product {
		use super::*;
		#[test]
		fn euclidean() {
			let m = SeaMetric::new(3, 0, 0);
			
			let a = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("1").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("3e12").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("3e12").unwrap();
			assert_eq!(
				a.inner_product(&b, &m),
				out
			);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("e1 + e2").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("e2 + e3").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("1").unwrap();
			assert_eq!(
				a.inner_product(&b, &m),
				out
			);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("1 + e1 + e3").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("     e1  + 3e12").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("1 + e1 + 3e2 + 3e12").unwrap();
			assert_eq!(
				a.inner_product(&b, &m),
				out
			);

			let a = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("e12").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("e2").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("e1").unwrap();
			assert_eq!(
				a.inner_product(&b, &m),
				out
			);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("1 + 2e1 + 1e12 + 2e3").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("  + 1e2 + 2e13 + 1e23").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("-3e1 -e2 +4e3 +2e13 + e23").unwrap();
			assert_eq!(
				a.inner_product(&b, &m),
				out
			);
	
			let m = SeaMetric::new(4, 0, 0);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str(" + 1 + e1 + 2e12 + e13 + e234 ").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("  + 1e4 + 2e124 + 3e13 +4e34").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("-3 -4e2 +3e3 -3e4 +3e13 +e23 +2e24 +4e34 +2e124​").unwrap();
			assert_eq!(
				a.inner_product(&b, &m),
				out
			);
		}
	
		#[test]
		fn projective() {
			let m = SeaMetric::new(3, 0, 1);
			let a = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("  + e1 + 2e2 - e234 + e1234 ").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("  - 3e13 + 1e4 - 2e124  +4e34").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("-3e3 +4e14 -5e24").unwrap();
			assert_eq!(
				a.inner_product(&b, &m),
				out
			);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str(" + 1 + e1 + 2e12 + e13 + e234 ").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("  + 3e13 + 1e4 + 2e124  +4e34").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("-3 -3e4 +3e3 +2e24 +4e34 +3e13 +2e124").unwrap();
			assert_eq!(
				a.inner_product(&b, &m),
				out
			);
		}
	
		#[test]
		fn hyperbolic() {
			let m = SeaMetric::new(3, 1, 0);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("1 + 2e1 + 1e12 + 2e3").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("  + 1e2 + 2e13 + 1e23").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("-3e1 -e2 +4e3 +2e13 +e23").unwrap();
			assert_eq!(
				a.inner_product(&b, &m),
				out
			);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str(" + e1 + 2e2 - e234 + e1234  ").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("  - 3e13 + 1e4 - 2e124  +4e34").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("-4e2 -5e3 +4e12 +4e14 +e23 -5e24 -e123").unwrap();
			assert_eq!(
				a.inner_product(&b, &m),
				out
			);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str(" + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e34 + 13e134 + 14e234 + 15e1234").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("  + 3e13 + 1e4 + 2e124  +4e34").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("47 +51e1 +49e2 +53e3 -28e4 +33e12 -10e13 -4e14 -14e23 +47e24 +4e34 -15e123 +2e124").unwrap();
			assert_eq!(
				a.inner_product(&b, &m),
				out
			);
		}
	
		#[test]
		fn other() {
			let m = SeaMetric::new(3, 2, 0);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u32>::from_str(" + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u32>::from_str("  + 3e13 + 1e4 + 2e124  +4e345").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u32>::from_str(" -49 -53e1 -63e2 +3e3 -5e4 -87e12 +3e13 -4e14 +2e24 +14e35 -23e45 +2e124 +13e135 +14e235 +45e245 +4e345 +15e1235").unwrap();
			assert_eq!(
				a.inner_product(&b, &m),
				out
			);
	
			let m = SeaMetric::new(3, 1, 1);
			
			let a = SeaMultivectorLowSparseCoefficents::<f32, u32>::from_str(" + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u32>::from_str("  + 3e13 + 1e4 + 2e124  + 4e345").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u32>::from_str(" -1 -e1 -7e2 +3e3 -5e4 +14e35 -23e45 -27e12 +3e13 -4e14 +2e24 +13e135 +14e235 +45e245 +4e345 +2e124 +15e1235​").unwrap();
			assert_eq!(
				a.inner_product(&b, &m),
				out
			);
		}
	
		#[test]
		fn dimensions_6() {
			let m = SeaMetric::new(3, 3, 0);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u64>::from_str(" + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345 +e6 + 12e126 + 8e46").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u64>::from_str("  + 3e13 + 1e4 + 2e124  +4e345 +83e1235 + 70e1245 + 4e1345 + 8e2345 + 12e123456").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u64>::from_str("115+115e1 -27e2 +1053e3 +1240e4 +189e5 +188e6 -87e12 +3e13 -4e14 +202e15 +168e16 +2e24 -215e25 -156e26 -351e35 -301e45 +2e124 -228e125 -144e126 -185e135 -156e145 +33e235 +83e245 -120e345 -132e356 -84e456 +194e1235 +70e1245 +4e1345 -120e1356 -72e1456 +8e2345 +108e2356 +60e2456 -36e3456 +12e12345 +96e12356 +48e12456 -24e13456 +12e23456 +12e123456​").unwrap();
			// 32.862µs
			assert_eq!(
				a.inner_product(&b, &m),
				out
			);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u64>::from_str("1 + 2e1 + 1e12 + 2e3").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u64>::from_str("  + 1e2 + 2e13 + 1e23").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u64>::from_str("-3e1 -e2 +4e3 +2e13 +e23").unwrap();
			assert_eq!(
				a.inner_product(&b, &m),
				out
			);

			// let now = Instant::now();
			// let num_iters = 10_000;
			// let mut t = a.clone();
			// for _ in 0..num_iters {
			// 	t = t.inner_product(&b, &m);
			// }
			// let elapsed_time = now.elapsed();
			// println!("time time time {:?}", elapsed_time/num_iters);
		}
	}

	mod scalar_product {
		use super::*;
		use num_traits::One;

		#[test]
		fn euclidean() {
			let m = SeaMetric::new(3, 0, 0);
			
			let a = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("1").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("3e12").unwrap();
			assert_eq!(
				a.scalar_product(&b, &m),
				0.0
			);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("e1 + e2").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("e2 + e3").unwrap();
			assert_eq!(
				a.scalar_product(&b, &m),
				1.0
			);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("1 + e1 + e3").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("     e1  + 3e12").unwrap();
			assert_eq!(
				a.scalar_product(&b, &m),
				1.0
			);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("1 + 2e1 + 1e12 + 2e3").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("  + 1e2 + 2e13 + 1e23").unwrap();
			assert_eq!(
				a.scalar_product(&b, &m),
				0.0
			);
	
			let m = SeaMetric::new(4, 0, 0);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str(" + 1 + e1 + 2e12 + e13 + e234 ").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str(" + 3e13 + 1e4 + 2e124 +4e34").unwrap();
			assert_eq!(
				a.scalar_product(&b, &m),
				-3.0
			);
		}
	
		#[test]
		fn projective() {
			let m = SeaMetric::new(3, 0, 1);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("1 + 2e1 + 1e12 + 2e3").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("  + 1e2 + 2e13 + 1e23").unwrap();
			assert_eq!(
				a.scalar_product(&b, &m),
				0.0
			);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str(" + 2e3 + 2e14 + 1e124 ").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("  + 1e24 + 2e134").unwrap();
			assert_eq!(
				a.scalar_product(&b, &m),
				0.0
			);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str(" + 1 + e1 + 2e12 + e13 + e234 ").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("  + 3e13 + 1e4 + 2e124  +4e34").unwrap();
			assert_eq!(
				a.scalar_product(&b, &m),
				-3.0
			);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("  + e1 + 2e2 - e234 + e1234 ").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("  - 3e13 + 1e4 - 2e124  +4e34").unwrap();
			assert_eq!(
				a.scalar_product(&b, &m),
				0.0
			);
		}
	
		#[test]
		fn hyperbolic() {
			let m = SeaMetric::new(3, 1, 0);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("1 + 2e1 + 1e12 + 2e3").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("  + 1e2 + 2e13 + 1e23").unwrap();
			assert_eq!(
				a.scalar_product(&b, &m),
				0.0
			);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str(" + e1 + 2e2 - e234 + e1234  ").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("  - 3e13 + 1e4 - 2e124  +4e34").unwrap();
			assert_eq!(
				a.scalar_product(&b, &m),
				0.0
			);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str(" + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e34 + 13e134 + 14e234 + 15e1234").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("  + 3e13 + 1e4 + 2e124  +4e34").unwrap();
			assert_eq!(
				a.scalar_product(&b, &m),
				47.0
			);
		}
	
		#[test]
		fn other() {
			let m = SeaMetric::new(3, 2, 0);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u32>::from_str(" + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u32>::from_str("  + 3e13 + 1e4 + 2e124  +4e345").unwrap();
	
			assert_eq!(
				a.scalar_product(&b, &m),
				-49.0
			);
	
			let m = SeaMetric::new(3, 1, 1);
			
			let a = SeaMultivectorLowSparseCoefficents::<f32, u32>::from_str(" + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u32>::from_str("  + 3e13 + 1e4 + 2e124  +4e345").unwrap();
	
			assert_eq!(
				a.scalar_product(&b, &m),
				-1.0
			);
		}
	
		#[test]
		fn dimensions_6() {
			let m = SeaMetric::new(3, 3, 0);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u64>::from_str(" + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345 +e6 + 12e126 + 8e46").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u64>::from_str("  + 3e13 + 1e4 + 2e124  +4e345 +83e1235 + 70e1245 + 4e1345 + 8e2345 + 12e123456").unwrap();
	
			assert_eq!(
				a.scalar_product(&b, &m),
				115.0
			);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u64>::from_str("1 + 2e1 + 1e12 + 2e3").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u64>::from_str("  + 1e2 + 2e13 + 1e23").unwrap();
			assert_eq!(
				a.scalar_product(&b, &m),
				0.0
			);
		}
	
	}

	mod outer_product {
		use super::*;
		#[test]
		fn works() {
			let a = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("1").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("3e12").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("3e12").unwrap();
			assert_eq!(
				a.outer_product(&b),
				out
			);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("e1 + e2").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("e2 + e3").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("e12 + e13 + e23").unwrap();
			assert_eq!(
				a.outer_product(&b),
				out
			);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("1 + e1 + e3").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("     e1  + 3e12").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("e1 + 3e12 - e13 + 3e123").unwrap();
			assert_eq!(
				a.outer_product(&b),
				out
			);

			let a = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("e12").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("e2").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("").unwrap();
			assert_eq!(
				a.outer_product(&b),
				out
			);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("1 + 2e1 + 1e12 + 2e3").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("  + 1e2 + 2e13 + 1e23").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("e2 +2e12 +2e13 -e23 +2e123").unwrap();
			assert_eq!(
				a.outer_product(&b),
				out
			);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str(" + 1 + e1 + 2e12 + e13 + e234 ").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("  + 1e4 + 2e124 + 3e13 +4e34").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("e4 +3e13 +e14 +4e34 +4e124 +5e134 +8e1234").unwrap();
			assert_eq!(
				a.outer_product(&b),
				out
			);

			let a = SeaMultivectorLowSparseCoefficents::<f32, u32>::from_str(" + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u32>::from_str("  + 3e13 + 1e4 + 2e124  +4e345").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u32>::from_str("e4 +3e13 +e14 +2e24 +4e34 -6e123 +5e124 +29e134 +6e234 +4e345 -15e1234 +4e1345 +8e2345 +12e12345").unwrap();
			assert_eq!(
				a.outer_product(&b),
				out
			);

			let a = SeaMultivectorLowSparseCoefficents::<f32, u32>::from_str("15 +14e1 -13e2 +7e45 -11e35 +12e12 +6e145 -5e245 +3e345 -10e135 +9e235 +4e1245 -2e1345 +e2345 -8e1235 +e12345").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u32>::from_str("-2e35 +4e12 +3e345 -e1235").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u32>::from_str("-30e35 +60e12 +45e345 -28e135 +26e235 +28e1245 +42e1345 -39e2345 -83e1235 +48e12345").unwrap();
			assert_eq!(
				a.outer_product(&b),
				out
			);

			let a = SeaMultivectorLowSparseCoefficents::<f32, u64>::from_str(" + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345 +e6 + 12e126 + 8e46").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u64>::from_str("  + 3e13 + 1e4 + 2e124  +4e345 +83e1235 + 70e1245 + 4e1345 + 8e2345 + 12e123456").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u64>::from_str("e4 +3e13 +e14 +2e24 +4e34 -e46 -6e123 +5e124 +29e134 +3e136 +6e234 +4e345 -15e1234 +83e1235 +70e1245 -14e1246 +8e1345 +24e1346 +16e2345 -4e3456 -372e12345 +83e12356 +70e12456 +4e13456 +8e23456 -700e123456").unwrap();
			//6.69µs
			assert_eq!(
				a.outer_product(&b),
				out
			);

			// let now = Instant::now();
			// let num_iters = 10_000;
			// let mut t = a.clone();
			// for _ in 0..num_iters {
			// 	t = a.outer_product(&b);
			// }
			// let elapsed_time = now.elapsed();
			// println!("time time time {:?}", elapsed_time/num_iters);
		}
	}

	mod regressive_product{
		use super::*;

		#[test]
		fn euclidean() {
			let m = SeaMetric::new(3, 0, 0);
			
			let a = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("1").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("3e12").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("").unwrap();
			assert_eq!(
				a.regressive_product(&b, &m),
				out
			);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("e1 + e2").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("e2 + e3").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("").unwrap();
			assert_eq!(
				a.regressive_product(&b, &m),
				out
			);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("1 + e1 + e3").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("     e1  + 3e12").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("3").unwrap();
			assert_eq!(
				a.regressive_product(&b, &m),
				out
			);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("1 + 2e1 + 1e12 + 2e3").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("  + 1e2 + 2e13 + 1e23").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("2 -2e1 -e2").unwrap();
			assert_eq!(
				a.regressive_product(&b, &m),
				out
			);
	
			let m = SeaMetric::new(4, 0, 0);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str(" + 1 + e1 + 2e12 + e13 + e234 ").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("  + 1e4 + 2e124 + 3e13 +4e34").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("8 -2e1 +3e3 +2e24").unwrap();
			assert_eq!(
				a.regressive_product(&b, &m),
				out
			);
		}
	
		#[test]
		fn projective() {
			let m = SeaMetric::new(3, 0, 1);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("1 + 2e1 + 1e12 + 2e3").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("  + 1e2 + 2e13 + 1e23").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("").unwrap();
			assert_eq!(
				a.regressive_product(&b, &m),
				out
			);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str(" + 2e3 + 2e14 + 1e124 ").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("  + 1e24 + 2e134").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("2e14").unwrap();
			assert_eq!(
				a.regressive_product(&b, &m),
				out
			);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str(" + 1 + e1 + 2e12 + e13 + e234 ").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("  + 3e13 + 1e4 + 2e124  +4e34").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("+8 -2e1 +3e3 -2e24").unwrap();
			assert_eq!(
				a.regressive_product(&b, &m),
				out
			);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("  + e1 + 2e2 - e234 + e1234 ").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("  - 3e13 + 1e4 - 2e124  +4e34").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("+e4 +3e3 -2e24 +4e34 -3e13 -2e124").unwrap();
			assert_eq!(
				a.regressive_product(&b, &m),
				out
			);

			let m = SeaMetric::new(3, 0, 2);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u32>::from_str(" + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u32>::from_str("  + 3e12 + 1e4 + 2e124  +4e345").unwrap();
			
			let out = SeaMultivectorLowSparseCoefficents::<f32, u32>::from_str("48 +83e4 -39e1 -42e2 +28e3 +26e14 +28e24 +45e12 +60e345 +30e124").unwrap();
			assert_eq!(
				a.regressive_product(&b, &m),
				out
			);
		}
	
		#[test]
		fn hyperbolic() {
			let m = SeaMetric::new(3, 1, 0);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("1 + 2e1 + 1e12 + 2e3").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("  + 1e2 + 2e13 + 1e23").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("").unwrap();
			assert_eq!(
				a.regressive_product(&b, &m),
				out
			);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str(" + e1 + 2e2 - e234 + e1234  ").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("  - 3e13 + 1e4 - 2e124  +4e34").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("3e3 +e4 -3e13 +2e24 +4e34 -2e124").unwrap();
			assert_eq!(
				a.regressive_product(&b, &m),
				out
			);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str(" + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e34 + 13e134 + 14e234 + 15e1234").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("  + 3e13 + 1e4 + 2e124  +4e34").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str(" -33 -43e1 -12e2 +70e3 +83e4 -14e12 +45e13 +26e14 +28e24 +60e34 +30e124").unwrap();
			assert_eq!(
				a.regressive_product(&b, &m),
				out
			);
		}
	
		#[test]
		fn other() {
			let m = SeaMetric::new(3, 2, 0);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u32>::from_str(" + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u32>::from_str("  + 3e13 + 1e4 + 2e124  +4e345").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u32>::from_str("12 +70e3 +83e4 +45e13 +26e14 +28e24 +30e124 +60e345").unwrap();
			assert_eq!(
				a.regressive_product(&b, &m),
				out
			);
	
			let m = SeaMetric::new(3, 1, 1);
			
			let a = SeaMultivectorLowSparseCoefficents::<f32, u32>::from_str(" + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u32>::from_str("  + 3e13 + 1e4 + 2e124  +4e345").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u32>::from_str("12 -14e3 +83e4 +45e13 +26e14 +28e24 +60e345 +30e124").unwrap();
			assert_eq!(
				a.regressive_product(&b, &m),
				out
			);
		}
	
		#[test]
		fn dimensions_6() {
			let m = SeaMetric::new(3, 3, 0);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u64>::from_str(" + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345 +e6 + 12e126 + 8e46").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u64>::from_str("  + 3e13 + 1e4 + 2e124  +4e345 +83e1235 + 70e1245 + 4e1345 + 8e2345 + 12e123456").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u64>::from_str(" -604 -36e1 -72e2 +48e3 +96e4 +12e6 +36e12 +60e13 +108e14 +72e23 +120e24 +96e46 +84e123 +132e124 +144e126 +144e345 +156e1345 +168e2345 +180e12345").unwrap();
			assert_eq!(
				a.regressive_product(&b, &m),
				out
			);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u64>::from_str("1 + 2e1 + 1e12 + 2e3").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u64>::from_str("  + 1e2 + 2e13 + 1e23").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u64>::from_str("").unwrap();
			assert_eq!(
				a.regressive_product(&b, &m),
				out
			);

			// let now = Instant::now();
			// let num_iters = 10_000;
			// let mut t = a.clone();
			// for _ in 0..num_iters {
			// 	t = t.geometric_product(&b, &m);
			// }
			// let elapsed_time = now.elapsed();
			// println!("time time time {:?}", elapsed_time/num_iters);
		}
	
	}
}