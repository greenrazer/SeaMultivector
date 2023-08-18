use std::iter::repeat;
use std::ops::{Add, AddAssign};
use std::str::FromStr;

use lazy_static::lazy_static;
use num_traits::Zero;
use regex::Regex;

use crate::coefficient::SeaCoefficient;
use crate::errors::ParseMultivectorError;
use crate::metric::SeaMetric;
use crate::multivector_traits::{Multivector, MultivectorMut};
use crate::multivector_product_traits::{LeftContration, RightContraction, InnerProduct, OuterProduct, GeometricProduct};

use super::basis_blade_bits::SeaBasisBladeBitsInner;
use super::basis_blade_bits::SeaBasisBladeBits;
use super::basis_blade::SeaBasisBladeLow;
use super::basis_blade_index::SeaBasisBladeIndexLow;


#[derive(Debug, Clone, PartialEq)]
pub struct SeaMultivectorLowSparseCoefficents<T:SeaCoefficient, I:SeaBasisBladeBitsInner>{
	pub basis_blade_bits: SeaBasisBladeBits<I>,
	pub coefficients: Vec<T>
}
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

	fn special_dot_product(self: &Self, rhs: &Self, dot_mask: &SeaBasisBladeBits<I>, neg_mask: &SeaBasisBladeBits<I>, output_index: &SeaBasisBladeIndexLow) -> T {
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

impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> Multivector for SeaMultivectorLowSparseCoefficents<T, I>{
	fn new() -> Self {
		Self { 
			basis_blade_bits: SeaBasisBladeBits::ZERO, 
			coefficients: vec![]
		}
	}

	fn reversion_conjugated(&self) -> Self {
		self.negated_basis_blades(&SeaBasisBladeBits::<I>::REVERSION_NEGATIONS)
	}

	fn involution_conjugated(&self) -> Self {
		self.negated_basis_blades(&SeaBasisBladeBits::<I>::ODD_GRADES)
	}

	fn clifford_conjugated(&self) -> Self {
		self.negated_basis_blades(&SeaBasisBladeBits::<I>::CONJUGATION_NEGATIONS)
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> MultivectorMut for SeaMultivectorLowSparseCoefficents<T, I>{
	fn reversion_conjugation(&mut self) {
		self.negate_basis_blades(&SeaBasisBladeBits::<I>::REVERSION_NEGATIONS)
	}

	fn involution_conjugation(&mut self) {
		self.negate_basis_blades(&SeaBasisBladeBits::<I>::ODD_GRADES)
	}

	fn clifford_conjugation(&mut self) {
		self.negate_basis_blades(&SeaBasisBladeBits::<I>::CONJUGATION_NEGATIONS)
	}
}

impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> LeftContration<(&Self, &SeaMetric)> for SeaMultivectorLowSparseCoefficents<T, I> {
	type Output = Self;

	fn left_contraction(&self, (rhs, metric): (&Self, &SeaMetric)) -> Self::Output {
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
		let mut out_idx = 0;
		for o in &output_basis_blade_bits {
			Self::update_zero_cache(&mut zero_cache, &o, &metric_non_zeros);

			let total_neg_mask = Self::total_neg_mask(&o, &metric_negs);

			let mut curr_rhs_bbb = &rhs.basis_blade_bits >> &o;
			curr_rhs_bbb.remove_basis_blades_with_any(&o);

			let mut dot_mask = &self.basis_blade_bits & &curr_rhs_bbb;
			dot_mask &= &zero_cache.1;
			
			output_coefficients[out_idx] = Self::special_dot_product(
				&self,
				&rhs,
				&dot_mask,
				&total_neg_mask,
				&o
			);

			out_idx += 1;
		}

		Self {
			basis_blade_bits: output_basis_blade_bits,
			coefficients: output_coefficients
		}
	}
}

impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> RightContraction<(&Self, &SeaMetric)> for SeaMultivectorLowSparseCoefficents<T, I> {
	type Output = Self;

	fn right_contraction(&self, (rhs, metric): (&Self, &SeaMetric)) -> Self::Output {
		let (num_metric_non_zeros, metric_non_zeros, metric_negs) = Self::metric_product_constants(metric);

		let output_basis_blade_bits = {
			self.basis_blade_bits.right_contraction(
				(&rhs.basis_blade_bits, num_metric_non_zeros)
			)
		};

		let output_num_coeffs = output_basis_blade_bits.count_ones() as usize;
		let mut output_coefficients = vec![T::zero(); output_num_coeffs];

		let mut zero_cache: (SeaBasisBladeIndexLow, SeaBasisBladeBits<I>) = (SeaBasisBladeIndexLow::ZERO, SeaBasisBladeBits::<I>::ONE);
		let mut out_idx = 0;
		for o in &output_basis_blade_bits {
			Self::update_zero_cache(&mut zero_cache, &o, &metric_non_zeros);

			let total_neg_mask = Self::total_neg_mask(&o, &metric_negs);

			let mut curr_rhs_bbb = &rhs.basis_blade_bits << &o;
			curr_rhs_bbb.remove_basis_blades_without_all(&o);

			let mut dot_mask = &self.basis_blade_bits & &curr_rhs_bbb;
			dot_mask &= &zero_cache.1;

			output_coefficients[out_idx] = Self::special_dot_product(
				&self,
				&rhs,
				&dot_mask,
				&total_neg_mask,
				&o
			);

			out_idx += 1;
		}

		Self {
			basis_blade_bits: output_basis_blade_bits,
			coefficients: output_coefficients
		}
	}
}

impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> InnerProduct<(&Self, &SeaMetric)> for SeaMultivectorLowSparseCoefficents<T, I> {
	type Output = Self;
	fn inner_product(&self, (rhs, metric): (&Self, &SeaMetric)) -> Self::Output {
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
		let mut out_idx = 0;
		for o in &output_basis_blade_bits {
			Self::update_zero_cache(&mut zero_cache, &o, &metric_non_zeros);

			let total_neg_mask = Self::total_neg_mask(&o, &metric_negs);

			let mut left_contraction_rhs_bbb = &rhs.basis_blade_bits >> &o;
			left_contraction_rhs_bbb.remove_basis_blades_with_any(&o);

			let mut dot_mask_left = &self.basis_blade_bits & &left_contraction_rhs_bbb;
			dot_mask_left &= &zero_cache.1;

			output_coefficients[out_idx] = Self::special_dot_product(
				&self,
				&rhs,
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
					&self,
					&rhs,
					&dot_mask_right,
					&total_neg_mask,
					&o
				);
			}

			
			out_idx += 1;
		}

		Self {
			basis_blade_bits: output_basis_blade_bits,
			coefficients: output_coefficients
		}
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

		let mut out_idx = 0;
		for o in &output_basis_blade_bits {
			let cannonical_neg_mask = SeaBasisBladeBits::<I>::cannonical_reordering_neg_mask(&o);

			let mut curr_rhs_bbb = &reverse_rhs_bbbs >> &(&initial_shift - &o);
			curr_rhs_bbb &= &SeaBasisBladeBits::<I>::zero_mask_for_output_index(&o);

			let dot_mask = &self.basis_blade_bits & &curr_rhs_bbb;

			output_coefficients[out_idx] = Self::special_dot_product(
				&self,
				&rhs,
				&dot_mask,
				&cannonical_neg_mask,
				&o
			);
			
			out_idx += 1;
		}

		Self {
			basis_blade_bits: output_basis_blade_bits,
			coefficients: output_coefficients
		}
	}
}

impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> GeometricProduct<(&Self, &SeaMetric)> for SeaMultivectorLowSparseCoefficents<T, I>{ 
	type Output = Self;

	fn geometric_product(&self, (rhs, metric): (&Self, &SeaMetric)) -> Self::Output {
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
		let mut out_idx = 0;
		for o in &output_basis_blade_bits {
			Self::update_zero_cache(&mut zero_cache, &o, &metric_non_zeros);

			let total_neg_mask = Self::total_neg_mask(&o, &metric_negs);

			let curr_rhs_bbb = rhs.basis_blade_bits.geometric_product_swap(&o);

			let mut dot_mask = &self.basis_blade_bits & &curr_rhs_bbb;
			dot_mask &= &zero_cache.1;

			output_coefficients[out_idx] = Self::special_dot_product(
				&self,
				&rhs,
				&dot_mask,
				&total_neg_mask,
				&o
			);

			out_idx += 1;
		}

		Self {
			basis_blade_bits: output_basis_blade_bits,
			coefficients: output_coefficients
		}
	}

}

#[cfg(test)]
mod tests{
	use std::time::Instant;

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
				a.geometric_product((&b, &m)),
				out
			);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("e1 + e2").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("e2 + e3").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("1 + e12 + e13 + e23").unwrap();
			//3.382µs
			assert_eq!(
				a.geometric_product((&b, &m)),
				out
			);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("1 + e1 + e3").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("     e1  + 3e12").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("1 + e1 + 3e2 +3e12 -e13 + 3e123").unwrap();
			assert_eq!(
				a.geometric_product((&b, &m)),
				out
			);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("1 + 2e1 + 1e12 + 2e3").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("  + 1e2 + 2e13 + 1e23").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("-3e1 - e2 + 2e12 + 4e3 + 3e13 -3e23 + 2e123").unwrap();
			//7.542µs
			assert_eq!(
				a.geometric_product((&b, &m)),
				out
			);
	
			let m = SeaMetric::new(4, 0, 0);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str(" + 1 + e1 + 2e12 + e13 + e234 ").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("  + 1e4 + 2e124 + 3e13 +4e34").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("-3 -4e2 + 3e3 -3e4 +e13 +5e14 -5e23 + 2e24 + 4e34 + 7e124 + 5e134 + 2e234 + 8e1234").unwrap();
			//21.11µs
			//10.495µs
			assert_eq!(
				a.geometric_product((&b, &m)),
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
				a.geometric_product((&b, &m)),
				out
			);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str(" + 2e3 + 2e14 + 1e124 ").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("  + 1e24 + 2e134").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("-4e14 -2e234").unwrap();
			assert_eq!(
				a.geometric_product((&b, &m)),
				out
			);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str(" + 1 + e1 + 2e12 + e13 + e234 ").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("  + 3e13 + 1e4 + 2e124  +4e34").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("-3 +3e3 + 3e13 -6e23 -3e4  + 5e14 +2e24 + 7e124 + 4e34  +5e134 + 2e234 + 8e1234").unwrap();
			//12.67µs
			assert_eq!(
				a.geometric_product((&b, &m)),
				out
			);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("  + e1 + 2e2 - e234 + e1234 ").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("  - 3e13 + 1e4 - 2e124  +4e34").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("-3e3 +6e123 + 5e14 - 3e24 +3e124 + 4e134 + 8e234 ").unwrap();
			//8.403µs
			assert_eq!(
				a.geometric_product((&b, &m)),
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
				a.geometric_product((&b, &m)),
				out
			);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str(" + e1 + 2e2 - e234 + e1234  ").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("  - 3e13 + 1e4 - 2e124  +4e34").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("-4e2 -5e3 + 4e12 + 2e13 + 5e14 + e23 -3e24 + 5e123 + 3e124 + 4e134 + 8e234").unwrap();
			assert_eq!(
				a.geometric_product((&b, &m)),
				out
			);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str(" + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e34 + 13e134 + 14e234 + 15e1234").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("  + 3e13 + 1e4 + 2e124  +4e34").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("47 +51e1 + 49e2 + 53e3 - 28e4 + 51e12 + 54e13 - 19e14 - 9e23 + 73e24 + 21e34 - e123 + 75e124 + 21e134 -9e234 -3e1234").unwrap();
			//21.547µs
			assert_eq!(
				a.geometric_product((&b, &m)),
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
				a.geometric_product((&b, &m)),
				out
			);
	
			let m = SeaMetric::new(3, 1, 1);
			
			let a = SeaMultivectorLowSparseCoefficents::<f32, u32>::from_str(" + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u32>::from_str("  + 3e13 + 1e4 + 2e124  +4e345").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u32>::from_str("-1 -e1 -7e2 +3e3 - 5e4 +14e35 -23e45 -9e12 +3e13 - 3e14 -9e23 +4e24 +17e34 +21e135 - 16e145 +80e235 +69e245 + 4e345 -6e123 + 5e124 + 17e134 - 17e234 +83e1235 +70e1245 +4e1345 + 8e2345 -15e1234 +12e12345").unwrap();
			// 45.752µs
			// new 15.492µs
			assert_eq!(
				a.geometric_product((&b, &m)),
				out
			);
		}
	
		#[test]
		fn dimensions_6() {
			let m = SeaMetric::new(3, 3, 0);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u64>::from_str(" + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345 +e6 + 12e126 + 8e46").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u64>::from_str("  + 3e13 + 1e4 + 2e124  +4e345 +83e1235 + 70e1245 + 4e1345 + 8e2345 + 12e123456").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u64>::from_str("115+115e1-27e2+1053e3+1240e4+189e5+188e6-21e12+983e13+1159e14+202e15+168e16-919e23-1075e24-215e25-156e26+17e34-351e35-301e45+23e46-846e123-991e124-228e125-128e126+17e134-265e135+3e136-228e145-17e234+143e235-36e236+135e245+303e345-1160e356-924e456-15e1234+230e1235+124e1245-14e1246-560e1256+442e1345+24e1346-152e1356-72e1456-393e2345+44e2356+60e2456-40e3456-360e12345+179e12356+118e12456+76e13456-28e23456-700e123456​​").unwrap();
			// 171.677µs
			// new: 58.328µs
			// with unsafe: 56.79µs
			// 49.433µs
			assert_eq!(
				a.geometric_product((&b, &m)),
				out
			);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u64>::from_str("1 + 2e1 + 1e12 + 2e3").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u64>::from_str("  + 1e2 + 2e13 + 1e23").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u64>::from_str("-3e1 - e2 + 2e12 + 4e3 + 3e13 -3e23 + 2e123").unwrap();
			//13.991µs
			//3.334µs
			assert_eq!(
				a.geometric_product((&b, &m)),
				out
			);

			// let now = Instant::now();
			// let num_iters = 10_000;
			// let mut t = a.clone();
			// for _ in 0..num_iters {
			// 	t = t.geometric_product((&b, &m));
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
				a.left_contraction((&b, &m)),
				out
			);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("e1 + e2").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("e2 + e3").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("1").unwrap();
			assert_eq!(
				a.left_contraction((&b, &m)),
				out
			);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("1 + e1 + e3").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("     e1  + 3e12").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("1 + e1 + 3e2 +3e12").unwrap();
			assert_eq!(
				a.left_contraction((&b, &m)),
				out
			);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("1 + 2e1 + 1e12 + 2e3").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("  + 1e2 + 2e13 + 1e23").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("-4e1 - e2 + 4e3 + 2e13 +1e23").unwrap();
			assert_eq!(
				a.left_contraction((&b, &m)),
				out
			);
	
			let m = SeaMetric::new(4, 0, 0);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str(" + 1 + e1 + 2e12 + e13 + e234 ").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("  + 1e4 + 2e124 + 3e13 +4e34").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("-3 + 3e3 -3e4 + 3e13 + 2e24 + 4e34 + 2e124").unwrap();
			assert_eq!(
				a.left_contraction((&b, &m)),
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
				a.left_contraction((&b, &m)),
				out
			);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str(" + 1 + e1 + 2e12 + e13 + e234 ").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("  + 3e13 + 1e4 + 2e124  +4e34").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("-3 -3e4 + 3e3 +2e24 + 4e34 + 3e13 + 2e124").unwrap();
			assert_eq!(
				a.left_contraction((&b, &m)),
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
				a.left_contraction((&b, &m)),
				out
			);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str(" + e1 + 2e2 - e234 + e1234  ").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("  - 3e13 + 1e4 - 2e124  +4e34").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("-3e3 + 4e14 -2e24").unwrap();
			assert_eq!(
				a.left_contraction((&b, &m)),
				out
			);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str(" + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e34 + 13e134 + 14e234 + 15e1234").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("  + 3e13 + 1e4 + 2e124  +4e34").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("47 +8e1 -18e2 +35e3 + 11e4 -16e12 + 3e13 -4e14 +2e24 +4e34 + 2e124").unwrap();
			//4.371µs
			assert_eq!(
				a.left_contraction((&b, &m)),
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
				a.left_contraction((&b, &m)),
				out
			);
	
			let m = SeaMetric::new(3, 1, 1);
			
			let a = SeaMultivectorLowSparseCoefficents::<f32, u32>::from_str(" + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u32>::from_str("  + 3e13 + 1e4 + 2e124  +4e345").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u32>::from_str("-1 +8e1 -18e2 + 3e3 -5e4 +32e35 +16e45 -16e12 + 3e13 -4e14 +2e24 +4e345 +2e124").unwrap();
			//3.867µs
			assert_eq!(
				a.left_contraction((&b, &m)),
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
				a.left_contraction((&b, &m)),
				out
			);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u64>::from_str("1 + 2e1 + 1e12 + 2e3").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u64>::from_str("  + 1e2 + 2e13 + 1e23").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u64>::from_str("-4e1 -e2 +4e3 +2e13 +e23").unwrap();
			//2.081µs
			assert_eq!(
				a.left_contraction((&b, &m)),
				out
			);

			// let now = Instant::now();
			// let num_iters = 10_000;
			// let mut t = a.clone();
			// for _ in 0..num_iters {
			// 	t = t.left_contraction((&b, &m));
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
				a.right_contraction((&b, &m)),
				out
			);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("e1 + e2").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("e2 + e3").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("1").unwrap();
			assert_eq!(
				a.right_contraction((&b, &m)),
				out
			);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("1 + e1 + e3").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("     e1  + 3e12").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("1").unwrap();
			assert_eq!(
				a.right_contraction((&b, &m)),
				out
			);

			let a = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("e12").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("e2").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("e1").unwrap();
			assert_eq!(
				a.right_contraction((&b, &m)),
				out
			);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("1 + 2e1 + 1e12 + 2e3").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("  + 1e2 + 2e13 + 1e23").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("e1").unwrap();
			assert_eq!(
				a.right_contraction((&b, &m)),
				out
			);
	
			let m = SeaMetric::new(4, 0, 0);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str(" + 1 + e1 + 2e12 + e13 + e234 ").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("  + 1e4 + 2e124 + 3e13 +4e34").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("-3 -4e2 + e23").unwrap();
			assert_eq!(
				a.right_contraction((&b, &m)),
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
				a.right_contraction((&b, &m)),
				out
			);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str(" + 1 + e1 + 2e12 + e13 + e234 ").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("  + 3e13 + 1e4 + 2e124  +4e34").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("-3").unwrap();
			assert_eq!(
				a.right_contraction((&b, &m)),
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
				a.right_contraction((&b, &m)),
				out
			);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str(" + e1 + 2e2 - e234 + e1234  ").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("  - 3e13 + 1e4 - 2e124  +4e34").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("-4e2 - 2e3 +4e12 + e23 - 3e24 - e123").unwrap();
			assert_eq!(
				a.right_contraction((&b, &m)),
				out
			);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str(" + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e34 + 13e134 + 14e234 + 15e1234").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("  + 3e13 + 1e4 + 2e124  +4e34").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("47+43e1+67e2+18e3-39e4+49e12-13e13-14e23+45e24-15e123").unwrap();
			assert_eq!(
				a.right_contraction((&b, &m)),
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
				a.right_contraction((&b, &m)),
				out
			);
	
			let m = SeaMetric::new(3, 1, 1);
			
			let a = SeaMultivectorLowSparseCoefficents::<f32, u32>::from_str(" + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u32>::from_str("  + 3e13 + 1e4 + 2e124  +4e345").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u32>::from_str("-1 -9e1 +11e2 -18e35 -39e45 -11e12 +13e135 +14e235 +45e245 +15e1235").unwrap();
			assert_eq!(
				a.right_contraction((&b, &m)),
				out
			);
		}
	
		#[test]
		fn dimensions_6() {
			let m = SeaMetric::new(3, 3, 0);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u64>::from_str(" + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345 +e6 + 12e126 + 8e46").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u64>::from_str("  + 3e13 + 1e4 + 2e124  +4e345 +83e1235 + 70e1245 + 4e1345 + 8e2345 + 12e123456").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u64>::from_str("115+59e1-105e2+1050e3+1245e4+8e6-71e12-18e35-39e45+13e135+14e235+45e245+15e1235​").unwrap();
			//40ns
			assert_eq!(
				a.right_contraction((&b, &m)),
				out
			);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u64>::from_str("1 + 2e1 + 1e12 + 2e3").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u64>::from_str("  + 1e2 + 2e13 + 1e23").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u64>::from_str("e1").unwrap();
			assert_eq!(
				a.right_contraction((&b, &m)),
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
				a.inner_product((&b, &m)),
				out
			);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("e1 + e2").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("e2 + e3").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("1").unwrap();
			assert_eq!(
				a.inner_product((&b, &m)),
				out
			);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("1 + e1 + e3").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("     e1  + 3e12").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("1 + e1 + 3e2 + 3e12").unwrap();
			assert_eq!(
				a.inner_product((&b, &m)),
				out
			);

			let a = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("e12").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("e2").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("e1").unwrap();
			assert_eq!(
				a.inner_product((&b, &m)),
				out
			);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("1 + 2e1 + 1e12 + 2e3").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("  + 1e2 + 2e13 + 1e23").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("-3e1 -e2 +4e3 +2e13 + e23").unwrap();
			assert_eq!(
				a.inner_product((&b, &m)),
				out
			);
	
			let m = SeaMetric::new(4, 0, 0);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str(" + 1 + e1 + 2e12 + e13 + e234 ").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("  + 1e4 + 2e124 + 3e13 +4e34").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("-3 -4e2 +3e3 -3e4 +3e13 +e23 +2e24 +4e34 +2e124​").unwrap();
			assert_eq!(
				a.inner_product((&b, &m)),
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
				a.inner_product((&b, &m)),
				out
			);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str(" + 1 + e1 + 2e12 + e13 + e234 ").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("  + 3e13 + 1e4 + 2e124  +4e34").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("-3 -3e4 +3e3 +2e24 +4e34 +3e13 +2e124").unwrap();
			assert_eq!(
				a.inner_product((&b, &m)),
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
				a.inner_product((&b, &m)),
				out
			);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str(" + e1 + 2e2 - e234 + e1234  ").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("  - 3e13 + 1e4 - 2e124  +4e34").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("-4e2 -5e3 +4e12 +4e14 +e23 -5e24 -e123").unwrap();
			assert_eq!(
				a.inner_product((&b, &m)),
				out
			);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str(" + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e34 + 13e134 + 14e234 + 15e1234").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("  + 3e13 + 1e4 + 2e124  +4e34").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("47 +51e1 +49e2 +53e3 -28e4 +33e12 -10e13 -4e14 -14e23 +47e24 +4e34 -15e123 +2e124").unwrap();
			assert_eq!(
				a.inner_product((&b, &m)),
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
				a.inner_product((&b, &m)),
				out
			);
	
			let m = SeaMetric::new(3, 1, 1);
			
			let a = SeaMultivectorLowSparseCoefficents::<f32, u32>::from_str(" + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u32>::from_str("  + 3e13 + 1e4 + 2e124  + 4e345").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u32>::from_str(" -1 -e1 -7e2 +3e3 -5e4 +14e35 -23e45 -27e12 +3e13 -4e14 +2e24 +13e135 +14e235 +45e245 +4e345 +2e124 +15e1235​").unwrap();
			assert_eq!(
				a.inner_product((&b, &m)),
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
				a.inner_product((&b, &m)),
				out
			);
	
			let a = SeaMultivectorLowSparseCoefficents::<f32, u64>::from_str("1 + 2e1 + 1e12 + 2e3").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u64>::from_str("  + 1e2 + 2e13 + 1e23").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u64>::from_str("-3e1 -e2 +4e3 +2e13 +e23").unwrap();
			assert_eq!(
				a.inner_product((&b, &m)),
				out
			);

			// let now = Instant::now();
			// let num_iters = 10_000;
			// let mut t = a.clone();
			// for _ in 0..num_iters {
			// 	t = t.inner_product((&b, &m));
			// }
			// let elapsed_time = now.elapsed();
			// println!("time time time {:?}", elapsed_time/num_iters);
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

			let a = SeaMultivectorLowSparseCoefficents::<f32, u64>::from_str(" + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345 +e6 + 12e126 + 8e46").unwrap();
			let b = SeaMultivectorLowSparseCoefficents::<f32, u64>::from_str("  + 3e13 + 1e4 + 2e124  +4e345 +83e1235 + 70e1245 + 4e1345 + 8e2345 + 12e123456").unwrap();
	
			let out = SeaMultivectorLowSparseCoefficents::<f32, u64>::from_str("e4 +3e13 +e14 +2e24 +4e34 -e46 -6e123 +5e124 +29e134 +3e136 +6e234 +4e345 -15e1234 +83e1235 +70e1245 -14e1246 +8e1345 +24e1346 +16e2345 -4e3456 -372e12345 +83e12356 +70e12456 +4e13456 +8e23456 -700e123456").unwrap();
			//39ns
			//xor 37ns
			assert_eq!(
				a.outer_product(&b),
				out
			);

			// let now = Instant::now();
			// let num_iters = 10_000;
			// let mut t = a.clone();
			// for _ in 0..num_iters {
			// 	t = t.outer_product(&b);
			// }
			// let elapsed_time = now.elapsed();
			// println!("time time time {:?}", elapsed_time/num_iters);
		}
	}

}