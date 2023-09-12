use super::*;

// Scalar Norm
impl<T:SeaCoefficient, I:SeaBasisBladeBits> SeaMultivectorLowSparseCoefficents<T, I> {
	pub fn scalar_norm(&self, metric: &SeaMetric) -> T {
		let (_, metric_non_zeros, metric_negs) = Self::metric_product_constants(metric);

		
		let mut total_neg_mask = I::total_negations(&SeaBasisBladeIndexLow::SCALAR, &metric_negs);
		total_neg_mask.symmetric_difference_mut(&I::reversion_negations());
		
		let zeros = I::outer_product_zeros(&metric_non_zeros);
		let dot_mask = self.basis_blade_bits.intersection(&zeros);

		let mut out = T::zero();
		for i in dot_mask.iter() {
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

#[cfg(test)]
mod tests{
	use std::time::Instant;

	use crate::tests::catch_unwind_silent;

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