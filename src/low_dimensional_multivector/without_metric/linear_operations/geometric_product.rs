use super::*;

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