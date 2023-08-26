use super::*;

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
