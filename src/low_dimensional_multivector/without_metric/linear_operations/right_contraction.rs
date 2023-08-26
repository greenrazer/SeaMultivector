use super::*;

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
