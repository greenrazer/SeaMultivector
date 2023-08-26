use super::*;

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
