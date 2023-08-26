use super::*;

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