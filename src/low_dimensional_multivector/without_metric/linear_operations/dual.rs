use super::*;

// Dual
impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> SeaMultivectorLowSparseCoefficents<T, I> {
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

//Dual Mut
impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> SeaMultivectorLowSparseCoefficents<T, I> {
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
