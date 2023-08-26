use super::*;

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
}

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
}