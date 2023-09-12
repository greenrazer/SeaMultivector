use super::*;

impl<T:SeaCoefficient, I:SeaBasisBladeSet> SeaMultivectorLowSparseCoefficents<T, I>{
	pub fn reversion_conjugated(&self) -> Self {
		self.negated_basis_blades(&I::reversion_negations())
	}

	pub fn involution_conjugated(&self) -> Self {
		self.negated_basis_blades(&I::involution_negations())
	}

	pub fn clifford_conjugated(&self) -> Self {
		self.negated_basis_blades(&I::clifford_negations())
	}
}

impl<T:SeaCoefficient, I:SeaBasisBladeSet> SeaMultivectorLowSparseCoefficents<T, I>{
	pub fn reversion_conjugation(&mut self) {
		self.negate_basis_blades(&I::reversion_negations())
	}

	pub fn involution_conjugation(&mut self) {
		self.negate_basis_blades(&I::involution_negations())
	}

	pub fn clifford_conjugation(&mut self) {
		self.negate_basis_blades(&I::clifford_negations())
	}
}