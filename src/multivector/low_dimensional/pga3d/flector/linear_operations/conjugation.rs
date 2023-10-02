use super::*;

impl<T:SeaCoefficient> MultivectorConjugation for SeaFlectorPGA3D<T> {
	fn reversion_conjugated(&self) -> Self {
		Self((**self).reversion_conjugated())
	}

	fn involution_conjugated(&self) -> Self {
		Self((**self).involution_conjugated())
	}

	fn clifford_conjugated(&self) -> Self {
		Self((**self).clifford_conjugated())
	}

	fn negated(&self) -> Self {
		Self((**self).negated())	
	}
}
impl<T:SeaCoefficient> MultivectorConjugationMut for SeaFlectorPGA3D<T> {	
	fn reversion_conjugation(&mut self) {
		(**self).reversion_conjugation()
	}

	fn involution_conjugation(&mut self) {
		(**self).involution_conjugation()
	}

	fn clifford_conjugation(&mut self) {
		(**self).clifford_conjugation()
	}

	fn negate(&mut self) {
		(**self).negate()
	}
}

impl<T:SeaCoefficient> Neg for SeaFlectorPGA3D<T> {
	type Output = Self;

	fn neg(self) -> Self::Output {
		self.negated()
	}
}
impl<T:SeaCoefficient> Neg for &SeaFlectorPGA3D<T> {
	type Output = SeaFlectorPGA3D<T>;

	fn neg(self) -> Self::Output {
		self.negated()
	}
}