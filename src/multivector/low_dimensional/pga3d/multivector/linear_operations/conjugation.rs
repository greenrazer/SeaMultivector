use super::*;

impl<T:SeaCoefficient> MultivectorConjugation for SeaMultivectorPGA3D<T> {
	fn reversion_conjugated(&self) -> Self {
		Self(self.0.reversion_conjugated())
	}

	fn involution_conjugated(&self) -> Self {
		Self(self.0.involution_conjugated())
	}

	fn clifford_conjugated(&self) -> Self {
		Self(self.0.clifford_conjugated())
	}

	fn negated(&self) -> Self {
		Self(self.0.negated())
	}
}
impl<T:SeaCoefficient> MultivectorConjugationMut for SeaMultivectorPGA3D<T> {
	fn reversion_conjugation(&mut self) {
		self.0.reversion_conjugation()
	}

	fn involution_conjugation(&mut self) {
		self.0.involution_conjugation()
	}

	fn clifford_conjugation(&mut self) {
		self.0.clifford_conjugation()
	}

	fn negate(&mut self) {
		self.0.negate()
	}
}

impl<T:SeaCoefficient> Neg for SeaMultivectorPGA3D<T> {
	type Output = Self;

	fn neg(self) -> Self::Output {
		Self(-&self.0)
	}
}
impl<T:SeaCoefficient> Neg for &SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;

	fn neg(self) -> Self::Output {
		SeaMultivectorPGA3D::<T>(-&self.0)
	}
}