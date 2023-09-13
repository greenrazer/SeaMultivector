use super::*;

impl<T:SeaCoefficient> MultivectorConjugation for SeaMotorPGA3D<T> {
	fn reversion_conjugated(&self) -> Self {
		Self((**self).reversion_conjugated())
	}

	fn involution_conjugated(&self) -> Self {
		Self((**self).involution_conjugated())
	}

	fn clifford_conjugated(&self) -> Self {
		Self((**self).clifford_conjugated())
	}
}
impl<T:SeaCoefficient> MultivectorConjugationMut for SeaMotorPGA3D<T> {	
	fn reversion_conjugation(&mut self) {
		(**self).reversion_conjugation()
	}

	fn involution_conjugation(&mut self) {
		(**self).involution_conjugation()
	}

	fn clifford_conjugation(&mut self) {
		(**self).clifford_conjugation()
	}
}