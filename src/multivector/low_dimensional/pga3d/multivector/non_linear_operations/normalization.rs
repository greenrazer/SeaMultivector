use super::*;

impl<T:SeaCoefficient> MultivectorNorm<T> for SeaMultivectorPGA3D<T> {
	fn scalar_norm(&self) -> T {
		self.0.scalar_norm(&Self::METRIC)
	}

	fn scalar_normalized(&self) -> Self {
		Self(self.0.scalar_normalized(&Self::METRIC))
	}
}
impl<T:SeaCoefficient> MultivectorNormMut for SeaMultivectorPGA3D<T> {
	fn scalar_normalize(&mut self) {
		self.0.scalar_normalize(&Self::METRIC)
	}
}