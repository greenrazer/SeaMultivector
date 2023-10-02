use super::*;

impl<T:SeaCoefficient> MultivectorNorm<T> for SeaFlectorPGA3D<T> {
	fn scalar_norm(&self) -> T {
		(**self).scalar_norm()
	}

	fn scalar_normalized(&self) -> Self {
		Self((**self).scalar_normalized())
	}
}
impl<T:SeaCoefficient> MultivectorNormMut for SeaFlectorPGA3D<T> {
	fn scalar_normalize(&mut self) {
		(**self).scalar_normalize()
	}
}