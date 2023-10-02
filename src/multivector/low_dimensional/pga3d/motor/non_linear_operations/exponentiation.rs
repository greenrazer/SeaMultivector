use super::*;

impl<T:SeaCoefficient> MultivectorExponentiation for SeaMotorPGA3D<T> {
	fn exp(&self) -> Self {
		Self((**self).exp())
	}
}