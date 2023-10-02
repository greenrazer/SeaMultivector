use super::*;

impl<T:SeaCoefficient> MultivectorExponentiation for SeaMotorPGA3D<T> {
	type ExponentialOutput = Self;
	fn exp(&self) -> Self {
		Self((**self).exp())
	}
}