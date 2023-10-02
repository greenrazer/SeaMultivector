use super::*;

impl<T:SeaCoefficient> MultivectorExponentiation for SeaVersorPGA3D<T> {
	type ExponentialOutput = SeaVersorPGA3D<T>;
	fn exp(&self) -> Self::ExponentialOutput {
		Self((**self).exp())
	}
}