use super::*;

impl<T:SeaCoefficient> MultivectorExponentiation for SeaFlectorPGA3D<T> {
	type ExponentialOutput = SeaMultivectorPGA3D<T>;
	fn exp(&self) -> Self::ExponentialOutput {
		(**self).exp()
	}
}