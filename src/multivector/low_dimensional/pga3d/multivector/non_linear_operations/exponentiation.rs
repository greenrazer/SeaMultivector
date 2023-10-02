use super::*;

impl<T:SeaCoefficient> MultivectorExponentiation for SeaMultivectorPGA3D<T> {
	type ExponentialOutput = Self;
	fn exp(&self) -> Self {
		Self(self.0.exp(&Self::METRIC))
	}
}