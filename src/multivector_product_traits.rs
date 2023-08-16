use crate::coefficient::SeaCoefficient;

pub trait GeometricProduct<Rhs=Self> {
	type Output;

	fn geometric_product(&self, rhs: Rhs) -> Self::Output;
}

impl<T: SeaCoefficient> GeometricProduct for T {
	type Output = Self;

	fn geometric_product(&self, rhs: Self) -> Self::Output {
		*self * rhs
	}
}