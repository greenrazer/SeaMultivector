use super::*;

impl<T:SeaCoefficient> MulAssign<T> for SeaMultivectorPGA3D<T> {
	fn mul_assign(&mut self, rhs: T) {
		self.0 *= rhs
	}
}
impl<T:SeaCoefficient> Mul<T> for SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;

	fn mul(self, rhs: T) -> Self::Output {
		Self(self.0 * rhs)
	}
}
impl<T:SeaCoefficient> Mul<T> for &SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;

	fn mul(self, rhs: T) -> Self::Output {
		SeaMultivectorPGA3D::<T>(&self.0 * rhs)
	}
}
