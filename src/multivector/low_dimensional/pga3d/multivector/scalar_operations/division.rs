use super::*;

impl<T:SeaCoefficient> DivAssign<T> for SeaMultivectorPGA3D<T> {
	fn div_assign(&mut self, rhs: T) {
		self.0 /= rhs
	}
}
impl<T:SeaCoefficient> Div<T> for SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;

	fn div(self, rhs: T) -> Self::Output {
		Self(self.0/rhs)
	}
}
impl<T:SeaCoefficient> Div<T> for &SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;

	fn div(self, rhs: T) -> Self::Output {
		SeaMultivectorPGA3D::<T>(&self.0/rhs)
	}
}
