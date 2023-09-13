use super::*;

impl<T:SeaCoefficient> SubAssign<T> for SeaMultivectorPGA3D<T> {
	fn sub_assign(&mut self, rhs: T) {
		self.0 -= rhs
	}
}
impl<T:SeaCoefficient> Sub<T> for SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;

	fn sub(self, rhs: T) -> Self::Output {
		Self(self.0 - rhs)
	}
	
}
impl<T:SeaCoefficient> Sub<T> for &SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	
	fn sub(self, rhs: T) -> Self::Output {
		SeaMultivectorPGA3D::<T>(&self.0 - rhs)
	}
}
