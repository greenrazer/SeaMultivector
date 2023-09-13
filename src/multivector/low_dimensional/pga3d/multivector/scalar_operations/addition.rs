use super::*;

impl<T:SeaCoefficient> AddAssign<T> for SeaMultivectorPGA3D<T> {
	fn add_assign(&mut self, rhs: T) {
		self.0 += rhs
	}
}
impl<T:SeaCoefficient> Add<T> for SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	
	fn add(self, rhs: T) -> Self::Output {
		Self(self.0 + rhs)
	}
}
impl<T:SeaCoefficient> Add<T> for &SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	
	fn add(self, rhs: T) -> Self::Output {
		SeaMultivectorPGA3D::<T>(&self.0 + rhs)
	}
}
