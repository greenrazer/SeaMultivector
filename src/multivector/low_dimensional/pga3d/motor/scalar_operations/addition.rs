use super::*;

impl<T:SeaCoefficient> AddAssign<T> for SeaMotorPGA3D<T> {
	fn add_assign(&mut self, rhs: T) {
		(**self) += rhs
	}
}
impl<T:SeaCoefficient> Add<T> for SeaMotorPGA3D<T> {
	type Output = SeaMotorPGA3D<T>;
	
	fn add(self, rhs: T) -> Self::Output {
		Self(&*self + rhs)
	}
}
impl<T:SeaCoefficient> Add<T> for &SeaMotorPGA3D<T> {
	type Output = SeaMotorPGA3D<T>;
	
	fn add(self, rhs: T) -> Self::Output {
		SeaMotorPGA3D::<T>(&(**self) + rhs)
	}
}
