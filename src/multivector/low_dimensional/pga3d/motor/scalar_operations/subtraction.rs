use super::*;

impl<T:SeaCoefficient> SubAssign<T> for SeaMotorPGA3D<T> {
	fn sub_assign(&mut self, rhs: T) {
		(**self) -= rhs
	}
}
impl<T:SeaCoefficient> Sub<T> for SeaMotorPGA3D<T> {
	type Output = SeaMotorPGA3D<T>;

	fn sub(self, rhs: T) -> Self::Output {
		Self(&*self - rhs)
	}
	
}
impl<T:SeaCoefficient> Sub<T> for &SeaMotorPGA3D<T> {
	type Output = SeaMotorPGA3D<T>;
	
	fn sub(self, rhs: T) -> Self::Output {
		SeaMotorPGA3D::<T>(&(**self) - rhs)
	}
}
