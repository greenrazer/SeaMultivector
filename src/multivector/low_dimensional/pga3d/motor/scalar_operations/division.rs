use super::*;

impl<T:SeaCoefficient> DivAssign<T> for SeaMotorPGA3D<T> {
	fn div_assign(&mut self, rhs: T) {
		(**self) /= rhs
	}
}
impl<T:SeaCoefficient> Div<T> for SeaMotorPGA3D<T> {
	type Output = SeaMotorPGA3D<T>;

	fn div(self, rhs: T) -> Self::Output {
		Self(&*self/rhs)
	}
}
impl<T:SeaCoefficient> Div<T> for &SeaMotorPGA3D<T> {
	type Output = SeaMotorPGA3D<T>;

	fn div(self, rhs: T) -> Self::Output {
		SeaMotorPGA3D::<T>(&(**self)/rhs)
	}
}
