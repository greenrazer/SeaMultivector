use super::*;


impl<T:SeaCoefficient> MulAssign<T> for SeaMotorPGA3D<T> {
	fn mul_assign(&mut self, rhs: T) {
		(**self) *= rhs
	}
}
impl<T:SeaCoefficient> Mul<T> for SeaMotorPGA3D<T> {
	type Output = SeaMotorPGA3D<T>;

	fn mul(self, rhs: T) -> Self::Output {
		Self(&*self * rhs)
	}
}
impl<T:SeaCoefficient> Mul<T> for &SeaMotorPGA3D<T> {
	type Output = SeaMotorPGA3D<T>;

	fn mul(self, rhs: T) -> Self::Output {
		SeaMotorPGA3D::<T>(&(**self) * rhs)
	}
}
