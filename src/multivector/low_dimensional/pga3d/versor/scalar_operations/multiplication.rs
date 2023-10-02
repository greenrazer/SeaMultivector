use super::*;


impl<T:SeaCoefficient> MulAssign<T> for SeaVersorPGA3D<T> {
	fn mul_assign(&mut self, rhs: T) {
		(**self) *= rhs
	}
}
impl<T:SeaCoefficient> Mul<T> for SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;

	fn mul(self, rhs: T) -> Self::Output {
		Self(&*self * rhs)
	}
}
impl<T:SeaCoefficient> Mul<T> for &SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;

	fn mul(self, rhs: T) -> Self::Output {
		SeaVersorPGA3D::<T>(&(**self) * rhs)
	}
}
