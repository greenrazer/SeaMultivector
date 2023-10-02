use super::*;


impl<T:SeaCoefficient> MulAssign<T> for SeaFlectorPGA3D<T> {
	fn mul_assign(&mut self, rhs: T) {
		(**self) *= rhs
	}
}
impl<T:SeaCoefficient> Mul<T> for SeaFlectorPGA3D<T> {
	type Output = SeaFlectorPGA3D<T>;

	fn mul(self, rhs: T) -> Self::Output {
		Self(&*self * rhs)
	}
}
impl<T:SeaCoefficient> Mul<T> for &SeaFlectorPGA3D<T> {
	type Output = SeaFlectorPGA3D<T>;

	fn mul(self, rhs: T) -> Self::Output {
		SeaFlectorPGA3D::<T>(&(**self) * rhs)
	}
}
