use super::*;

impl<T:SeaCoefficient> SubAssign<T> for SeaVersorPGA3D<T> {
	fn sub_assign(&mut self, rhs: T) {
		(**self) -= rhs
	}
}
impl<T:SeaCoefficient> Sub<T> for SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;

	fn sub(self, rhs: T) -> Self::Output {
		Self(&*self - rhs)
	}
	
}
impl<T:SeaCoefficient> Sub<T> for &SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;
	
	fn sub(self, rhs: T) -> Self::Output {
		SeaVersorPGA3D::<T>(&(**self) - rhs)
	}
}
