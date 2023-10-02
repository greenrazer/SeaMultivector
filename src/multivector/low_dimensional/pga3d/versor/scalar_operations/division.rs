use super::*;

impl<T:SeaCoefficient> DivAssign<T> for SeaVersorPGA3D<T> {
	fn div_assign(&mut self, rhs: T) {
		(**self) /= rhs
	}
}
impl<T:SeaCoefficient> Div<T> for SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;

	fn div(self, rhs: T) -> Self::Output {
		Self(&*self/rhs)
	}
}
impl<T:SeaCoefficient> Div<T> for &SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;

	fn div(self, rhs: T) -> Self::Output {
		SeaVersorPGA3D::<T>(&(**self)/rhs)
	}
}
