use super::*;

impl<T:SeaCoefficient> DivAssign<T> for SeaFlectorPGA3D<T> {
	fn div_assign(&mut self, rhs: T) {
		(**self) /= rhs
	}
}
impl<T:SeaCoefficient> Div<T> for SeaFlectorPGA3D<T> {
	type Output = SeaFlectorPGA3D<T>;

	fn div(self, rhs: T) -> Self::Output {
		Self(&*self/rhs)
	}
}
impl<T:SeaCoefficient> Div<T> for &SeaFlectorPGA3D<T> {
	type Output = SeaFlectorPGA3D<T>;

	fn div(self, rhs: T) -> Self::Output {
		SeaFlectorPGA3D::<T>(&(**self)/rhs)
	}
}
