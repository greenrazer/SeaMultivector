use super::*;

impl<T:SeaCoefficient> SubAssign<T> for SeaFlectorPGA3D<T> {
	fn sub_assign(&mut self, _rhs: T) {
		unimplemented!("Cannot sub_assign scalar from flector");
	}
}
impl<T:SeaCoefficient> Sub<T> for SeaFlectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;

	fn sub(self, rhs: T) -> Self::Output {
		&*self - rhs
	}
	
}
impl<T:SeaCoefficient> Sub<T> for &SeaFlectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	
	fn sub(self, rhs: T) -> Self::Output {
		&(**self) - rhs
	}
}
