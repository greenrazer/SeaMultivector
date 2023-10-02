use super::*;

impl<T:SeaCoefficient> AddAssign<T> for SeaFlectorPGA3D<T> {
	fn add_assign(&mut self, _rhs: T) {
		unimplemented!("Cannot add_assign scalar to flector");
	}
}

impl<T:SeaCoefficient> Add<T> for SeaFlectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	
	fn add(self, rhs: T) -> Self::Output {
		&*self + rhs
	}
}
impl<T:SeaCoefficient> Add<T> for &SeaFlectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	
	fn add(self, rhs: T) -> Self::Output {
		&(**self) + rhs
	}
}
