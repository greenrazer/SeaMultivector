use super::*;

impl<T:SeaCoefficient> AddAssign<T> for SeaVersorPGA3D<T> {
	fn add_assign(&mut self, rhs: T) {
		(**self) += rhs
	}
}

impl<T:SeaCoefficient> Add<T> for SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;
	
	fn add(self, rhs: T) -> Self::Output {
		Self(&*self + rhs)
	}
}
impl<T:SeaCoefficient> Add<T> for &SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;
	
	fn add(self, rhs: T) -> Self::Output {
		SeaVersorPGA3D::<T>(&(**self) + rhs)
	}
}
