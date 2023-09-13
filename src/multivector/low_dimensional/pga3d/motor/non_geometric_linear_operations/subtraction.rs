use super::*;

impl<T:SeaCoefficient> SubAssign<Self> for SeaMotorPGA3D<T> {
	fn sub_assign(&mut self, rhs: Self) {
		**self -= &*rhs;
	}
}
impl<T:SeaCoefficient> SubAssign<&Self> for SeaMotorPGA3D<T> {
	fn sub_assign(&mut self, rhs: &Self) {
		**self -= &**rhs;
	}
}
impl<T:SeaCoefficient> Sub<Self> for SeaMotorPGA3D<T> {
	type Output = SeaMotorPGA3D<T>;
	fn sub(self, rhs: Self) -> Self::Output {
		Self(&*self - &*rhs)
	}
}
impl<T:SeaCoefficient> Sub<&Self> for SeaMotorPGA3D<T> {
	type Output = SeaMotorPGA3D<T>;
	fn sub(self, rhs: &Self) -> Self::Output {
		Self(&*self - &**rhs)
	}
}
impl<T:SeaCoefficient> Sub<SeaMotorPGA3D<T>> for &SeaMotorPGA3D<T> {
	type Output = SeaMotorPGA3D<T>;
	fn sub(self, rhs: SeaMotorPGA3D<T>) -> Self::Output {
		SeaMotorPGA3D::<T>(&**self - &*rhs)
	}
}
impl<T:SeaCoefficient> Sub<&SeaMotorPGA3D<T>> for &SeaMotorPGA3D<T> {
	type Output = SeaMotorPGA3D<T>;
	fn sub(self, rhs: &SeaMotorPGA3D<T>) -> Self::Output {
		SeaMotorPGA3D::<T>(&**self - &**rhs)
	}
}

//Multivector

impl<T:SeaCoefficient> SubAssign<SeaMultivectorPGA3D<T>> for SeaMotorPGA3D<T> {
	fn sub_assign(&mut self, rhs: SeaMultivectorPGA3D<T>) {
		**self -= rhs;
	}
}
impl<T:SeaCoefficient> SubAssign<&SeaMultivectorPGA3D<T>> for SeaMotorPGA3D<T> {
	fn sub_assign(&mut self, rhs: &SeaMultivectorPGA3D<T>) {
		**self -= &*rhs;
	}
}
impl<T:SeaCoefficient> Sub<SeaMultivectorPGA3D<T>> for SeaMotorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn sub(self, rhs: SeaMultivectorPGA3D<T>) -> Self::Output {
		&*self - rhs
	}
}
impl<T:SeaCoefficient> Sub<&SeaMultivectorPGA3D<T>> for SeaMotorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn sub(self, rhs: &SeaMultivectorPGA3D<T>) -> Self::Output {
		&*self - &*rhs
	}
}
impl<T:SeaCoefficient> Sub<SeaMultivectorPGA3D<T>> for &SeaMotorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn sub(self, rhs: SeaMultivectorPGA3D<T>) -> Self::Output {
		&**self - rhs
	}
}
impl<T:SeaCoefficient> Sub<&SeaMultivectorPGA3D<T>> for &SeaMotorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn sub(self, rhs: &SeaMultivectorPGA3D<T>) -> Self::Output {
		&**self - &*rhs
	}
}
