use super::*;

impl<T:SeaCoefficient> SubAssign<Self> for SeaVersorPGA3D<T> {
	fn sub_assign(&mut self, rhs: Self) {
		**self -= &*rhs;
	}
}
impl<T:SeaCoefficient> SubAssign<&Self> for SeaVersorPGA3D<T> {
	fn sub_assign(&mut self, rhs: &Self) {
		**self -= &**rhs;
	}
}
impl<T:SeaCoefficient> Sub<Self> for SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;
	fn sub(self, rhs: Self) -> Self::Output {
		Self(&*self - &*rhs)
	}
}
impl<T:SeaCoefficient> Sub<&Self> for SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;
	fn sub(self, rhs: &Self) -> Self::Output {
		Self(&*self - &**rhs)
	}
}
impl<T:SeaCoefficient> Sub<SeaVersorPGA3D<T>> for &SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;
	fn sub(self, rhs: SeaVersorPGA3D<T>) -> Self::Output {
		SeaVersorPGA3D::<T>(&**self - &*rhs)
	}
}
impl<T:SeaCoefficient> Sub<&SeaVersorPGA3D<T>> for &SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;
	fn sub(self, rhs: &SeaVersorPGA3D<T>) -> Self::Output {
		SeaVersorPGA3D::<T>(&**self - &**rhs)
	}
}

//Multivector
impl<T:SeaCoefficient> SubAssign<SeaMultivectorPGA3D<T>> for SeaVersorPGA3D<T> {
	fn sub_assign(&mut self, rhs: SeaMultivectorPGA3D<T>) {
		**self -= rhs;
	}
}
impl<T:SeaCoefficient> SubAssign<&SeaMultivectorPGA3D<T>> for SeaVersorPGA3D<T> {
	fn sub_assign(&mut self, rhs: &SeaMultivectorPGA3D<T>) {
		**self -= &*rhs;
	}
}
impl<T:SeaCoefficient> Sub<SeaMultivectorPGA3D<T>> for SeaVersorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn sub(self, rhs: SeaMultivectorPGA3D<T>) -> Self::Output {
		&*self - rhs
	}
}
impl<T:SeaCoefficient> Sub<&SeaMultivectorPGA3D<T>> for SeaVersorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn sub(self, rhs: &SeaMultivectorPGA3D<T>) -> Self::Output {
		&*self - &*rhs
	}
}
impl<T:SeaCoefficient> Sub<SeaMultivectorPGA3D<T>> for &SeaVersorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn sub(self, rhs: SeaMultivectorPGA3D<T>) -> Self::Output {
		&**self - rhs
	}
}
impl<T:SeaCoefficient> Sub<&SeaMultivectorPGA3D<T>> for &SeaVersorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn sub(self, rhs: &SeaMultivectorPGA3D<T>) -> Self::Output {
		&**self - &*rhs
	}
}

//Motor
impl<T:SeaCoefficient> SubAssign<SeaMotorPGA3D<T>> for SeaVersorPGA3D<T> {
	fn sub_assign(&mut self, rhs: SeaMotorPGA3D<T>) {
		**self -= rhs;
	}
}
impl<T:SeaCoefficient> SubAssign<&SeaMotorPGA3D<T>> for SeaVersorPGA3D<T> {
	fn sub_assign(&mut self, rhs: &SeaMotorPGA3D<T>) {
		**self -= &*rhs;
	}
}
impl<T:SeaCoefficient> Sub<SeaMotorPGA3D<T>> for SeaVersorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn sub(self, rhs: SeaMotorPGA3D<T>) -> Self::Output {
		&*self - rhs
	}
}
impl<T:SeaCoefficient> Sub<&SeaMotorPGA3D<T>> for SeaVersorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn sub(self, rhs: &SeaMotorPGA3D<T>) -> Self::Output {
		&*self - &*rhs
	}
}
impl<T:SeaCoefficient> Sub<SeaMotorPGA3D<T>> for &SeaVersorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn sub(self, rhs: SeaMotorPGA3D<T>) -> Self::Output {
		&**self - rhs
	}
}
impl<T:SeaCoefficient> Sub<&SeaMotorPGA3D<T>> for &SeaVersorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn sub(self, rhs: &SeaMotorPGA3D<T>) -> Self::Output {
		&**self - &*rhs
	}
}


//Flector
impl<T:SeaCoefficient> SubAssign<SeaFlectorPGA3D<T>> for SeaVersorPGA3D<T> {
	fn sub_assign(&mut self, rhs: SeaFlectorPGA3D<T>) {
		**self -= rhs;
	}
}
impl<T:SeaCoefficient> SubAssign<&SeaFlectorPGA3D<T>> for SeaVersorPGA3D<T> {
	fn sub_assign(&mut self, rhs: &SeaFlectorPGA3D<T>) {
		**self -= &*rhs;
	}
}
impl<T:SeaCoefficient> Sub<SeaFlectorPGA3D<T>> for SeaVersorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn sub(self, rhs: SeaFlectorPGA3D<T>) -> Self::Output {
		&*self - rhs
	}
}
impl<T:SeaCoefficient> Sub<&SeaFlectorPGA3D<T>> for SeaVersorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn sub(self, rhs: &SeaFlectorPGA3D<T>) -> Self::Output {
		&*self - &*rhs
	}
}
impl<T:SeaCoefficient> Sub<SeaFlectorPGA3D<T>> for &SeaVersorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn sub(self, rhs: SeaFlectorPGA3D<T>) -> Self::Output {
		&**self - rhs
	}
}
impl<T:SeaCoefficient> Sub<&SeaFlectorPGA3D<T>> for &SeaVersorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn sub(self, rhs: &SeaFlectorPGA3D<T>) -> Self::Output {
		&**self - &*rhs
	}
}

