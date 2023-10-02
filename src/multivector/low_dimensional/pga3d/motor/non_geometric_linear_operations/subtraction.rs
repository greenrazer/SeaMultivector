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

// Flector
impl<T:SeaCoefficient> SubAssign<SeaFlectorPGA3D<T>> for SeaMotorPGA3D<T> {
	fn sub_assign(&mut self, _rhs: SeaFlectorPGA3D<T>) {
		unimplemented!("cannot sub_assign flector from motor")
	}
}
impl<T:SeaCoefficient> SubAssign<&SeaFlectorPGA3D<T>> for SeaMotorPGA3D<T> {
	fn sub_assign(&mut self, _rhs: &SeaFlectorPGA3D<T>) {
		unimplemented!("cannot sub_assign flector from motor")
	}
}
impl<T:SeaCoefficient> Sub<SeaFlectorPGA3D<T>> for SeaMotorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn sub(self, rhs: SeaFlectorPGA3D<T>) -> Self::Output {
		&*self - rhs
	}
}
impl<T:SeaCoefficient> Sub<&SeaFlectorPGA3D<T>> for SeaMotorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn sub(self, rhs: &SeaFlectorPGA3D<T>) -> Self::Output {
		&*self - &*rhs
	}
}
impl<T:SeaCoefficient> Sub<SeaFlectorPGA3D<T>> for &SeaMotorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn sub(self, rhs: SeaFlectorPGA3D<T>) -> Self::Output {
		&**self - rhs
	}
}
impl<T:SeaCoefficient> Sub<&SeaFlectorPGA3D<T>> for &SeaMotorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn sub(self, rhs: &SeaFlectorPGA3D<T>) -> Self::Output {
		&**self - &*rhs
	}
}

// Versor
impl<T:SeaCoefficient> SubAssign<SeaVersorPGA3D<T>> for SeaMotorPGA3D<T> {
	fn sub_assign(&mut self, _rhs: SeaVersorPGA3D<T>) {
		unimplemented!("cannot sub_assign versor from motor")
	}
}
impl<T:SeaCoefficient> SubAssign<&SeaVersorPGA3D<T>> for SeaMotorPGA3D<T> {
	fn sub_assign(&mut self, _rhs: &SeaVersorPGA3D<T>) {
		unimplemented!("cannot sub_assign versor from motor")
	}
}
impl<T:SeaCoefficient> Sub<SeaVersorPGA3D<T>> for SeaMotorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn sub(self, rhs: SeaVersorPGA3D<T>) -> Self::Output {
		&*self - rhs
	}
}
impl<T:SeaCoefficient> Sub<&SeaVersorPGA3D<T>> for SeaMotorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn sub(self, rhs: &SeaVersorPGA3D<T>) -> Self::Output {
		&*self - &*rhs
	}
}
impl<T:SeaCoefficient> Sub<SeaVersorPGA3D<T>> for &SeaMotorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn sub(self, rhs: SeaVersorPGA3D<T>) -> Self::Output {
		&**self - rhs
	}
}
impl<T:SeaCoefficient> Sub<&SeaVersorPGA3D<T>> for &SeaMotorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn sub(self, rhs: &SeaVersorPGA3D<T>) -> Self::Output {
		&**self - &*rhs
	}
}