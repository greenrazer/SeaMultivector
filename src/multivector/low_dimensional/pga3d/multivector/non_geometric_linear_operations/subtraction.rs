use super::*;

// Self
impl<T:SeaCoefficient> SubAssign<Self> for SeaMultivectorPGA3D<T> {
	fn sub_assign(&mut self, rhs: Self) {
		self.0 -= &rhs.0;
	}
}
impl<T:SeaCoefficient> SubAssign<&Self> for SeaMultivectorPGA3D<T> {
	fn sub_assign(&mut self, rhs: &Self) {
		self.0 -= &rhs.0;
	}
}


impl<T:SeaCoefficient> Sub<Self> for SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn sub(self, rhs: Self) -> Self::Output {
		Self(&self.0 - &rhs.0)
	}
}
impl<T:SeaCoefficient> Sub<&Self> for SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn sub(self, rhs: &Self) -> Self::Output {
		Self(&self.0 - &rhs.0)
	}
}
impl<T:SeaCoefficient> Sub<SeaMultivectorPGA3D<T>> for &SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn sub(self, rhs: SeaMultivectorPGA3D<T>) -> Self::Output {
		SeaMultivectorPGA3D::<T>(&self.0 - &rhs.0)
	}
}
impl<T:SeaCoefficient> Sub<&SeaMultivectorPGA3D<T>> for &SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn sub(self, rhs: &SeaMultivectorPGA3D<T>) -> Self::Output {
		SeaMultivectorPGA3D::<T>(&self.0 - &rhs.0)
	}
}


// Motor
impl<T:SeaCoefficient> SubAssign<SeaMotorPGA3D<T>> for SeaMultivectorPGA3D<T> {
	fn sub_assign(&mut self, rhs: SeaMotorPGA3D<T>) {
		*self -= &*rhs;
	}
}
impl<T:SeaCoefficient> SubAssign<&SeaMotorPGA3D<T>> for SeaMultivectorPGA3D<T> {
	fn sub_assign(&mut self, rhs: &SeaMotorPGA3D<T>) {
		*self -= &**rhs;
	}
}

impl<T:SeaCoefficient> Sub<SeaMotorPGA3D<T>> for SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn sub(self, rhs: SeaMotorPGA3D<T>) -> Self::Output {
		self - &*rhs
	}
}
impl<T:SeaCoefficient> Sub<&SeaMotorPGA3D<T>> for SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn sub(self, rhs: &SeaMotorPGA3D<T>) -> Self::Output {
		self - &**rhs
	}
}
impl<T:SeaCoefficient> Sub<SeaMotorPGA3D<T>> for &SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn sub(self, rhs: SeaMotorPGA3D<T>) -> Self::Output {
		self - &*rhs
	}
}
impl<T:SeaCoefficient> Sub<&SeaMotorPGA3D<T>> for &SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn sub(self, rhs: &SeaMotorPGA3D<T>) -> Self::Output {
		self - &**rhs
	}
}

// Flector
impl<T:SeaCoefficient> SubAssign<SeaFlectorPGA3D<T>> for SeaMultivectorPGA3D<T> {
	fn sub_assign(&mut self, rhs: SeaFlectorPGA3D<T>) {
		*self -= &*rhs;
	}
}
impl<T:SeaCoefficient> SubAssign<&SeaFlectorPGA3D<T>> for SeaMultivectorPGA3D<T> {
	fn sub_assign(&mut self, rhs: &SeaFlectorPGA3D<T>) {
		*self -= &**rhs;
	}
}

impl<T:SeaCoefficient> Sub<SeaFlectorPGA3D<T>> for SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn sub(self, rhs: SeaFlectorPGA3D<T>) -> Self::Output {
		self - &*rhs
	}
}
impl<T:SeaCoefficient> Sub<&SeaFlectorPGA3D<T>> for SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn sub(self, rhs: &SeaFlectorPGA3D<T>) -> Self::Output {
		self - &**rhs
	}
}
impl<T:SeaCoefficient> Sub<SeaFlectorPGA3D<T>> for &SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn sub(self, rhs: SeaFlectorPGA3D<T>) -> Self::Output {
		self - &*rhs
	}
}
impl<T:SeaCoefficient> Sub<&SeaFlectorPGA3D<T>> for &SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn sub(self, rhs: &SeaFlectorPGA3D<T>) -> Self::Output {
		self - &**rhs
	}
}

// Versor
impl<T:SeaCoefficient> SubAssign<SeaVersorPGA3D<T>> for SeaMultivectorPGA3D<T> {
	fn sub_assign(&mut self, rhs: SeaVersorPGA3D<T>) {
		*self -= &*rhs;
	}
}
impl<T:SeaCoefficient> SubAssign<&SeaVersorPGA3D<T>> for SeaMultivectorPGA3D<T> {
	fn sub_assign(&mut self, rhs: &SeaVersorPGA3D<T>) {
		*self -= &**rhs;
	}
}

impl<T:SeaCoefficient> Sub<SeaVersorPGA3D<T>> for SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn sub(self, rhs: SeaVersorPGA3D<T>) -> Self::Output {
		self - &*rhs
	}
}
impl<T:SeaCoefficient> Sub<&SeaVersorPGA3D<T>> for SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn sub(self, rhs: &SeaVersorPGA3D<T>) -> Self::Output {
		self - &**rhs
	}
}
impl<T:SeaCoefficient> Sub<SeaVersorPGA3D<T>> for &SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn sub(self, rhs: SeaVersorPGA3D<T>) -> Self::Output {
		self - &*rhs
	}
}
impl<T:SeaCoefficient> Sub<&SeaVersorPGA3D<T>> for &SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn sub(self, rhs: &SeaVersorPGA3D<T>) -> Self::Output {
		self - &**rhs
	}
}