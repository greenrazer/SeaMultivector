use super::*;

impl<T:SeaCoefficient> AddAssign<Self> for SeaMotorPGA3D<T> {
	fn add_assign(&mut self, rhs: Self) {
		**self += &*rhs;
	}
}
impl<T:SeaCoefficient> AddAssign<&Self> for SeaMotorPGA3D<T> {
	fn add_assign(&mut self, rhs: &Self) {
		**self += &**rhs;
	}
}
impl<T:SeaCoefficient> Add<Self> for SeaMotorPGA3D<T> {
	type Output = SeaMotorPGA3D<T>;
	fn add(self, rhs: Self) -> Self::Output {
		Self(&*self + &*rhs)
	}
}
impl<T:SeaCoefficient> Add<&Self> for SeaMotorPGA3D<T> {
	type Output = SeaMotorPGA3D<T>;
	fn add(self, rhs: &Self) -> Self::Output {
		Self(&*self + &**rhs)
	}
}
impl<T:SeaCoefficient> Add<SeaMotorPGA3D<T>> for &SeaMotorPGA3D<T> {
	type Output = SeaMotorPGA3D<T>;
	fn add(self, rhs: SeaMotorPGA3D<T>) -> Self::Output {
		SeaMotorPGA3D::<T>(&**self + &*rhs)
	}
}
impl<T:SeaCoefficient> Add<&SeaMotorPGA3D<T>> for &SeaMotorPGA3D<T> {
	type Output = SeaMotorPGA3D<T>;
	fn add(self, rhs: &SeaMotorPGA3D<T>) -> Self::Output {
		SeaMotorPGA3D::<T>(&**self + &**rhs)
	}
}

// Multivector
impl<T:SeaCoefficient> AddAssign<SeaMultivectorPGA3D<T>> for SeaMotorPGA3D<T> {
	fn add_assign(&mut self, rhs: SeaMultivectorPGA3D<T>) {
		**self += rhs;
	}
}
impl<T:SeaCoefficient> AddAssign<&SeaMultivectorPGA3D<T>> for SeaMotorPGA3D<T> {
	fn add_assign(&mut self, rhs: &SeaMultivectorPGA3D<T>) {
		**self += rhs;
	}
}
impl<T:SeaCoefficient> Add<SeaMultivectorPGA3D<T>> for SeaMotorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn add(self, rhs: SeaMultivectorPGA3D<T>) -> Self::Output {
		&*self + rhs
	}
}
impl<T:SeaCoefficient> Add<&SeaMultivectorPGA3D<T>> for SeaMotorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn add(self, rhs: &SeaMultivectorPGA3D<T>) -> Self::Output {
		&*self + &*rhs
	}
}
impl<T:SeaCoefficient> Add<SeaMultivectorPGA3D<T>> for &SeaMotorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn add(self, rhs: SeaMultivectorPGA3D<T>) -> Self::Output {
		&**self + rhs
	}
}
impl<T:SeaCoefficient> Add<&SeaMultivectorPGA3D<T>> for &SeaMotorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn add(self, rhs: &SeaMultivectorPGA3D<T>) -> Self::Output {
		&**self + &*rhs
	}
}

// Flector
impl<T:SeaCoefficient> AddAssign<SeaFlectorPGA3D<T>> for SeaMotorPGA3D<T> {
	fn add_assign(&mut self, _rhs: SeaFlectorPGA3D<T>) {
		unimplemented!("cannot add_assign flector to motor")
	}
}
impl<T:SeaCoefficient> AddAssign<&SeaFlectorPGA3D<T>> for SeaMotorPGA3D<T> {
	fn add_assign(&mut self, _rhs: &SeaFlectorPGA3D<T>) {
		unimplemented!("cannot add_assign flector to motor")
	}
}
impl<T:SeaCoefficient> Add<SeaFlectorPGA3D<T>> for SeaMotorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn add(self, rhs: SeaFlectorPGA3D<T>) -> Self::Output {
		&*self + rhs
	}
}
impl<T:SeaCoefficient> Add<&SeaFlectorPGA3D<T>> for SeaMotorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn add(self, rhs: &SeaFlectorPGA3D<T>) -> Self::Output {
		&*self + &*rhs
	}
}
impl<T:SeaCoefficient> Add<SeaFlectorPGA3D<T>> for &SeaMotorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn add(self, rhs: SeaFlectorPGA3D<T>) -> Self::Output {
		&**self + rhs
	}
}
impl<T:SeaCoefficient> Add<&SeaFlectorPGA3D<T>> for &SeaMotorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn add(self, rhs: &SeaFlectorPGA3D<T>) -> Self::Output {
		&**self + &*rhs
	}
}

// Versor
impl<T:SeaCoefficient> AddAssign<SeaVersorPGA3D<T>> for SeaMotorPGA3D<T> {
	fn add_assign(&mut self, _rhs: SeaVersorPGA3D<T>) {
		unimplemented!("cannot add_assign versor to motor")
	}
}
impl<T:SeaCoefficient> AddAssign<&SeaVersorPGA3D<T>> for SeaMotorPGA3D<T> {
	fn add_assign(&mut self, _rhs: &SeaVersorPGA3D<T>) {
		unimplemented!("cannot add_assign versor to motor")
	}
}
impl<T:SeaCoefficient> Add<SeaVersorPGA3D<T>> for SeaMotorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn add(self, rhs: SeaVersorPGA3D<T>) -> Self::Output {
		&*self + rhs
	}
}
impl<T:SeaCoefficient> Add<&SeaVersorPGA3D<T>> for SeaMotorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn add(self, rhs: &SeaVersorPGA3D<T>) -> Self::Output {
		&*self + &*rhs
	}
}
impl<T:SeaCoefficient> Add<SeaVersorPGA3D<T>> for &SeaMotorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn add(self, rhs: SeaVersorPGA3D<T>) -> Self::Output {
		&**self + rhs
	}
}
impl<T:SeaCoefficient> Add<&SeaVersorPGA3D<T>> for &SeaMotorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn add(self, rhs: &SeaVersorPGA3D<T>) -> Self::Output {
		&**self + &*rhs
	}
}
