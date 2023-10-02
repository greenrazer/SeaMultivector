use super::*;

// Self
impl<T:SeaCoefficient> AddAssign<Self> for SeaMultivectorPGA3D<T> {
	fn add_assign(&mut self, rhs: Self) {
		self.0 += &rhs.0;
	}
}
impl<T:SeaCoefficient> AddAssign<&Self> for SeaMultivectorPGA3D<T> {
	fn add_assign(&mut self, rhs: &Self) {
		self.0 += &rhs.0;
	}
}


impl<T:SeaCoefficient> Add<Self> for SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn add(self, rhs: Self) -> Self::Output {
		Self(&self.0 + &rhs.0)
	}
}
impl<T:SeaCoefficient> Add<&Self> for SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn add(self, rhs: &Self) -> Self::Output {
		Self(&self.0 + &rhs.0)
	}
}
impl<T:SeaCoefficient> Add<SeaMultivectorPGA3D<T>> for &SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn add(self, rhs: SeaMultivectorPGA3D<T>) -> Self::Output {
		SeaMultivectorPGA3D::<T>(&self.0 + &rhs.0)
	}
}
impl<T:SeaCoefficient> Add<&SeaMultivectorPGA3D<T>> for &SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn add(self, rhs: &SeaMultivectorPGA3D<T>) -> Self::Output {
		SeaMultivectorPGA3D::<T>(&self.0 + &rhs.0)
	}
}

// Motor
impl<T:SeaCoefficient> AddAssign<SeaMotorPGA3D<T>> for SeaMultivectorPGA3D<T> {
	fn add_assign(&mut self, rhs: SeaMotorPGA3D<T>) {
		*self += &*rhs;
	}
}
impl<T:SeaCoefficient> AddAssign<&SeaMotorPGA3D<T>> for SeaMultivectorPGA3D<T> {
	fn add_assign(&mut self, rhs: &SeaMotorPGA3D<T>) {
		*self += &**rhs;
	}
}

impl<T:SeaCoefficient> Add<SeaMotorPGA3D<T>> for SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn add(self, rhs: SeaMotorPGA3D<T>) -> Self::Output {
		self + &*rhs
	}
}
impl<T:SeaCoefficient> Add<&SeaMotorPGA3D<T>> for SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn add(self, rhs: &SeaMotorPGA3D<T>) -> Self::Output {
		self + &**rhs
	}
}
impl<T:SeaCoefficient> Add<SeaMotorPGA3D<T>> for &SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn add(self, rhs: SeaMotorPGA3D<T>) -> Self::Output {
		self + &*rhs
	}
}
impl<T:SeaCoefficient> Add<&SeaMotorPGA3D<T>> for &SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn add(self, rhs: &SeaMotorPGA3D<T>) -> Self::Output {
		self + &**rhs
	}
}

// Flector
impl<T:SeaCoefficient> AddAssign<SeaFlectorPGA3D<T>> for SeaMultivectorPGA3D<T> {
	fn add_assign(&mut self, rhs: SeaFlectorPGA3D<T>) {
		*self += &*rhs;
	}
}
impl<T:SeaCoefficient> AddAssign<&SeaFlectorPGA3D<T>> for SeaMultivectorPGA3D<T> {
	fn add_assign(&mut self, rhs: &SeaFlectorPGA3D<T>) {
		*self += &**rhs;
	}
}

impl<T:SeaCoefficient> Add<SeaFlectorPGA3D<T>> for SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn add(self, rhs: SeaFlectorPGA3D<T>) -> Self::Output {
		self + &*rhs
	}
}
impl<T:SeaCoefficient> Add<&SeaFlectorPGA3D<T>> for SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn add(self, rhs: &SeaFlectorPGA3D<T>) -> Self::Output {
		self + &**rhs
	}
}
impl<T:SeaCoefficient> Add<SeaFlectorPGA3D<T>> for &SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn add(self, rhs: SeaFlectorPGA3D<T>) -> Self::Output {
		self + &*rhs
	}
}
impl<T:SeaCoefficient> Add<&SeaFlectorPGA3D<T>> for &SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn add(self, rhs: &SeaFlectorPGA3D<T>) -> Self::Output {
		self + &**rhs
	}
}

// Versor
impl<T:SeaCoefficient> AddAssign<SeaVersorPGA3D<T>> for SeaMultivectorPGA3D<T> {
	fn add_assign(&mut self, rhs: SeaVersorPGA3D<T>) {
		*self += &*rhs;
	}
}
impl<T:SeaCoefficient> AddAssign<&SeaVersorPGA3D<T>> for SeaMultivectorPGA3D<T> {
	fn add_assign(&mut self, rhs: &SeaVersorPGA3D<T>) {
		*self += &**rhs;
	}
}

impl<T:SeaCoefficient> Add<SeaVersorPGA3D<T>> for SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn add(self, rhs: SeaVersorPGA3D<T>) -> Self::Output {
		self + &*rhs
	}
}
impl<T:SeaCoefficient> Add<&SeaVersorPGA3D<T>> for SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn add(self, rhs: &SeaVersorPGA3D<T>) -> Self::Output {
		self + &**rhs
	}
}
impl<T:SeaCoefficient> Add<SeaVersorPGA3D<T>> for &SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn add(self, rhs: SeaVersorPGA3D<T>) -> Self::Output {
		self + &*rhs
	}
}
impl<T:SeaCoefficient> Add<&SeaVersorPGA3D<T>> for &SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn add(self, rhs: &SeaVersorPGA3D<T>) -> Self::Output {
		self + &**rhs
	}
}