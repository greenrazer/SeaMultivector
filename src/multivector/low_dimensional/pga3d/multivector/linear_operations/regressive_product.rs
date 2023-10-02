use super::*;

impl<T:SeaCoefficient> RegressiveProduct<Self> for SeaMultivectorPGA3D<T> {
	type Output = Self;

	fn regressive_product(&self, rhs: Self) -> Self::Output {
		Self(self.0.regressive_product(&rhs.0, &Self::METRIC))
	}
}
impl<T:SeaCoefficient> RegressiveProduct<&Self> for SeaMultivectorPGA3D<T> {
	type Output = Self;

	fn regressive_product(&self, rhs: &Self) -> Self::Output {
		Self(self.0.regressive_product(&rhs.0, &Self::METRIC))
	}
}
impl<T:SeaCoefficient> BitAnd<Self> for SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	
	fn bitand(self, rhs: Self) -> Self::Output {
		self.regressive_product(rhs)
	}
}
impl<T:SeaCoefficient> BitAnd<&Self> for SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	
	fn bitand(self, rhs: &Self) -> Self::Output {
		self.regressive_product(rhs)
	}
}
impl<T:SeaCoefficient> BitAnd<SeaMultivectorPGA3D<T>> for &SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	
	fn bitand(self, rhs: SeaMultivectorPGA3D<T>) -> Self::Output {
		self.regressive_product(rhs)
	}
}
impl<T:SeaCoefficient> BitAnd<Self> for &SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	
	fn bitand(self, rhs: Self) -> Self::Output {
		self.regressive_product(rhs)
	}
}

//Motor
impl<T:SeaCoefficient> RegressiveProduct<SeaMotorPGA3D<T>> for SeaMultivectorPGA3D<T> {
	type Output = Self;

	fn regressive_product(&self, rhs: SeaMotorPGA3D<T>) -> Self::Output {
		self.regressive_product(&*rhs)
	}
}
impl<T:SeaCoefficient> RegressiveProduct<&SeaMotorPGA3D<T>> for SeaMultivectorPGA3D<T> {
	type Output = Self;

	fn regressive_product(&self, rhs: &SeaMotorPGA3D<T>) -> Self::Output {
		self.regressive_product(&**rhs)
	}
}
impl<T:SeaCoefficient> BitAnd<SeaMotorPGA3D<T>> for SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	
	fn bitand(self, rhs: SeaMotorPGA3D<T>) -> Self::Output {
		self.regressive_product(rhs)
	}
}
impl<T:SeaCoefficient> BitAnd<&SeaMotorPGA3D<T>> for SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	
	fn bitand(self, rhs: &SeaMotorPGA3D<T>) -> Self::Output {
		self.regressive_product(rhs)
	}
}
impl<T:SeaCoefficient> BitAnd<SeaMotorPGA3D<T>> for &SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	
	fn bitand(self, rhs: SeaMotorPGA3D<T>) -> Self::Output {
		self.regressive_product(rhs)
	}
}
impl<T:SeaCoefficient> BitAnd<&SeaMotorPGA3D<T>> for &SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	
	fn bitand(self, rhs: &SeaMotorPGA3D<T>) -> Self::Output {
		self.regressive_product(rhs)
	}
}

//Flector
impl<T:SeaCoefficient> RegressiveProduct<SeaFlectorPGA3D<T>> for SeaMultivectorPGA3D<T> {
	type Output = Self;

	fn regressive_product(&self, rhs: SeaFlectorPGA3D<T>) -> Self::Output {
		self.regressive_product(&*rhs)
	}
}
impl<T:SeaCoefficient> RegressiveProduct<&SeaFlectorPGA3D<T>> for SeaMultivectorPGA3D<T> {
	type Output = Self;

	fn regressive_product(&self, rhs: &SeaFlectorPGA3D<T>) -> Self::Output {
		self.regressive_product(&**rhs)
	}
}
impl<T:SeaCoefficient> BitAnd<SeaFlectorPGA3D<T>> for SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	
	fn bitand(self, rhs: SeaFlectorPGA3D<T>) -> Self::Output {
		self.regressive_product(rhs)
	}
}
impl<T:SeaCoefficient> BitAnd<&SeaFlectorPGA3D<T>> for SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	
	fn bitand(self, rhs: &SeaFlectorPGA3D<T>) -> Self::Output {
		self.regressive_product(rhs)
	}
}
impl<T:SeaCoefficient> BitAnd<SeaFlectorPGA3D<T>> for &SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	
	fn bitand(self, rhs: SeaFlectorPGA3D<T>) -> Self::Output {
		self.regressive_product(rhs)
	}
}
impl<T:SeaCoefficient> BitAnd<&SeaFlectorPGA3D<T>> for &SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	
	fn bitand(self, rhs: &SeaFlectorPGA3D<T>) -> Self::Output {
		self.regressive_product(rhs)
	}
}

//Versor
impl<T:SeaCoefficient> RegressiveProduct<SeaVersorPGA3D<T>> for SeaMultivectorPGA3D<T> {
	type Output = Self;

	fn regressive_product(&self, rhs: SeaVersorPGA3D<T>) -> Self::Output {
		self.regressive_product(&*rhs)
	}
}
impl<T:SeaCoefficient> RegressiveProduct<&SeaVersorPGA3D<T>> for SeaMultivectorPGA3D<T> {
	type Output = Self;

	fn regressive_product(&self, rhs: &SeaVersorPGA3D<T>) -> Self::Output {
		self.regressive_product(&**rhs)
	}
}
impl<T:SeaCoefficient> BitAnd<SeaVersorPGA3D<T>> for SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	
	fn bitand(self, rhs: SeaVersorPGA3D<T>) -> Self::Output {
		self.regressive_product(rhs)
	}
}
impl<T:SeaCoefficient> BitAnd<&SeaVersorPGA3D<T>> for SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	
	fn bitand(self, rhs: &SeaVersorPGA3D<T>) -> Self::Output {
		self.regressive_product(rhs)
	}
}
impl<T:SeaCoefficient> BitAnd<SeaVersorPGA3D<T>> for &SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	
	fn bitand(self, rhs: SeaVersorPGA3D<T>) -> Self::Output {
		self.regressive_product(rhs)
	}
}
impl<T:SeaCoefficient> BitAnd<&SeaVersorPGA3D<T>> for &SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	
	fn bitand(self, rhs: &SeaVersorPGA3D<T>) -> Self::Output {
		self.regressive_product(rhs)
	}
}

