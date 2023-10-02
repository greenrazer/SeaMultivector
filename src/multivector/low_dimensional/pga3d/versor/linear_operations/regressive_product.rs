use super::*;

impl<T:SeaCoefficient> RegressiveProduct<Self> for SeaVersorPGA3D<T> {
	type Output = Self;

	fn regressive_product(&self, rhs: Self) -> Self::Output {
		Self((**self).regressive_product(&*rhs))
	}
}
impl<T:SeaCoefficient> RegressiveProduct<&Self> for SeaVersorPGA3D<T> {
	type Output = Self;

	fn regressive_product(&self, rhs: &Self) -> Self::Output {
		Self((**self).regressive_product(&**rhs))
	}
}
impl<T:SeaCoefficient> BitAnd<Self> for SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;
	
	fn bitand(self, rhs: Self) -> Self::Output {
		self.regressive_product(rhs)
	}
}
impl<T:SeaCoefficient> BitAnd<&Self> for SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;
	
	fn bitand(self, rhs: &Self) -> Self::Output {
		self.regressive_product(rhs)
	}
}
impl<T:SeaCoefficient> BitAnd<SeaVersorPGA3D<T>> for &SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;
	
	fn bitand(self, rhs: SeaVersorPGA3D<T>) -> Self::Output {
		self.regressive_product(rhs)
	}
}
impl<T:SeaCoefficient> BitAnd<Self> for &SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;
	
	fn bitand(self, rhs: Self) -> Self::Output {
		self.regressive_product(rhs)
	}
}

// Multivector
impl<T:SeaCoefficient> RegressiveProduct<SeaMultivectorPGA3D<T>> for SeaVersorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;

	fn regressive_product(&self, rhs: SeaMultivectorPGA3D<T>) -> Self::Output {
		(**self).regressive_product(rhs)
	}
}
impl<T:SeaCoefficient> RegressiveProduct<&SeaMultivectorPGA3D<T>> for SeaVersorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;

	fn regressive_product(&self, rhs: &SeaMultivectorPGA3D<T>) -> Self::Output {
		(**self).regressive_product(rhs)
	}
}
impl<T:SeaCoefficient> BitAnd<SeaMultivectorPGA3D<T>> for SeaVersorPGA3D<T> {
	type Output =SeaMultivectorPGA3D<T>;
	
	fn bitand(self, rhs: SeaMultivectorPGA3D<T>) -> Self::Output {
		self.regressive_product(rhs)
	}
}
impl<T:SeaCoefficient> BitAnd<&SeaMultivectorPGA3D<T>> for SeaVersorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	
	fn bitand(self, rhs: &SeaMultivectorPGA3D<T>) -> Self::Output {
		self.regressive_product(rhs)
	}
}
impl<T:SeaCoefficient> BitAnd<SeaMultivectorPGA3D<T>> for &SeaVersorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	
	fn bitand(self, rhs: SeaMultivectorPGA3D<T>) -> Self::Output {
		self.regressive_product(rhs)
	}
}
impl<T:SeaCoefficient> BitAnd<&SeaMultivectorPGA3D<T>> for &SeaVersorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	
	fn bitand(self, rhs: &SeaMultivectorPGA3D<T>) -> Self::Output {
		self.regressive_product(rhs)
	}
}

// Motor
impl<T:SeaCoefficient> RegressiveProduct<SeaMotorPGA3D<T>> for SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;

	fn regressive_product(&self, rhs: SeaMotorPGA3D<T>) -> Self::Output {
		SeaVersorPGA3D::<T>((**self).regressive_product(rhs))
	}
}
impl<T:SeaCoefficient> RegressiveProduct<&SeaMotorPGA3D<T>> for SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;

	fn regressive_product(&self, rhs: &SeaMotorPGA3D<T>) -> Self::Output {
		SeaVersorPGA3D::<T>((**self).regressive_product(rhs))
	}
}
impl<T:SeaCoefficient> BitAnd<SeaMotorPGA3D<T>> for SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;
	fn bitand(self, rhs: SeaMotorPGA3D<T>) -> Self::Output {
		self.regressive_product(rhs)
	}
}
impl<T:SeaCoefficient> BitAnd<&SeaMotorPGA3D<T>> for SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;
	fn bitand(self, rhs: &SeaMotorPGA3D<T>) -> Self::Output {
		self.regressive_product(rhs)
	}
}
impl<T:SeaCoefficient> BitAnd<SeaMotorPGA3D<T>> for &SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;
	fn bitand(self, rhs: SeaMotorPGA3D<T>) -> Self::Output {
		self.regressive_product(rhs)
	}
}
impl<T:SeaCoefficient> BitAnd<&SeaMotorPGA3D<T>> for &SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;
	fn bitand(self, rhs: &SeaMotorPGA3D<T>) -> Self::Output {
		self.regressive_product(rhs)
	}
}

// Flector
impl<T:SeaCoefficient> RegressiveProduct<SeaFlectorPGA3D<T>> for SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;

	fn regressive_product(&self, rhs: SeaFlectorPGA3D<T>) -> Self::Output {
		SeaVersorPGA3D::<T>((**self).regressive_product(rhs))
	}
}
impl<T:SeaCoefficient> RegressiveProduct<&SeaFlectorPGA3D<T>> for SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;

	fn regressive_product(&self, rhs: &SeaFlectorPGA3D<T>) -> Self::Output {
		SeaVersorPGA3D::<T>((**self).regressive_product(rhs))
	}
}
impl<T:SeaCoefficient> BitAnd<SeaFlectorPGA3D<T>> for SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;
	fn bitand(self, rhs: SeaFlectorPGA3D<T>) -> Self::Output {
		self.regressive_product(rhs)
	}
}
impl<T:SeaCoefficient> BitAnd<&SeaFlectorPGA3D<T>> for SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;
	fn bitand(self, rhs: &SeaFlectorPGA3D<T>) -> Self::Output {
		self.regressive_product(rhs)
	}
}
impl<T:SeaCoefficient> BitAnd<SeaFlectorPGA3D<T>> for &SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;
	fn bitand(self, rhs: SeaFlectorPGA3D<T>) -> Self::Output {
		self.regressive_product(rhs)
	}
}
impl<T:SeaCoefficient> BitAnd<&SeaFlectorPGA3D<T>> for &SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;
	fn bitand(self, rhs: &SeaFlectorPGA3D<T>) -> Self::Output {
		self.regressive_product(rhs)
	}
}


