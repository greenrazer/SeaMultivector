use super::*;

impl<T:SeaCoefficient> RegressiveProduct<Self> for SeaMotorPGA3D<T> {
	type Output = Self;

	fn regressive_product(&self, rhs: Self) -> Self::Output {
		Self((**self).regressive_product(&*rhs))
	}
}
impl<T:SeaCoefficient> RegressiveProduct<&Self> for SeaMotorPGA3D<T> {
	type Output = Self;

	fn regressive_product(&self, rhs: &Self) -> Self::Output {
		Self((**self).regressive_product(&**rhs))
	}
}
impl<T:SeaCoefficient> BitAnd<Self> for SeaMotorPGA3D<T> {
	type Output = SeaMotorPGA3D<T>;
	
	fn bitand(self, rhs: Self) -> Self::Output {
		self.regressive_product(rhs)
	}
}
impl<T:SeaCoefficient> BitAnd<&Self> for SeaMotorPGA3D<T> {
	type Output = SeaMotorPGA3D<T>;
	
	fn bitand(self, rhs: &Self) -> Self::Output {
		self.regressive_product(rhs)
	}
}
impl<T:SeaCoefficient> BitAnd<SeaMotorPGA3D<T>> for &SeaMotorPGA3D<T> {
	type Output = SeaMotorPGA3D<T>;
	
	fn bitand(self, rhs: SeaMotorPGA3D<T>) -> Self::Output {
		self.regressive_product(rhs)
	}
}
impl<T:SeaCoefficient> BitAnd<Self> for &SeaMotorPGA3D<T> {
	type Output = SeaMotorPGA3D<T>;
	
	fn bitand(self, rhs: Self) -> Self::Output {
		self.regressive_product(rhs)
	}
}

// Multivector
impl<T:SeaCoefficient> RegressiveProduct<SeaMultivectorPGA3D<T>> for SeaMotorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;

	fn regressive_product(&self, rhs: SeaMultivectorPGA3D<T>) -> Self::Output {
		(**self).regressive_product(rhs)
	}
}
impl<T:SeaCoefficient> RegressiveProduct<&SeaMultivectorPGA3D<T>> for SeaMotorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;

	fn regressive_product(&self, rhs: &SeaMultivectorPGA3D<T>) -> Self::Output {
		(**self).regressive_product(rhs)
	}
}
impl<T:SeaCoefficient> BitAnd<SeaMultivectorPGA3D<T>> for SeaMotorPGA3D<T> {
	type Output =SeaMultivectorPGA3D<T>;
	
	fn bitand(self, rhs: SeaMultivectorPGA3D<T>) -> Self::Output {
		self.regressive_product(rhs)
	}
}
impl<T:SeaCoefficient> BitAnd<&SeaMultivectorPGA3D<T>> for SeaMotorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	
	fn bitand(self, rhs: &SeaMultivectorPGA3D<T>) -> Self::Output {
		self.regressive_product(rhs)
	}
}
impl<T:SeaCoefficient> BitAnd<SeaMultivectorPGA3D<T>> for &SeaMotorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	
	fn bitand(self, rhs: SeaMultivectorPGA3D<T>) -> Self::Output {
		self.regressive_product(rhs)
	}
}
impl<T:SeaCoefficient> BitAnd<&SeaMultivectorPGA3D<T>> for &SeaMotorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	
	fn bitand(self, rhs: &SeaMultivectorPGA3D<T>) -> Self::Output {
		self.regressive_product(rhs)
	}
}

