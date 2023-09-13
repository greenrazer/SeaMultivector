use super::*;

impl<T:SeaCoefficient> RegressiveProduct<Self> for SeaMultivectorPGA3D<T> {
	type Output = Self;

	fn regressive_product(&self, rhs: Self) -> Self::Output {
		Self(self.0.regressive_product(&rhs.0, &METRIC))
	}
}
impl<T:SeaCoefficient> RegressiveProduct<&Self> for SeaMultivectorPGA3D<T> {
	type Output = Self;

	fn regressive_product(&self, rhs: &Self) -> Self::Output {
		Self(self.0.regressive_product(&rhs.0, &METRIC))
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
