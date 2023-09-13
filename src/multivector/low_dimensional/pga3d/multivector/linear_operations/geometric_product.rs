use super::*;

// Self
impl<T:SeaCoefficient> GeometricProduct<Self> for SeaMultivectorPGA3D<T> {
	type Output = Self;

	fn geometric_product(&self, rhs: Self) -> Self::Output {
		Self(self.0.geometric_product(&rhs.0, &METRIC))
	}
}
impl<T:SeaCoefficient> GeometricProduct<&Self> for SeaMultivectorPGA3D<T> {
	type Output = Self;

	fn geometric_product(&self, rhs: &Self) -> Self::Output {
		Self(self.0.geometric_product(&rhs.0, &METRIC))
	}
}
impl<T:SeaCoefficient> Mul<Self> for SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn mul(self, rhs: Self) -> Self::Output {
		self.geometric_product(rhs)
	}
}
impl<T:SeaCoefficient> Mul<&Self> for SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn mul(self, rhs: &Self) -> Self::Output {
		self.geometric_product(rhs)
	}
}
impl<T:SeaCoefficient> Mul<SeaMultivectorPGA3D<T>> for &SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn mul(self, rhs: SeaMultivectorPGA3D<T>) -> Self::Output {
		self.geometric_product(rhs)
	}
}
impl<T:SeaCoefficient> Mul<Self> for &SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn mul(self, rhs: Self) -> Self::Output {
		self.geometric_product(rhs)
	}
}

// Motor
impl<T:SeaCoefficient> GeometricProduct<SeaMotorPGA3D<T>> for SeaMultivectorPGA3D<T> {
	type Output = Self;

	fn geometric_product(&self, rhs: SeaMotorPGA3D<T>) -> Self::Output {
		self.geometric_product(&*rhs)
	}
}
impl<T:SeaCoefficient> GeometricProduct<&SeaMotorPGA3D<T>> for SeaMultivectorPGA3D<T> {
	type Output = Self;

	fn geometric_product(&self, rhs: &SeaMotorPGA3D<T>) -> Self::Output {
		self.geometric_product(&**rhs)
	}
}
impl<T:SeaCoefficient> Mul<SeaMotorPGA3D<T>> for SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn mul(self, rhs: SeaMotorPGA3D<T>) -> Self::Output {
		self.geometric_product(rhs)
	}
}
impl<T:SeaCoefficient> Mul<&SeaMotorPGA3D<T>> for SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn mul(self, rhs: &SeaMotorPGA3D<T>) -> Self::Output {
		self.geometric_product(rhs)
	}
}
impl<T:SeaCoefficient> Mul<SeaMotorPGA3D<T>> for &SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn mul(self, rhs: SeaMotorPGA3D<T>) -> Self::Output {
		self.geometric_product(rhs)
	}
}
impl<T:SeaCoefficient> Mul<&SeaMotorPGA3D<T>> for &SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn mul(self, rhs: &SeaMotorPGA3D<T>) -> Self::Output {
		self.geometric_product(rhs)
	}
}

