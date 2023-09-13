use super::*;

impl<T:SeaCoefficient> GeometricProduct<Self> for SeaMotorPGA3D<T> {
	type Output = Self;

	fn geometric_product(&self, rhs: Self) -> Self::Output {
		Self((**self).geometric_product(&*rhs))
	}
}
impl<T:SeaCoefficient> GeometricProduct<&Self> for SeaMotorPGA3D<T> {
	type Output = Self;

	fn geometric_product(&self, rhs: &Self) -> Self::Output {
		Self((**self).geometric_product(&**rhs))
	}
}
impl<T:SeaCoefficient> Mul<Self> for SeaMotorPGA3D<T> {
	type Output = SeaMotorPGA3D<T>;
	fn mul(self, rhs: Self) -> Self::Output {
		self.geometric_product(rhs)
	}
}
impl<T:SeaCoefficient> Mul<&Self> for SeaMotorPGA3D<T> {
	type Output = SeaMotorPGA3D<T>;
	fn mul(self, rhs: &Self) -> Self::Output {
		self.geometric_product(rhs)
	}
}
impl<T:SeaCoefficient> Mul<SeaMotorPGA3D<T>> for &SeaMotorPGA3D<T> {
	type Output = SeaMotorPGA3D<T>;
	fn mul(self, rhs: SeaMotorPGA3D<T>) -> Self::Output {
		self.geometric_product(rhs)
	}
}
impl<T:SeaCoefficient> Mul<Self> for &SeaMotorPGA3D<T> {
	type Output = SeaMotorPGA3D<T>;
	fn mul(self, rhs: Self) -> Self::Output {
		self.geometric_product(rhs)
	}
}

// Multivector
impl<T:SeaCoefficient> GeometricProduct<SeaMultivectorPGA3D<T>> for SeaMotorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;

	fn geometric_product(&self, rhs: SeaMultivectorPGA3D<T>) -> Self::Output {
		(**self).geometric_product(rhs)
	}
}
impl<T:SeaCoefficient> GeometricProduct<&SeaMultivectorPGA3D<T>> for SeaMotorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;

	fn geometric_product(&self, rhs: &SeaMultivectorPGA3D<T>) -> Self::Output {
		(**self).geometric_product(rhs)
	}
}
impl<T:SeaCoefficient> Mul<SeaMultivectorPGA3D<T>> for SeaMotorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn mul(self, rhs: SeaMultivectorPGA3D<T>) -> Self::Output {
		self.geometric_product(rhs)
	}
}
impl<T:SeaCoefficient> Mul<&SeaMultivectorPGA3D<T>> for SeaMotorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn mul(self, rhs: &SeaMultivectorPGA3D<T>) -> Self::Output {
		self.geometric_product(rhs)
	}
}
impl<T:SeaCoefficient> Mul<SeaMultivectorPGA3D<T>> for &SeaMotorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn mul(self, rhs: SeaMultivectorPGA3D<T>) -> Self::Output {
		self.geometric_product(rhs)
	}
}
impl<T:SeaCoefficient> Mul<&SeaMultivectorPGA3D<T>> for &SeaMotorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn mul(self, rhs: &SeaMultivectorPGA3D<T>) -> Self::Output {
		self.geometric_product(rhs)
	}
}

