use super::*;

impl<T:SeaCoefficient> GeometricProduct<Self> for SeaVersorPGA3D<T> {
	type Output = Self;

	fn geometric_product(&self, rhs: Self) -> Self::Output {
		Self((**self).geometric_product(&*rhs))
	}
}
impl<T:SeaCoefficient> GeometricProduct<&Self> for SeaVersorPGA3D<T> {
	type Output = Self;

	fn geometric_product(&self, rhs: &Self) -> Self::Output {
		Self((**self).geometric_product(&**rhs))
	}
}
impl<T:SeaCoefficient> Mul<Self> for SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;
	fn mul(self, rhs: Self) -> Self::Output {
		self.geometric_product(rhs)
	}
}
impl<T:SeaCoefficient> Mul<&Self> for SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;
	fn mul(self, rhs: &Self) -> Self::Output {
		self.geometric_product(rhs)
	}
}
impl<T:SeaCoefficient> Mul<SeaVersorPGA3D<T>> for &SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;
	fn mul(self, rhs: SeaVersorPGA3D<T>) -> Self::Output {
		self.geometric_product(rhs)
	}
}
impl<T:SeaCoefficient> Mul<Self> for &SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;
	fn mul(self, rhs: Self) -> Self::Output {
		self.geometric_product(rhs)
	}
}

// Multivector
impl<T:SeaCoefficient> GeometricProduct<SeaMultivectorPGA3D<T>> for SeaVersorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;

	fn geometric_product(&self, rhs: SeaMultivectorPGA3D<T>) -> Self::Output {
		(**self).geometric_product(rhs)
	}
}
impl<T:SeaCoefficient> GeometricProduct<&SeaMultivectorPGA3D<T>> for SeaVersorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;

	fn geometric_product(&self, rhs: &SeaMultivectorPGA3D<T>) -> Self::Output {
		(**self).geometric_product(rhs)
	}
}
impl<T:SeaCoefficient> Mul<SeaMultivectorPGA3D<T>> for SeaVersorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn mul(self, rhs: SeaMultivectorPGA3D<T>) -> Self::Output {
		self.geometric_product(rhs)
	}
}
impl<T:SeaCoefficient> Mul<&SeaMultivectorPGA3D<T>> for SeaVersorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn mul(self, rhs: &SeaMultivectorPGA3D<T>) -> Self::Output {
		self.geometric_product(rhs)
	}
}
impl<T:SeaCoefficient> Mul<SeaMultivectorPGA3D<T>> for &SeaVersorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn mul(self, rhs: SeaMultivectorPGA3D<T>) -> Self::Output {
		self.geometric_product(rhs)
	}
}
impl<T:SeaCoefficient> Mul<&SeaMultivectorPGA3D<T>> for &SeaVersorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn mul(self, rhs: &SeaMultivectorPGA3D<T>) -> Self::Output {
		self.geometric_product(rhs)
	}
}

// Motor
impl<T:SeaCoefficient> GeometricProduct<SeaMotorPGA3D<T>> for SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;

	fn geometric_product(&self, rhs: SeaMotorPGA3D<T>) -> Self::Output {
		SeaVersorPGA3D::<T>((**self).geometric_product(rhs))
	}
}
impl<T:SeaCoefficient> GeometricProduct<&SeaMotorPGA3D<T>> for SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;

	fn geometric_product(&self, rhs: &SeaMotorPGA3D<T>) -> Self::Output {
		SeaVersorPGA3D::<T>((**self).geometric_product(rhs))
	}
}
impl<T:SeaCoefficient> Mul<SeaMotorPGA3D<T>> for SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;
	fn mul(self, rhs: SeaMotorPGA3D<T>) -> Self::Output {
		self.geometric_product(rhs)
	}
}
impl<T:SeaCoefficient> Mul<&SeaMotorPGA3D<T>> for SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;
	fn mul(self, rhs: &SeaMotorPGA3D<T>) -> Self::Output {
		self.geometric_product(rhs)
	}
}
impl<T:SeaCoefficient> Mul<SeaMotorPGA3D<T>> for &SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;
	fn mul(self, rhs: SeaMotorPGA3D<T>) -> Self::Output {
		self.geometric_product(rhs)
	}
}
impl<T:SeaCoefficient> Mul<&SeaMotorPGA3D<T>> for &SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;
	fn mul(self, rhs: &SeaMotorPGA3D<T>) -> Self::Output {
		self.geometric_product(rhs)
	}
}

// Flector
impl<T:SeaCoefficient> GeometricProduct<SeaFlectorPGA3D<T>> for SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;

	fn geometric_product(&self, rhs: SeaFlectorPGA3D<T>) -> Self::Output {
		SeaVersorPGA3D::<T>((**self).geometric_product(rhs))
	}
}
impl<T:SeaCoefficient> GeometricProduct<&SeaFlectorPGA3D<T>> for SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;

	fn geometric_product(&self, rhs: &SeaFlectorPGA3D<T>) -> Self::Output {
		SeaVersorPGA3D::<T>((**self).geometric_product(rhs))
	}
}
impl<T:SeaCoefficient> Mul<SeaFlectorPGA3D<T>> for SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;
	fn mul(self, rhs: SeaFlectorPGA3D<T>) -> Self::Output {
		self.geometric_product(rhs)
	}
}
impl<T:SeaCoefficient> Mul<&SeaFlectorPGA3D<T>> for SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;
	fn mul(self, rhs: &SeaFlectorPGA3D<T>) -> Self::Output {
		self.geometric_product(rhs)
	}
}
impl<T:SeaCoefficient> Mul<SeaFlectorPGA3D<T>> for &SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;
	fn mul(self, rhs: SeaFlectorPGA3D<T>) -> Self::Output {
		self.geometric_product(rhs)
	}
}
impl<T:SeaCoefficient> Mul<&SeaFlectorPGA3D<T>> for &SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;
	fn mul(self, rhs: &SeaFlectorPGA3D<T>) -> Self::Output {
		self.geometric_product(rhs)
	}
}

