use super::*;

impl<T:SeaCoefficient> InnerProduct<Self> for SeaMultivectorPGA3D<T> {
	type Output = Self;

	fn inner_product(&self, rhs: Self) -> Self::Output {
		Self(self.0.inner_product(&rhs.0, &Self::METRIC))
	}
}
impl<T:SeaCoefficient> InnerProduct<&Self> for SeaMultivectorPGA3D<T> {
	type Output = Self;

	fn inner_product(&self, rhs: &Self) -> Self::Output {
		Self(self.0.inner_product(&rhs.0, &Self::METRIC))
	}
}
impl<T:SeaCoefficient> BitOr<Self> for SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn bitor(self, rhs: Self) -> Self::Output {
		self.inner_product(rhs)
	}
}
impl<T:SeaCoefficient> BitOr<&Self> for SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn bitor(self, rhs: &Self) -> Self::Output {
		self.inner_product(rhs)
	}
}
impl<T:SeaCoefficient> BitOr<SeaMultivectorPGA3D<T>> for &SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn bitor(self, rhs: SeaMultivectorPGA3D<T>) -> Self::Output {
		self.inner_product(rhs)
	}
}
impl<T:SeaCoefficient> BitOr<Self> for &SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn bitor(self, rhs: Self) -> Self::Output {
		self.inner_product(rhs)
	}
}

// Motor
impl<T:SeaCoefficient> InnerProduct<SeaMotorPGA3D<T>> for SeaMultivectorPGA3D<T> {
	type Output = Self;

	fn inner_product(&self, rhs: SeaMotorPGA3D<T>) -> Self::Output {
		self.inner_product(&*rhs)
	}
}
impl<T:SeaCoefficient> InnerProduct<&SeaMotorPGA3D<T>> for SeaMultivectorPGA3D<T> {
	type Output = Self;

	fn inner_product(&self, rhs: &SeaMotorPGA3D<T>) -> Self::Output {
		self.inner_product(&**rhs)
	}
}
impl<T:SeaCoefficient> BitOr<SeaMotorPGA3D<T>> for SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn bitor(self, rhs: SeaMotorPGA3D<T>) -> Self::Output {
		self.inner_product(rhs)
	}
}
impl<T:SeaCoefficient> BitOr<&SeaMotorPGA3D<T>> for SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn bitor(self, rhs: &SeaMotorPGA3D<T>) -> Self::Output {
		self.inner_product(rhs)
	}
}
impl<T:SeaCoefficient> BitOr<SeaMotorPGA3D<T>> for &SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn bitor(self, rhs: SeaMotorPGA3D<T>) -> Self::Output {
		self.inner_product(rhs)
	}
}
impl<T:SeaCoefficient> BitOr<&SeaMotorPGA3D<T>> for &SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn bitor(self, rhs: &SeaMotorPGA3D<T>) -> Self::Output {
		self.inner_product(rhs)
	}
}

// Flector
impl<T:SeaCoefficient> InnerProduct<SeaFlectorPGA3D<T>> for SeaMultivectorPGA3D<T> {
	type Output = Self;

	fn inner_product(&self, rhs: SeaFlectorPGA3D<T>) -> Self::Output {
		self.inner_product(&*rhs)
	}
}
impl<T:SeaCoefficient> InnerProduct<&SeaFlectorPGA3D<T>> for SeaMultivectorPGA3D<T> {
	type Output = Self;

	fn inner_product(&self, rhs: &SeaFlectorPGA3D<T>) -> Self::Output {
		self.inner_product(&**rhs)
	}
}
impl<T:SeaCoefficient> BitOr<SeaFlectorPGA3D<T>> for SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn bitor(self, rhs: SeaFlectorPGA3D<T>) -> Self::Output {
		self.inner_product(rhs)
	}
}
impl<T:SeaCoefficient> BitOr<&SeaFlectorPGA3D<T>> for SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn bitor(self, rhs: &SeaFlectorPGA3D<T>) -> Self::Output {
		self.inner_product(rhs)
	}
}
impl<T:SeaCoefficient> BitOr<SeaFlectorPGA3D<T>> for &SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn bitor(self, rhs: SeaFlectorPGA3D<T>) -> Self::Output {
		self.inner_product(rhs)
	}
}
impl<T:SeaCoefficient> BitOr<&SeaFlectorPGA3D<T>> for &SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn bitor(self, rhs: &SeaFlectorPGA3D<T>) -> Self::Output {
		self.inner_product(rhs)
	}
}

// Versor
impl<T:SeaCoefficient> InnerProduct<SeaVersorPGA3D<T>> for SeaMultivectorPGA3D<T> {
	type Output = Self;

	fn inner_product(&self, rhs: SeaVersorPGA3D<T>) -> Self::Output {
		self.inner_product(&*rhs)
	}
}
impl<T:SeaCoefficient> InnerProduct<&SeaVersorPGA3D<T>> for SeaMultivectorPGA3D<T> {
	type Output = Self;

	fn inner_product(&self, rhs: &SeaVersorPGA3D<T>) -> Self::Output {
		self.inner_product(&**rhs)
	}
}
impl<T:SeaCoefficient> BitOr<SeaVersorPGA3D<T>> for SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn bitor(self, rhs: SeaVersorPGA3D<T>) -> Self::Output {
		self.inner_product(rhs)
	}
}
impl<T:SeaCoefficient> BitOr<&SeaVersorPGA3D<T>> for SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn bitor(self, rhs: &SeaVersorPGA3D<T>) -> Self::Output {
		self.inner_product(rhs)
	}
}
impl<T:SeaCoefficient> BitOr<SeaVersorPGA3D<T>> for &SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn bitor(self, rhs: SeaVersorPGA3D<T>) -> Self::Output {
		self.inner_product(rhs)
	}
}
impl<T:SeaCoefficient> BitOr<&SeaVersorPGA3D<T>> for &SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn bitor(self, rhs: &SeaVersorPGA3D<T>) -> Self::Output {
		self.inner_product(rhs)
	}
}

