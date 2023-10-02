use super::*;

impl<T:SeaCoefficient> InnerProduct<Self> for SeaVersorPGA3D<T> {
	type Output = Self;

	fn inner_product(&self, rhs: Self) -> Self::Output {
		Self((**self).inner_product(&*rhs))
	}
}
impl<T:SeaCoefficient> InnerProduct<&Self> for SeaVersorPGA3D<T> {
	type Output = Self;

	fn inner_product(&self, rhs: &Self) -> Self::Output {
		Self((**self).inner_product(&**rhs))
	}
}
impl<T:SeaCoefficient> BitOr<Self> for SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;
	fn bitor(self, rhs: Self) -> Self::Output {
		self.inner_product(rhs)
	}
}
impl<T:SeaCoefficient> BitOr<&Self> for SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;
	fn bitor(self, rhs: &Self) -> Self::Output {
		self.inner_product(rhs)
	}
}
impl<T:SeaCoefficient> BitOr<SeaVersorPGA3D<T>> for &SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;
	fn bitor(self, rhs: SeaVersorPGA3D<T>) -> Self::Output {
		self.inner_product(rhs)
	}
}
impl<T:SeaCoefficient> BitOr<Self> for &SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;
	fn bitor(self, rhs: Self) -> Self::Output {
		self.inner_product(rhs)
	}
}

//Multivector
impl<T:SeaCoefficient> InnerProduct<SeaMultivectorPGA3D<T>> for SeaVersorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;

	fn inner_product(&self, rhs: SeaMultivectorPGA3D<T>) -> Self::Output {
		(**self).inner_product(rhs)
	}
}
impl<T:SeaCoefficient> InnerProduct<&SeaMultivectorPGA3D<T>> for SeaVersorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;

	fn inner_product(&self, rhs: &SeaMultivectorPGA3D<T>) -> Self::Output {
		(**self).inner_product(rhs)
	}
}
impl<T:SeaCoefficient> BitOr<SeaMultivectorPGA3D<T>> for SeaVersorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn bitor(self, rhs: SeaMultivectorPGA3D<T>) -> Self::Output {
		self.inner_product(rhs)
	}
}
impl<T:SeaCoefficient> BitOr<&SeaMultivectorPGA3D<T>> for SeaVersorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn bitor(self, rhs: &SeaMultivectorPGA3D<T>) -> Self::Output {
		self.inner_product(rhs)
	}
}
impl<T:SeaCoefficient> BitOr<SeaMultivectorPGA3D<T>> for &SeaVersorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn bitor(self, rhs: SeaMultivectorPGA3D<T>) -> Self::Output {
		self.inner_product(rhs)
	}
}
impl<T:SeaCoefficient> BitOr<&SeaMultivectorPGA3D<T>> for &SeaVersorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn bitor(self, rhs: &SeaMultivectorPGA3D<T>) -> Self::Output {
		self.inner_product(rhs)
	}
}

// Motor
impl<T:SeaCoefficient> InnerProduct<SeaMotorPGA3D<T>> for SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;

	fn inner_product(&self, rhs: SeaMotorPGA3D<T>) -> Self::Output {
		SeaVersorPGA3D::<T>((**self).inner_product(rhs))
	}
}
impl<T:SeaCoefficient> InnerProduct<&SeaMotorPGA3D<T>> for SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;

	fn inner_product(&self, rhs: &SeaMotorPGA3D<T>) -> Self::Output {
		SeaVersorPGA3D::<T>((**self).inner_product(rhs))
	}
}
impl<T:SeaCoefficient> BitOr<SeaMotorPGA3D<T>> for SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;
	fn bitor(self, rhs: SeaMotorPGA3D<T>) -> Self::Output {
		self.inner_product(rhs)
	}
}
impl<T:SeaCoefficient> BitOr<&SeaMotorPGA3D<T>> for SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;
	fn bitor(self, rhs: &SeaMotorPGA3D<T>) -> Self::Output {
		self.inner_product(rhs)
	}
}
impl<T:SeaCoefficient> BitOr<SeaMotorPGA3D<T>> for &SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;
	fn bitor(self, rhs: SeaMotorPGA3D<T>) -> Self::Output {
		self.inner_product(rhs)
	}
}
impl<T:SeaCoefficient> BitOr<&SeaMotorPGA3D<T>> for &SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;
	fn bitor(self, rhs: &SeaMotorPGA3D<T>) -> Self::Output {
		self.inner_product(rhs)
	}
}

// Flector
impl<T:SeaCoefficient> InnerProduct<SeaFlectorPGA3D<T>> for SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;

	fn inner_product(&self, rhs: SeaFlectorPGA3D<T>) -> Self::Output {
		SeaVersorPGA3D::<T>((**self).inner_product(rhs))
	}
}
impl<T:SeaCoefficient> InnerProduct<&SeaFlectorPGA3D<T>> for SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;

	fn inner_product(&self, rhs: &SeaFlectorPGA3D<T>) -> Self::Output {
		SeaVersorPGA3D::<T>((**self).inner_product(rhs))
	}
}
impl<T:SeaCoefficient> BitOr<SeaFlectorPGA3D<T>> for SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;
	fn bitor(self, rhs: SeaFlectorPGA3D<T>) -> Self::Output {
		self.inner_product(rhs)
	}
}
impl<T:SeaCoefficient> BitOr<&SeaFlectorPGA3D<T>> for SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;
	fn bitor(self, rhs: &SeaFlectorPGA3D<T>) -> Self::Output {
		self.inner_product(rhs)
	}
}
impl<T:SeaCoefficient> BitOr<SeaFlectorPGA3D<T>> for &SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;
	fn bitor(self, rhs: SeaFlectorPGA3D<T>) -> Self::Output {
		self.inner_product(rhs)
	}
}
impl<T:SeaCoefficient> BitOr<&SeaFlectorPGA3D<T>> for &SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;
	fn bitor(self, rhs: &SeaFlectorPGA3D<T>) -> Self::Output {
		self.inner_product(rhs)
	}
}

