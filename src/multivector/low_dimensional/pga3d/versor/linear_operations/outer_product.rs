use super::*;

impl<T:SeaCoefficient> OuterProduct<Self> for SeaVersorPGA3D<T> {
	type Output = Self;

	fn outer_product(&self, rhs: Self) -> Self::Output {
		Self((**self).outer_product(&*rhs))
	}
}
impl<T:SeaCoefficient> OuterProduct<&Self> for SeaVersorPGA3D<T> {
	type Output = Self;

	fn outer_product(&self, rhs: &Self) -> Self::Output {
		Self((**self).outer_product(&**rhs))
	}
}
impl<T:SeaCoefficient> BitXor<Self> for SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;
	fn bitxor(self, rhs: Self) -> Self::Output {
		self.outer_product(rhs)
	}
}
impl<T:SeaCoefficient> BitXor<&Self> for SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;
	fn bitxor(self, rhs: &Self) -> Self::Output {
		self.outer_product(rhs)
	}
}
impl<T:SeaCoefficient> BitXor<SeaVersorPGA3D<T>> for &SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;
	fn bitxor(self, rhs: SeaVersorPGA3D<T>) -> Self::Output {
		self.outer_product(rhs)
	}
}
impl<T:SeaCoefficient> BitXor<Self> for &SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;
	fn bitxor(self, rhs: Self) -> Self::Output {
		self.outer_product(rhs)
	}
}

// Multivector
impl<T:SeaCoefficient> OuterProduct<SeaMultivectorPGA3D<T>> for SeaVersorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;

	fn outer_product(&self, rhs: SeaMultivectorPGA3D<T>) -> Self::Output {
		(**self).outer_product(rhs)
	}
}
impl<T:SeaCoefficient> OuterProduct<&SeaMultivectorPGA3D<T>> for SeaVersorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;

	fn outer_product(&self, rhs: &SeaMultivectorPGA3D<T>) -> Self::Output {
		(**self).outer_product(rhs)
	}
}
impl<T:SeaCoefficient> BitXor<SeaMultivectorPGA3D<T>> for SeaVersorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn bitxor(self, rhs: SeaMultivectorPGA3D<T>) -> Self::Output {
		self.outer_product(rhs)
	}
}
impl<T:SeaCoefficient> BitXor<&SeaMultivectorPGA3D<T>> for SeaVersorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn bitxor(self, rhs: &SeaMultivectorPGA3D<T>) -> Self::Output {
		self.outer_product(rhs)
	}
}
impl<T:SeaCoefficient> BitXor<SeaMultivectorPGA3D<T>> for &SeaVersorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn bitxor(self, rhs: SeaMultivectorPGA3D<T>) -> Self::Output {
		self.outer_product(rhs)
	}
}
impl<T:SeaCoefficient> BitXor<&SeaMultivectorPGA3D<T>> for &SeaVersorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn bitxor(self, rhs: &SeaMultivectorPGA3D<T>) -> Self::Output {
		self.outer_product(rhs)
	}
}

// Motor
impl<T:SeaCoefficient> OuterProduct<SeaMotorPGA3D<T>> for SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;

	fn outer_product(&self, rhs: SeaMotorPGA3D<T>) -> Self::Output {
		SeaVersorPGA3D::<T>((**self).outer_product(rhs))
	}
}
impl<T:SeaCoefficient> OuterProduct<&SeaMotorPGA3D<T>> for SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;

	fn outer_product(&self, rhs: &SeaMotorPGA3D<T>) -> Self::Output {
		SeaVersorPGA3D::<T>((**self).outer_product(rhs))
	}
}
impl<T:SeaCoefficient> BitXor<SeaMotorPGA3D<T>> for SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;
	fn bitxor(self, rhs: SeaMotorPGA3D<T>) -> Self::Output {
		self.outer_product(rhs)
	}
}
impl<T:SeaCoefficient> BitXor<&SeaMotorPGA3D<T>> for SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;
	fn bitxor(self, rhs: &SeaMotorPGA3D<T>) -> Self::Output {
		self.outer_product(rhs)
	}
}
impl<T:SeaCoefficient> BitXor<SeaMotorPGA3D<T>> for &SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;
	fn bitxor(self, rhs: SeaMotorPGA3D<T>) -> Self::Output {
		self.outer_product(rhs)
	}
}
impl<T:SeaCoefficient> BitXor<&SeaMotorPGA3D<T>> for &SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;
	fn bitxor(self, rhs: &SeaMotorPGA3D<T>) -> Self::Output {
		self.outer_product(rhs)
	}
}

// Flector
impl<T:SeaCoefficient> OuterProduct<SeaFlectorPGA3D<T>> for SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;

	fn outer_product(&self, rhs: SeaFlectorPGA3D<T>) -> Self::Output {
		SeaVersorPGA3D::<T>((**self).outer_product(rhs))
	}
}
impl<T:SeaCoefficient> OuterProduct<&SeaFlectorPGA3D<T>> for SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;

	fn outer_product(&self, rhs: &SeaFlectorPGA3D<T>) -> Self::Output {
		SeaVersorPGA3D::<T>((**self).outer_product(rhs))
	}
}
impl<T:SeaCoefficient> BitXor<SeaFlectorPGA3D<T>> for SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;
	fn bitxor(self, rhs: SeaFlectorPGA3D<T>) -> Self::Output {
		self.outer_product(rhs)
	}
}
impl<T:SeaCoefficient> BitXor<&SeaFlectorPGA3D<T>> for SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;
	fn bitxor(self, rhs: &SeaFlectorPGA3D<T>) -> Self::Output {
		self.outer_product(rhs)
	}
}
impl<T:SeaCoefficient> BitXor<SeaFlectorPGA3D<T>> for &SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;
	fn bitxor(self, rhs: SeaFlectorPGA3D<T>) -> Self::Output {
		self.outer_product(rhs)
	}
}
impl<T:SeaCoefficient> BitXor<&SeaFlectorPGA3D<T>> for &SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;
	fn bitxor(self, rhs: &SeaFlectorPGA3D<T>) -> Self::Output {
		self.outer_product(rhs)
	}
}

