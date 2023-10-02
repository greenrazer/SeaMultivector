use super::*;

impl<T:SeaCoefficient> OuterProduct<Self> for SeaMultivectorPGA3D<T> {
	type Output = Self;

	fn outer_product(&self, rhs: Self) -> Self::Output {
		Self(self.0.outer_product(&rhs.0))
	}
}
impl<T:SeaCoefficient> OuterProduct<&Self> for SeaMultivectorPGA3D<T> {
	type Output = Self;

	fn outer_product(&self, rhs: &Self) -> Self::Output {
		Self(self.0.outer_product(&rhs.0))
	}
}
impl<T:SeaCoefficient> BitXor<Self> for SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn bitxor(self, rhs: Self) -> Self::Output {
		self.outer_product(rhs)
	}
}
impl<T:SeaCoefficient> BitXor<&Self> for SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn bitxor(self, rhs: &Self) -> Self::Output {
		self.outer_product(rhs)
	}
}
impl<T:SeaCoefficient> BitXor<SeaMultivectorPGA3D<T>> for &SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn bitxor(self, rhs: SeaMultivectorPGA3D<T>) -> Self::Output {
		self.outer_product(rhs)
	}
}
impl<T:SeaCoefficient> BitXor<Self> for &SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn bitxor(self, rhs: Self) -> Self::Output {
		self.outer_product(rhs)
	}
}

// Motor
impl<T:SeaCoefficient> OuterProduct<SeaMotorPGA3D<T>> for SeaMultivectorPGA3D<T> {
	type Output = Self;

	fn outer_product(&self, rhs: SeaMotorPGA3D<T>) -> Self::Output {
		self.outer_product(&*rhs)
	}
}
impl<T:SeaCoefficient> OuterProduct<&SeaMotorPGA3D<T>> for SeaMultivectorPGA3D<T> {
	type Output = Self;

	fn outer_product(&self, rhs: &SeaMotorPGA3D<T>) -> Self::Output {
		self.outer_product(&**rhs)
	}
}
impl<T:SeaCoefficient> BitXor<SeaMotorPGA3D<T>> for SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn bitxor(self, rhs: SeaMotorPGA3D<T>) -> Self::Output {
		self.outer_product(rhs)
	}
}
impl<T:SeaCoefficient> BitXor<&SeaMotorPGA3D<T>> for SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn bitxor(self, rhs: &SeaMotorPGA3D<T>) -> Self::Output {
		self.outer_product(rhs)
	}
}
impl<T:SeaCoefficient> BitXor<SeaMotorPGA3D<T>> for &SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn bitxor(self, rhs: SeaMotorPGA3D<T>) -> Self::Output {
		self.outer_product(rhs)
	}
}
impl<T:SeaCoefficient> BitXor<&SeaMotorPGA3D<T>> for &SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn bitxor(self, rhs: &SeaMotorPGA3D<T>) -> Self::Output {
		self.outer_product(rhs)
	}
}

// Flector
impl<T:SeaCoefficient> OuterProduct<SeaFlectorPGA3D<T>> for SeaMultivectorPGA3D<T> {
	type Output = Self;

	fn outer_product(&self, rhs: SeaFlectorPGA3D<T>) -> Self::Output {
		self.outer_product(&*rhs)
	}
}
impl<T:SeaCoefficient> OuterProduct<&SeaFlectorPGA3D<T>> for SeaMultivectorPGA3D<T> {
	type Output = Self;

	fn outer_product(&self, rhs: &SeaFlectorPGA3D<T>) -> Self::Output {
		self.outer_product(&**rhs)
	}
}
impl<T:SeaCoefficient> BitXor<SeaFlectorPGA3D<T>> for SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn bitxor(self, rhs: SeaFlectorPGA3D<T>) -> Self::Output {
		self.outer_product(rhs)
	}
}
impl<T:SeaCoefficient> BitXor<&SeaFlectorPGA3D<T>> for SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn bitxor(self, rhs: &SeaFlectorPGA3D<T>) -> Self::Output {
		self.outer_product(rhs)
	}
}
impl<T:SeaCoefficient> BitXor<SeaFlectorPGA3D<T>> for &SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn bitxor(self, rhs: SeaFlectorPGA3D<T>) -> Self::Output {
		self.outer_product(rhs)
	}
}
impl<T:SeaCoefficient> BitXor<&SeaFlectorPGA3D<T>> for &SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn bitxor(self, rhs: &SeaFlectorPGA3D<T>) -> Self::Output {
		self.outer_product(rhs)
	}
}

// Versor
impl<T:SeaCoefficient> OuterProduct<SeaVersorPGA3D<T>> for SeaMultivectorPGA3D<T> {
	type Output = Self;

	fn outer_product(&self, rhs: SeaVersorPGA3D<T>) -> Self::Output {
		self.outer_product(&*rhs)
	}
}
impl<T:SeaCoefficient> OuterProduct<&SeaVersorPGA3D<T>> for SeaMultivectorPGA3D<T> {
	type Output = Self;

	fn outer_product(&self, rhs: &SeaVersorPGA3D<T>) -> Self::Output {
		self.outer_product(&**rhs)
	}
}
impl<T:SeaCoefficient> BitXor<SeaVersorPGA3D<T>> for SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn bitxor(self, rhs: SeaVersorPGA3D<T>) -> Self::Output {
		self.outer_product(rhs)
	}
}
impl<T:SeaCoefficient> BitXor<&SeaVersorPGA3D<T>> for SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn bitxor(self, rhs: &SeaVersorPGA3D<T>) -> Self::Output {
		self.outer_product(rhs)
	}
}
impl<T:SeaCoefficient> BitXor<SeaVersorPGA3D<T>> for &SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn bitxor(self, rhs: SeaVersorPGA3D<T>) -> Self::Output {
		self.outer_product(rhs)
	}
}
impl<T:SeaCoefficient> BitXor<&SeaVersorPGA3D<T>> for &SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn bitxor(self, rhs: &SeaVersorPGA3D<T>) -> Self::Output {
		self.outer_product(rhs)
	}
}
