use super::*;

impl<T:SeaCoefficient> OuterProduct<Self> for SeaMotorPGA3D<T> {
	type Output = Self;

	fn outer_product(&self, rhs: Self) -> Self::Output {
		Self((**self).outer_product(&*rhs))
	}
}
impl<T:SeaCoefficient> OuterProduct<&Self> for SeaMotorPGA3D<T> {
	type Output = Self;

	fn outer_product(&self, rhs: &Self) -> Self::Output {
		Self((**self).outer_product(&**rhs))
	}
}
impl<T:SeaCoefficient> BitXor<Self> for SeaMotorPGA3D<T> {
	type Output = SeaMotorPGA3D<T>;
	fn bitxor(self, rhs: Self) -> Self::Output {
		self.outer_product(rhs)
	}
}
impl<T:SeaCoefficient> BitXor<&Self> for SeaMotorPGA3D<T> {
	type Output = SeaMotorPGA3D<T>;
	fn bitxor(self, rhs: &Self) -> Self::Output {
		self.outer_product(rhs)
	}
}
impl<T:SeaCoefficient> BitXor<SeaMotorPGA3D<T>> for &SeaMotorPGA3D<T> {
	type Output = SeaMotorPGA3D<T>;
	fn bitxor(self, rhs: SeaMotorPGA3D<T>) -> Self::Output {
		self.outer_product(rhs)
	}
}
impl<T:SeaCoefficient> BitXor<Self> for &SeaMotorPGA3D<T> {
	type Output = SeaMotorPGA3D<T>;
	fn bitxor(self, rhs: Self) -> Self::Output {
		self.outer_product(rhs)
	}
}

// Multivector
impl<T:SeaCoefficient> OuterProduct<SeaMultivectorPGA3D<T>> for SeaMotorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;

	fn outer_product(&self, rhs: SeaMultivectorPGA3D<T>) -> Self::Output {
		(**self).outer_product(rhs)
	}
}
impl<T:SeaCoefficient> OuterProduct<&SeaMultivectorPGA3D<T>> for SeaMotorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;

	fn outer_product(&self, rhs: &SeaMultivectorPGA3D<T>) -> Self::Output {
		(**self).outer_product(rhs)
	}
}
impl<T:SeaCoefficient> BitXor<SeaMultivectorPGA3D<T>> for SeaMotorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn bitxor(self, rhs: SeaMultivectorPGA3D<T>) -> Self::Output {
		self.outer_product(rhs)
	}
}
impl<T:SeaCoefficient> BitXor<&SeaMultivectorPGA3D<T>> for SeaMotorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn bitxor(self, rhs: &SeaMultivectorPGA3D<T>) -> Self::Output {
		self.outer_product(rhs)
	}
}
impl<T:SeaCoefficient> BitXor<SeaMultivectorPGA3D<T>> for &SeaMotorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn bitxor(self, rhs: SeaMultivectorPGA3D<T>) -> Self::Output {
		self.outer_product(rhs)
	}
}
impl<T:SeaCoefficient> BitXor<&SeaMultivectorPGA3D<T>> for &SeaMotorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn bitxor(self, rhs: &SeaMultivectorPGA3D<T>) -> Self::Output {
		self.outer_product(rhs)
	}
}
