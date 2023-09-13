use super::*;

impl<T:SeaCoefficient> InnerProduct<Self> for SeaMotorPGA3D<T> {
	type Output = Self;

	fn inner_product(&self, rhs: Self) -> Self::Output {
		Self((**self).inner_product(&*rhs))
	}
}
impl<T:SeaCoefficient> InnerProduct<&Self> for SeaMotorPGA3D<T> {
	type Output = Self;

	fn inner_product(&self, rhs: &Self) -> Self::Output {
		Self((**self).inner_product(&**rhs))
	}
}
impl<T:SeaCoefficient> BitOr<Self> for SeaMotorPGA3D<T> {
	type Output = SeaMotorPGA3D<T>;
	fn bitor(self, rhs: Self) -> Self::Output {
		self.inner_product(rhs)
	}
}
impl<T:SeaCoefficient> BitOr<&Self> for SeaMotorPGA3D<T> {
	type Output = SeaMotorPGA3D<T>;
	fn bitor(self, rhs: &Self) -> Self::Output {
		self.inner_product(rhs)
	}
}
impl<T:SeaCoefficient> BitOr<SeaMotorPGA3D<T>> for &SeaMotorPGA3D<T> {
	type Output = SeaMotorPGA3D<T>;
	fn bitor(self, rhs: SeaMotorPGA3D<T>) -> Self::Output {
		self.inner_product(rhs)
	}
}
impl<T:SeaCoefficient> BitOr<Self> for &SeaMotorPGA3D<T> {
	type Output = SeaMotorPGA3D<T>;
	fn bitor(self, rhs: Self) -> Self::Output {
		self.inner_product(rhs)
	}
}

//Multivector
impl<T:SeaCoefficient> InnerProduct<SeaMultivectorPGA3D<T>> for SeaMotorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;

	fn inner_product(&self, rhs: SeaMultivectorPGA3D<T>) -> Self::Output {
		(**self).inner_product(rhs)
	}
}
impl<T:SeaCoefficient> InnerProduct<&SeaMultivectorPGA3D<T>> for SeaMotorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;

	fn inner_product(&self, rhs: &SeaMultivectorPGA3D<T>) -> Self::Output {
		(**self).inner_product(rhs)
	}
}
impl<T:SeaCoefficient> BitOr<SeaMultivectorPGA3D<T>> for SeaMotorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn bitor(self, rhs: SeaMultivectorPGA3D<T>) -> Self::Output {
		self.inner_product(rhs)
	}
}
impl<T:SeaCoefficient> BitOr<&SeaMultivectorPGA3D<T>> for SeaMotorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn bitor(self, rhs: &SeaMultivectorPGA3D<T>) -> Self::Output {
		self.inner_product(rhs)
	}
}
impl<T:SeaCoefficient> BitOr<SeaMultivectorPGA3D<T>> for &SeaMotorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn bitor(self, rhs: SeaMultivectorPGA3D<T>) -> Self::Output {
		self.inner_product(rhs)
	}
}
impl<T:SeaCoefficient> BitOr<&SeaMultivectorPGA3D<T>> for &SeaMotorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn bitor(self, rhs: &SeaMultivectorPGA3D<T>) -> Self::Output {
		self.inner_product(rhs)
	}
}
