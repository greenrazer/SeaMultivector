use super::*;

impl<T:SeaCoefficient> LeftContraction<Self> for SeaMotorPGA3D<T> {
	type Output = Self;

	fn left_contraction(&self, rhs: Self) -> Self::Output {
		Self((**self).left_contraction(&*rhs))
	}
}
impl<T:SeaCoefficient> LeftContraction<&Self> for SeaMotorPGA3D<T> {
	type Output = Self;

	fn left_contraction(&self, rhs: &Self) -> Self::Output {
		Self((**self).left_contraction(&**rhs))
	}
}
impl<T:SeaCoefficient> Shl<Self> for SeaMotorPGA3D<T> {
	type Output = SeaMotorPGA3D<T>;
	fn shl(self, rhs: Self) -> Self::Output {
		self.left_contraction(rhs)
	}
}
impl<T:SeaCoefficient> Shl<&Self> for SeaMotorPGA3D<T> {
	type Output = SeaMotorPGA3D<T>;
	fn shl(self, rhs: &Self) -> Self::Output {
		self.left_contraction(rhs)
	}
}
impl<T:SeaCoefficient> Shl<SeaMotorPGA3D<T>> for &SeaMotorPGA3D<T> {
	type Output = SeaMotorPGA3D<T>;
	fn shl(self, rhs: SeaMotorPGA3D<T>) -> Self::Output {
		self.left_contraction(rhs)
	}
}
impl<T:SeaCoefficient> Shl<Self> for &SeaMotorPGA3D<T> {
	type Output = SeaMotorPGA3D<T>;
	fn shl(self, rhs: Self) -> Self::Output {
		self.left_contraction(rhs)
	}
}

// Multivector
impl<T:SeaCoefficient> LeftContraction<SeaMultivectorPGA3D<T>> for SeaMotorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;

	fn left_contraction(&self, rhs: SeaMultivectorPGA3D<T>) -> Self::Output {
		(**self).left_contraction(rhs)
	}
}
impl<T:SeaCoefficient> LeftContraction<&SeaMultivectorPGA3D<T>> for SeaMotorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;

	fn left_contraction(&self, rhs: &SeaMultivectorPGA3D<T>) -> Self::Output {
		(**self).left_contraction(rhs)
	}
}
impl<T:SeaCoefficient> Shl<SeaMultivectorPGA3D<T>> for SeaMotorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn shl(self, rhs: SeaMultivectorPGA3D<T>) -> Self::Output {
		self.left_contraction(rhs)
	}
}
impl<T:SeaCoefficient> Shl<&SeaMultivectorPGA3D<T>> for SeaMotorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn shl(self, rhs: &SeaMultivectorPGA3D<T>) -> Self::Output {
		self.left_contraction(rhs)
	}
}
impl<T:SeaCoefficient> Shl<SeaMultivectorPGA3D<T>> for &SeaMotorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn shl(self, rhs: SeaMultivectorPGA3D<T>) -> Self::Output {
		self.left_contraction(rhs)
	}
}
impl<T:SeaCoefficient> Shl<&SeaMultivectorPGA3D<T>> for &SeaMotorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn shl(self, rhs: &SeaMultivectorPGA3D<T>) -> Self::Output {
		self.left_contraction(rhs)
	}
}
