use super::*;

impl<T:SeaCoefficient> RightContraction<Self> for SeaMotorPGA3D<T> {
	type Output = Self;

	fn right_contraction(&self, rhs: Self) -> Self::Output {
		Self((**self).right_contraction(&*rhs))
	}
}
impl<T:SeaCoefficient> RightContraction<&Self> for SeaMotorPGA3D<T> {
	type Output = Self;

	fn right_contraction(&self, rhs: &Self) -> Self::Output {
		Self((**self).right_contraction(&**rhs))
	}
}
impl<T:SeaCoefficient> Shr<Self> for SeaMotorPGA3D<T> {
	type Output = SeaMotorPGA3D<T>;
	fn shr(self, rhs: Self) -> Self::Output {
		self.right_contraction(rhs)
	}
}
impl<T:SeaCoefficient> Shr<&Self> for SeaMotorPGA3D<T> {
	type Output = SeaMotorPGA3D<T>;
	fn shr(self, rhs: &Self) -> Self::Output {
		self.right_contraction(rhs)
	}
}
impl<T:SeaCoefficient> Shr<SeaMotorPGA3D<T>> for &SeaMotorPGA3D<T> {
	type Output = SeaMotorPGA3D<T>;
	fn shr(self, rhs: SeaMotorPGA3D<T>) -> Self::Output {
		self.right_contraction(rhs)
	}
}
impl<T:SeaCoefficient> Shr<Self> for &SeaMotorPGA3D<T> {
	type Output = SeaMotorPGA3D<T>;
	fn shr(self, rhs: Self) -> Self::Output {
		self.right_contraction(rhs)
	}
}

// Multivector
impl<T:SeaCoefficient> RightContraction<SeaMultivectorPGA3D<T>> for SeaMotorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;

	fn right_contraction(&self, rhs: SeaMultivectorPGA3D<T>) -> Self::Output {
		(**self).right_contraction(rhs)
	}
}
impl<T:SeaCoefficient> RightContraction<&SeaMultivectorPGA3D<T>> for SeaMotorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;

	fn right_contraction(&self, rhs: &SeaMultivectorPGA3D<T>) -> Self::Output {
		(**self).right_contraction(rhs)
	}
}
impl<T:SeaCoefficient> Shr<SeaMultivectorPGA3D<T>> for SeaMotorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn shr(self, rhs: SeaMultivectorPGA3D<T>) -> Self::Output {
		self.right_contraction(rhs)
	}
}
impl<T:SeaCoefficient> Shr<&SeaMultivectorPGA3D<T>> for SeaMotorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn shr(self, rhs: &SeaMultivectorPGA3D<T>) -> Self::Output {
		self.right_contraction(rhs)
	}
}
impl<T:SeaCoefficient> Shr<SeaMultivectorPGA3D<T>> for &SeaMotorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn shr(self, rhs: SeaMultivectorPGA3D<T>) -> Self::Output {
		self.right_contraction(rhs)
	}
}
impl<T:SeaCoefficient> Shr<&SeaMultivectorPGA3D<T>> for &SeaMotorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn shr(self, rhs: &SeaMultivectorPGA3D<T>) -> Self::Output {
		self.right_contraction(rhs)
	}
}
