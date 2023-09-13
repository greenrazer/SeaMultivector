use super::*;

impl<T:SeaCoefficient> LeftContraction<Self> for SeaMultivectorPGA3D<T> {
	type Output = Self;

	fn left_contraction(&self, rhs: Self) -> Self::Output {
		Self(self.0.left_contraction(&rhs.0, &METRIC))
	}
}
impl<T:SeaCoefficient> LeftContraction<&Self> for SeaMultivectorPGA3D<T> {
	type Output = Self;

	fn left_contraction(&self, rhs: &Self) -> Self::Output {
		Self(self.0.left_contraction(&rhs.0, &METRIC))
	}
}
impl<T:SeaCoefficient> Shl<Self> for SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn shl(self, rhs: Self) -> Self::Output {
		self.left_contraction(rhs)
	}
}
impl<T:SeaCoefficient> Shl<&Self> for SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn shl(self, rhs: &Self) -> Self::Output {
		self.left_contraction(rhs)
	}
}
impl<T:SeaCoefficient> Shl<SeaMultivectorPGA3D<T>> for &SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn shl(self, rhs: SeaMultivectorPGA3D<T>) -> Self::Output {
		self.left_contraction(rhs)
	}
}
impl<T:SeaCoefficient> Shl<Self> for &SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn shl(self, rhs: Self) -> Self::Output {
		self.left_contraction(rhs)
	}
}

// Motor
impl<T:SeaCoefficient> LeftContraction<SeaMotorPGA3D<T>> for SeaMultivectorPGA3D<T> {
	type Output = Self;

	fn left_contraction(&self, rhs: SeaMotorPGA3D<T>) -> Self::Output {
		self.left_contraction(&*rhs)
	}
}
impl<T:SeaCoefficient> LeftContraction<&SeaMotorPGA3D<T>> for SeaMultivectorPGA3D<T> {
	type Output = Self;

	fn left_contraction(&self, rhs: &SeaMotorPGA3D<T>) -> Self::Output {
		self.left_contraction(&**rhs)
	}
}
impl<T:SeaCoefficient> Shl<SeaMotorPGA3D<T>> for SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn shl(self, rhs: SeaMotorPGA3D<T>) -> Self::Output {
		self.left_contraction(rhs)
	}
}
impl<T:SeaCoefficient> Shl<&SeaMotorPGA3D<T>> for SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn shl(self, rhs: &SeaMotorPGA3D<T>) -> Self::Output {
		self.left_contraction(rhs)
	}
}
impl<T:SeaCoefficient> Shl<SeaMotorPGA3D<T>> for &SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn shl(self, rhs: SeaMotorPGA3D<T>) -> Self::Output {
		self.left_contraction(rhs)
	}
}
impl<T:SeaCoefficient> Shl<&SeaMotorPGA3D<T>> for &SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn shl(self, rhs: &SeaMotorPGA3D<T>) -> Self::Output {
		self.left_contraction(rhs)
	}
}