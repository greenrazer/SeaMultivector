use super::*;

impl<T:SeaCoefficient> LeftContraction<Self> for SeaVersorPGA3D<T> {
	type Output = Self;

	fn left_contraction(&self, rhs: Self) -> Self::Output {
		Self((**self).left_contraction(&*rhs))
	}
}
impl<T:SeaCoefficient> LeftContraction<&Self> for SeaVersorPGA3D<T> {
	type Output = Self;

	fn left_contraction(&self, rhs: &Self) -> Self::Output {
		Self((**self).left_contraction(&**rhs))
	}
}
impl<T:SeaCoefficient> Shl<Self> for SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;
	fn shl(self, rhs: Self) -> Self::Output {
		self.left_contraction(rhs)
	}
}
impl<T:SeaCoefficient> Shl<&Self> for SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;
	fn shl(self, rhs: &Self) -> Self::Output {
		self.left_contraction(rhs)
	}
}
impl<T:SeaCoefficient> Shl<SeaVersorPGA3D<T>> for &SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;
	fn shl(self, rhs: SeaVersorPGA3D<T>) -> Self::Output {
		self.left_contraction(rhs)
	}
}
impl<T:SeaCoefficient> Shl<Self> for &SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;
	fn shl(self, rhs: Self) -> Self::Output {
		self.left_contraction(rhs)
	}
}

// Multivector
impl<T:SeaCoefficient> LeftContraction<SeaMultivectorPGA3D<T>> for SeaVersorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;

	fn left_contraction(&self, rhs: SeaMultivectorPGA3D<T>) -> Self::Output {
		(**self).left_contraction(rhs)
	}
}
impl<T:SeaCoefficient> LeftContraction<&SeaMultivectorPGA3D<T>> for SeaVersorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;

	fn left_contraction(&self, rhs: &SeaMultivectorPGA3D<T>) -> Self::Output {
		(**self).left_contraction(rhs)
	}
}
impl<T:SeaCoefficient> Shl<SeaMultivectorPGA3D<T>> for SeaVersorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn shl(self, rhs: SeaMultivectorPGA3D<T>) -> Self::Output {
		self.left_contraction(rhs)
	}
}
impl<T:SeaCoefficient> Shl<&SeaMultivectorPGA3D<T>> for SeaVersorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn shl(self, rhs: &SeaMultivectorPGA3D<T>) -> Self::Output {
		self.left_contraction(rhs)
	}
}
impl<T:SeaCoefficient> Shl<SeaMultivectorPGA3D<T>> for &SeaVersorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn shl(self, rhs: SeaMultivectorPGA3D<T>) -> Self::Output {
		self.left_contraction(rhs)
	}
}
impl<T:SeaCoefficient> Shl<&SeaMultivectorPGA3D<T>> for &SeaVersorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn shl(self, rhs: &SeaMultivectorPGA3D<T>) -> Self::Output {
		self.left_contraction(rhs)
	}
}

// Motor
impl<T:SeaCoefficient> LeftContraction<SeaMotorPGA3D<T>> for SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;

	fn left_contraction(&self, rhs: SeaMotorPGA3D<T>) -> Self::Output {
		SeaVersorPGA3D::<T>((**self).left_contraction(rhs))
	}
}
impl<T:SeaCoefficient> LeftContraction<&SeaMotorPGA3D<T>> for SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;

	fn left_contraction(&self, rhs: &SeaMotorPGA3D<T>) -> Self::Output {
		SeaVersorPGA3D::<T>((**self).left_contraction(rhs))
	}
}
impl<T:SeaCoefficient> Shl<SeaMotorPGA3D<T>> for SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;
	fn shl(self, rhs: SeaMotorPGA3D<T>) -> Self::Output {
		self.left_contraction(rhs)
	}
}
impl<T:SeaCoefficient> Shl<&SeaMotorPGA3D<T>> for SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;
	fn shl(self, rhs: &SeaMotorPGA3D<T>) -> Self::Output {
		self.left_contraction(rhs)
	}
}
impl<T:SeaCoefficient> Shl<SeaMotorPGA3D<T>> for &SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;
	fn shl(self, rhs: SeaMotorPGA3D<T>) -> Self::Output {
		self.left_contraction(rhs)
	}
}
impl<T:SeaCoefficient> Shl<&SeaMotorPGA3D<T>> for &SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;
	fn shl(self, rhs: &SeaMotorPGA3D<T>) -> Self::Output {
		self.left_contraction(rhs)
	}
}

// Flector
impl<T:SeaCoefficient> LeftContraction<SeaFlectorPGA3D<T>> for SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;

	fn left_contraction(&self, rhs: SeaFlectorPGA3D<T>) -> Self::Output {
		SeaVersorPGA3D::<T>((**self).left_contraction(rhs))
	}
}
impl<T:SeaCoefficient> LeftContraction<&SeaFlectorPGA3D<T>> for SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;

	fn left_contraction(&self, rhs: &SeaFlectorPGA3D<T>) -> Self::Output {
		SeaVersorPGA3D::<T>((**self).left_contraction(rhs))
	}
}
impl<T:SeaCoefficient> Shl<SeaFlectorPGA3D<T>> for SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;
	fn shl(self, rhs: SeaFlectorPGA3D<T>) -> Self::Output {
		self.left_contraction(rhs)
	}
}
impl<T:SeaCoefficient> Shl<&SeaFlectorPGA3D<T>> for SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;
	fn shl(self, rhs: &SeaFlectorPGA3D<T>) -> Self::Output {
		self.left_contraction(rhs)
	}
}
impl<T:SeaCoefficient> Shl<SeaFlectorPGA3D<T>> for &SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;
	fn shl(self, rhs: SeaFlectorPGA3D<T>) -> Self::Output {
		self.left_contraction(rhs)
	}
}
impl<T:SeaCoefficient> Shl<&SeaFlectorPGA3D<T>> for &SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;
	fn shl(self, rhs: &SeaFlectorPGA3D<T>) -> Self::Output {
		self.left_contraction(rhs)
	}
}
