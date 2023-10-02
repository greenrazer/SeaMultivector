use super::*;

impl<T:SeaCoefficient> LeftContraction<Self> for SeaFlectorPGA3D<T> {
	type Output = Self;

	fn left_contraction(&self, rhs: Self) -> Self::Output {
		Self((**self).left_contraction(&*rhs))
	}
}
impl<T:SeaCoefficient> LeftContraction<&Self> for SeaFlectorPGA3D<T> {
	type Output = Self;

	fn left_contraction(&self, rhs: &Self) -> Self::Output {
		Self((**self).left_contraction(&**rhs))
	}
}
impl<T:SeaCoefficient> Shl<Self> for SeaFlectorPGA3D<T> {
	type Output = SeaFlectorPGA3D<T>;
	fn shl(self, rhs: Self) -> Self::Output {
		self.left_contraction(rhs)
	}
}
impl<T:SeaCoefficient> Shl<&Self> for SeaFlectorPGA3D<T> {
	type Output = SeaFlectorPGA3D<T>;
	fn shl(self, rhs: &Self) -> Self::Output {
		self.left_contraction(rhs)
	}
}
impl<T:SeaCoefficient> Shl<SeaFlectorPGA3D<T>> for &SeaFlectorPGA3D<T> {
	type Output = SeaFlectorPGA3D<T>;
	fn shl(self, rhs: SeaFlectorPGA3D<T>) -> Self::Output {
		self.left_contraction(rhs)
	}
}
impl<T:SeaCoefficient> Shl<Self> for &SeaFlectorPGA3D<T> {
	type Output = SeaFlectorPGA3D<T>;
	fn shl(self, rhs: Self) -> Self::Output {
		self.left_contraction(rhs)
	}
}

// Multivector
impl<T:SeaCoefficient> LeftContraction<SeaMultivectorPGA3D<T>> for SeaFlectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;

	fn left_contraction(&self, rhs: SeaMultivectorPGA3D<T>) -> Self::Output {
		(**self).left_contraction(rhs)
	}
}
impl<T:SeaCoefficient> LeftContraction<&SeaMultivectorPGA3D<T>> for SeaFlectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;

	fn left_contraction(&self, rhs: &SeaMultivectorPGA3D<T>) -> Self::Output {
		(**self).left_contraction(rhs)
	}
}
impl<T:SeaCoefficient> Shl<SeaMultivectorPGA3D<T>> for SeaFlectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn shl(self, rhs: SeaMultivectorPGA3D<T>) -> Self::Output {
		self.left_contraction(rhs)
	}
}
impl<T:SeaCoefficient> Shl<&SeaMultivectorPGA3D<T>> for SeaFlectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn shl(self, rhs: &SeaMultivectorPGA3D<T>) -> Self::Output {
		self.left_contraction(rhs)
	}
}
impl<T:SeaCoefficient> Shl<SeaMultivectorPGA3D<T>> for &SeaFlectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn shl(self, rhs: SeaMultivectorPGA3D<T>) -> Self::Output {
		self.left_contraction(rhs)
	}
}
impl<T:SeaCoefficient> Shl<&SeaMultivectorPGA3D<T>> for &SeaFlectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn shl(self, rhs: &SeaMultivectorPGA3D<T>) -> Self::Output {
		self.left_contraction(rhs)
	}
}

//Motor
impl<T:SeaCoefficient> LeftContraction<SeaMotorPGA3D<T>> for SeaFlectorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;

	fn left_contraction(&self, rhs: SeaMotorPGA3D<T>) -> Self::Output {
		SeaVersorPGA3D::<T>((**self).left_contraction(rhs))
	}
}
impl<T:SeaCoefficient> LeftContraction<&SeaMotorPGA3D<T>> for SeaFlectorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;

	fn left_contraction(&self, rhs: &SeaMotorPGA3D<T>) -> Self::Output {
		SeaVersorPGA3D::<T>((**self).left_contraction(rhs))
	}
}
impl<T:SeaCoefficient> Shl<SeaMotorPGA3D<T>> for SeaFlectorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;
	fn shl(self, rhs: SeaMotorPGA3D<T>) -> Self::Output {
		self.left_contraction(rhs)
	}
}
impl<T:SeaCoefficient> Shl<&SeaMotorPGA3D<T>> for SeaFlectorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;
	fn shl(self, rhs: &SeaMotorPGA3D<T>) -> Self::Output {
		self.left_contraction(rhs)
	}
}
impl<T:SeaCoefficient> Shl<SeaMotorPGA3D<T>> for &SeaFlectorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;
	fn shl(self, rhs: SeaMotorPGA3D<T>) -> Self::Output {
		self.left_contraction(rhs)
	}
}
impl<T:SeaCoefficient> Shl<&SeaMotorPGA3D<T>> for &SeaFlectorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;
	fn shl(self, rhs: &SeaMotorPGA3D<T>) -> Self::Output {
		self.left_contraction(rhs)
	}
}

//Versor
impl<T:SeaCoefficient> LeftContraction<SeaVersorPGA3D<T>> for SeaFlectorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;

	fn left_contraction(&self, rhs: SeaVersorPGA3D<T>) -> Self::Output {
		SeaVersorPGA3D::<T>((**self).left_contraction(rhs))
	}
}
impl<T:SeaCoefficient> LeftContraction<&SeaVersorPGA3D<T>> for SeaFlectorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;

	fn left_contraction(&self, rhs: &SeaVersorPGA3D<T>) -> Self::Output {
		SeaVersorPGA3D::<T>((**self).left_contraction(rhs))
	}
}
impl<T:SeaCoefficient> Shl<SeaVersorPGA3D<T>> for SeaFlectorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;
	fn shl(self, rhs: SeaVersorPGA3D<T>) -> Self::Output {
		self.left_contraction(rhs)
	}
}
impl<T:SeaCoefficient> Shl<&SeaVersorPGA3D<T>> for SeaFlectorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;
	fn shl(self, rhs: &SeaVersorPGA3D<T>) -> Self::Output {
		self.left_contraction(rhs)
	}
}
impl<T:SeaCoefficient> Shl<SeaVersorPGA3D<T>> for &SeaFlectorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;
	fn shl(self, rhs: SeaVersorPGA3D<T>) -> Self::Output {
		self.left_contraction(rhs)
	}
}
impl<T:SeaCoefficient> Shl<&SeaVersorPGA3D<T>> for &SeaFlectorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;
	fn shl(self, rhs: &SeaVersorPGA3D<T>) -> Self::Output {
		self.left_contraction(rhs)
	}
}
