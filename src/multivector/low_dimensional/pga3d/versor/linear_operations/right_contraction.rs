use super::*;

impl<T:SeaCoefficient> RightContraction<Self> for SeaVersorPGA3D<T> {
	type Output = Self;

	fn right_contraction(&self, rhs: Self) -> Self::Output {
		Self((**self).right_contraction(&*rhs))
	}
}
impl<T:SeaCoefficient> RightContraction<&Self> for SeaVersorPGA3D<T> {
	type Output = Self;

	fn right_contraction(&self, rhs: &Self) -> Self::Output {
		Self((**self).right_contraction(&**rhs))
	}
}
impl<T:SeaCoefficient> Shr<Self> for SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;
	fn shr(self, rhs: Self) -> Self::Output {
		self.right_contraction(rhs)
	}
}
impl<T:SeaCoefficient> Shr<&Self> for SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;
	fn shr(self, rhs: &Self) -> Self::Output {
		self.right_contraction(rhs)
	}
}
impl<T:SeaCoefficient> Shr<SeaVersorPGA3D<T>> for &SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;
	fn shr(self, rhs: SeaVersorPGA3D<T>) -> Self::Output {
		self.right_contraction(rhs)
	}
}
impl<T:SeaCoefficient> Shr<Self> for &SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;
	fn shr(self, rhs: Self) -> Self::Output {
		self.right_contraction(rhs)
	}
}

// Multivector
impl<T:SeaCoefficient> RightContraction<SeaMultivectorPGA3D<T>> for SeaVersorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;

	fn right_contraction(&self, rhs: SeaMultivectorPGA3D<T>) -> Self::Output {
		(**self).right_contraction(rhs)
	}
}
impl<T:SeaCoefficient> RightContraction<&SeaMultivectorPGA3D<T>> for SeaVersorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;

	fn right_contraction(&self, rhs: &SeaMultivectorPGA3D<T>) -> Self::Output {
		(**self).right_contraction(rhs)
	}
}
impl<T:SeaCoefficient> Shr<SeaMultivectorPGA3D<T>> for SeaVersorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn shr(self, rhs: SeaMultivectorPGA3D<T>) -> Self::Output {
		self.right_contraction(rhs)
	}
}
impl<T:SeaCoefficient> Shr<&SeaMultivectorPGA3D<T>> for SeaVersorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn shr(self, rhs: &SeaMultivectorPGA3D<T>) -> Self::Output {
		self.right_contraction(rhs)
	}
}
impl<T:SeaCoefficient> Shr<SeaMultivectorPGA3D<T>> for &SeaVersorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn shr(self, rhs: SeaMultivectorPGA3D<T>) -> Self::Output {
		self.right_contraction(rhs)
	}
}
impl<T:SeaCoefficient> Shr<&SeaMultivectorPGA3D<T>> for &SeaVersorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn shr(self, rhs: &SeaMultivectorPGA3D<T>) -> Self::Output {
		self.right_contraction(rhs)
	}
}

// Motor
impl<T:SeaCoefficient> RightContraction<SeaMotorPGA3D<T>> for SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;

	fn right_contraction(&self, rhs: SeaMotorPGA3D<T>) -> Self::Output {
		SeaVersorPGA3D::<T>((**self).right_contraction(rhs))
	}
}
impl<T:SeaCoefficient> RightContraction<&SeaMotorPGA3D<T>> for SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;

	fn right_contraction(&self, rhs: &SeaMotorPGA3D<T>) -> Self::Output {
		SeaVersorPGA3D::<T>((**self).right_contraction(rhs))
	}
}
impl<T:SeaCoefficient> Shr<SeaMotorPGA3D<T>> for SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;
	fn shr(self, rhs: SeaMotorPGA3D<T>) -> Self::Output {
		self.right_contraction(rhs)
	}
}
impl<T:SeaCoefficient> Shr<&SeaMotorPGA3D<T>> for SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;
	fn shr(self, rhs: &SeaMotorPGA3D<T>) -> Self::Output {
		self.right_contraction(rhs)
	}
}
impl<T:SeaCoefficient> Shr<SeaMotorPGA3D<T>> for &SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;
	fn shr(self, rhs: SeaMotorPGA3D<T>) -> Self::Output {
		self.right_contraction(rhs)
	}
}
impl<T:SeaCoefficient> Shr<&SeaMotorPGA3D<T>> for &SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;
	fn shr(self, rhs: &SeaMotorPGA3D<T>) -> Self::Output {
		self.right_contraction(rhs)
	}
}

// Flector
impl<T:SeaCoefficient> RightContraction<SeaFlectorPGA3D<T>> for SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;

	fn right_contraction(&self, rhs: SeaFlectorPGA3D<T>) -> Self::Output {
		SeaVersorPGA3D::<T>((**self).right_contraction(rhs))
	}
}
impl<T:SeaCoefficient> RightContraction<&SeaFlectorPGA3D<T>> for SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;

	fn right_contraction(&self, rhs: &SeaFlectorPGA3D<T>) -> Self::Output {
		SeaVersorPGA3D::<T>((**self).right_contraction(rhs))
	}
}
impl<T:SeaCoefficient> Shr<SeaFlectorPGA3D<T>> for SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;
	fn shr(self, rhs: SeaFlectorPGA3D<T>) -> Self::Output {
		self.right_contraction(rhs)
	}
}
impl<T:SeaCoefficient> Shr<&SeaFlectorPGA3D<T>> for SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;
	fn shr(self, rhs: &SeaFlectorPGA3D<T>) -> Self::Output {
		self.right_contraction(rhs)
	}
}
impl<T:SeaCoefficient> Shr<SeaFlectorPGA3D<T>> for &SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;
	fn shr(self, rhs: SeaFlectorPGA3D<T>) -> Self::Output {
		self.right_contraction(rhs)
	}
}
impl<T:SeaCoefficient> Shr<&SeaFlectorPGA3D<T>> for &SeaVersorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;
	fn shr(self, rhs: &SeaFlectorPGA3D<T>) -> Self::Output {
		self.right_contraction(rhs)
	}
}

