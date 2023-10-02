use super::*;

impl<T:SeaCoefficient> RightContraction<Self> for SeaMultivectorPGA3D<T> {
	type Output = Self;

	fn right_contraction(&self, rhs: Self) -> Self {
		Self(self.0.right_contraction(&rhs.0, &Self::METRIC))
	}
}
impl<T:SeaCoefficient> RightContraction<&Self> for SeaMultivectorPGA3D<T> {
	type Output = Self;

	fn right_contraction(&self, rhs: &Self) -> Self {
		Self(self.0.right_contraction(&rhs.0, &Self::METRIC))
	}
}
impl<T:SeaCoefficient> Shr<Self> for SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn shr(self, rhs: Self) -> Self::Output {
		self.right_contraction(rhs)
	}
}
impl<T:SeaCoefficient> Shr<&Self> for SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn shr(self, rhs: &Self) -> Self::Output {
		self.right_contraction(rhs)
	}
}
impl<T:SeaCoefficient> Shr<SeaMultivectorPGA3D<T>> for &SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn shr(self, rhs: SeaMultivectorPGA3D<T>) -> Self::Output {
		self.right_contraction(rhs)
	}
}
impl<T:SeaCoefficient> Shr<Self> for &SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn shr(self, rhs: Self) -> Self::Output {
		self.right_contraction(rhs)
	}
}

// Motor
impl<T:SeaCoefficient> RightContraction<SeaMotorPGA3D<T>> for SeaMultivectorPGA3D<T> {
	type Output = Self;

	fn right_contraction(&self, rhs: SeaMotorPGA3D<T>) -> Self {
		self.right_contraction(&*rhs)
	}
}
impl<T:SeaCoefficient> RightContraction<&SeaMotorPGA3D<T>> for SeaMultivectorPGA3D<T> {
	type Output = Self;
	
	fn right_contraction(&self, rhs: &SeaMotorPGA3D<T>) -> Self {
		self.right_contraction(&**rhs)
	}
}
impl<T:SeaCoefficient> Shr<SeaMotorPGA3D<T>> for SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn shr(self, rhs: SeaMotorPGA3D<T>) -> Self::Output {
		self.right_contraction(rhs)
	}
}
impl<T:SeaCoefficient> Shr<&SeaMotorPGA3D<T>> for SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn shr(self, rhs: &SeaMotorPGA3D<T>) -> Self::Output {
		self.right_contraction(rhs)
	}
}
impl<T:SeaCoefficient> Shr<SeaMotorPGA3D<T>> for &SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn shr(self, rhs: SeaMotorPGA3D<T>) -> Self::Output {
		self.right_contraction(rhs)
	}
}
impl<T:SeaCoefficient> Shr<&SeaMotorPGA3D<T>> for &SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn shr(self, rhs: &SeaMotorPGA3D<T>) -> Self::Output {
		self.right_contraction(rhs)
	}
}

// Flector
impl<T:SeaCoefficient> RightContraction<SeaFlectorPGA3D<T>> for SeaMultivectorPGA3D<T> {
	type Output = Self;

	fn right_contraction(&self, rhs: SeaFlectorPGA3D<T>) -> Self {
		self.right_contraction(&*rhs)
	}
}
impl<T:SeaCoefficient> RightContraction<&SeaFlectorPGA3D<T>> for SeaMultivectorPGA3D<T> {
	type Output = Self;
	
	fn right_contraction(&self, rhs: &SeaFlectorPGA3D<T>) -> Self {
		self.right_contraction(&**rhs)
	}
}
impl<T:SeaCoefficient> Shr<SeaFlectorPGA3D<T>> for SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn shr(self, rhs: SeaFlectorPGA3D<T>) -> Self::Output {
		self.right_contraction(rhs)
	}
}
impl<T:SeaCoefficient> Shr<&SeaFlectorPGA3D<T>> for SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn shr(self, rhs: &SeaFlectorPGA3D<T>) -> Self::Output {
		self.right_contraction(rhs)
	}
}
impl<T:SeaCoefficient> Shr<SeaFlectorPGA3D<T>> for &SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn shr(self, rhs: SeaFlectorPGA3D<T>) -> Self::Output {
		self.right_contraction(rhs)
	}
}
impl<T:SeaCoefficient> Shr<&SeaFlectorPGA3D<T>> for &SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn shr(self, rhs: &SeaFlectorPGA3D<T>) -> Self::Output {
		self.right_contraction(rhs)
	}
}

// Versor
impl<T:SeaCoefficient> RightContraction<SeaVersorPGA3D<T>> for SeaMultivectorPGA3D<T> {
	type Output = Self;

	fn right_contraction(&self, rhs: SeaVersorPGA3D<T>) -> Self {
		self.right_contraction(&*rhs)
	}
}
impl<T:SeaCoefficient> RightContraction<&SeaVersorPGA3D<T>> for SeaMultivectorPGA3D<T> {
	type Output = Self;
	
	fn right_contraction(&self, rhs: &SeaVersorPGA3D<T>) -> Self {
		self.right_contraction(&**rhs)
	}
}
impl<T:SeaCoefficient> Shr<SeaVersorPGA3D<T>> for SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn shr(self, rhs: SeaVersorPGA3D<T>) -> Self::Output {
		self.right_contraction(rhs)
	}
}
impl<T:SeaCoefficient> Shr<&SeaVersorPGA3D<T>> for SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn shr(self, rhs: &SeaVersorPGA3D<T>) -> Self::Output {
		self.right_contraction(rhs)
	}
}
impl<T:SeaCoefficient> Shr<SeaVersorPGA3D<T>> for &SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn shr(self, rhs: SeaVersorPGA3D<T>) -> Self::Output {
		self.right_contraction(rhs)
	}
}
impl<T:SeaCoefficient> Shr<&SeaVersorPGA3D<T>> for &SeaMultivectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn shr(self, rhs: &SeaVersorPGA3D<T>) -> Self::Output {
		self.right_contraction(rhs)
	}
}