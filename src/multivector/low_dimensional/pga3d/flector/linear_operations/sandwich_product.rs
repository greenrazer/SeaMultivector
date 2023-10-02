use super::*;

impl<T:SeaCoefficient> SandwichProduct<Self> for SeaFlectorPGA3D<T> {
	type Output = Self;

	fn sandwich_product(&self, rhs: Self) -> Self::Output {
		let a_mv: &SeaCoreMultivector4D<T> = (&**self).into();
		let b_mv: SeaCoreMultivector4D<T> = rhs.0.into();
		Self(
			SeaMultivectorPGA3D::from(
				a_mv.sandwich_product(&b_mv, &Self::METRIC)
			)
		)
	}
}
impl<T:SeaCoefficient> SandwichProduct<&Self> for SeaFlectorPGA3D<T> {
	type Output = Self;

	fn sandwich_product(&self, rhs: &Self) -> Self::Output {
		let a_mv: &SeaCoreMultivector4D<T> = (&**self).into();
		let b_mv: &SeaCoreMultivector4D<T> = (&**rhs).into();
		Self(
			SeaMultivectorPGA3D::from(
				a_mv.sandwich_product(&b_mv, &Self::METRIC)
			)
		)
	}
}
impl<T:SeaCoefficient> Div<Self> for SeaFlectorPGA3D<T> {
	type Output = SeaFlectorPGA3D<T>;
	fn div(self, rhs: Self) -> Self::Output {
		self.sandwich_product(rhs)
	}
}
impl<T:SeaCoefficient> Div<&Self> for SeaFlectorPGA3D<T> {
	type Output = SeaFlectorPGA3D<T>;
	fn div(self, rhs: &Self) -> Self::Output {
		self.sandwich_product(rhs)
	}
}
impl<T:SeaCoefficient> Div<SeaFlectorPGA3D<T>> for &SeaFlectorPGA3D<T> {
	type Output = SeaFlectorPGA3D<T>;
	fn div(self, rhs: SeaFlectorPGA3D<T>) -> Self::Output {
		self.sandwich_product(rhs)
	}
}
impl<T:SeaCoefficient> Div<Self> for &SeaFlectorPGA3D<T> {
	type Output = SeaFlectorPGA3D<T>;
	fn div(self, rhs: Self) -> Self::Output {
		self.sandwich_product(rhs)
	}
}

// Multivector
impl<T:SeaCoefficient> SandwichProduct<SeaMultivectorPGA3D<T>> for SeaFlectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;

	fn sandwich_product(&self, rhs: SeaMultivectorPGA3D<T>) -> Self::Output {
		let a_mv: &SeaCoreMultivector4D<T> = (&**self).into();
		let b_mv: SeaCoreMultivector4D<T> = rhs.into();
		let out_mv = a_mv.sandwich_product(&b_mv, &Self::METRIC);
		SeaMultivectorPGA3D::from(out_mv)
	}
}
impl<T:SeaCoefficient> SandwichProduct<&SeaMultivectorPGA3D<T>> for SeaFlectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;

	fn sandwich_product(&self, rhs: &SeaMultivectorPGA3D<T>) -> Self::Output {
		let a_mv: &SeaCoreMultivector4D<T> = (&**self).into();
		let b_mv: &SeaCoreMultivector4D<T> = (&*rhs).into();
		let out_mv = a_mv.sandwich_product(&b_mv, &Self::METRIC);
		SeaMultivectorPGA3D::from(out_mv)
	}
}
impl<T:SeaCoefficient> Div<SeaMultivectorPGA3D<T>> for SeaFlectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn div(self, rhs: SeaMultivectorPGA3D<T>) -> Self::Output {
		self.sandwich_product(rhs)
	}
}
impl<T:SeaCoefficient> Div<&SeaMultivectorPGA3D<T>> for SeaFlectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn div(self, rhs: &SeaMultivectorPGA3D<T>) -> Self::Output {
		self.sandwich_product(rhs)
	}
}
impl<T:SeaCoefficient> Div<SeaMultivectorPGA3D<T>> for &SeaFlectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn div(self, rhs: SeaMultivectorPGA3D<T>) -> Self::Output {
		self.sandwich_product(rhs)
	}
}
impl<T:SeaCoefficient> Div<&SeaMultivectorPGA3D<T>> for &SeaFlectorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn div(self, rhs: &SeaMultivectorPGA3D<T>) -> Self::Output {
		self.sandwich_product(rhs)
	}
}

//Motor
impl<T:SeaCoefficient> SandwichProduct<SeaMotorPGA3D<T>> for SeaFlectorPGA3D<T> {
	type Output = SeaMotorPGA3D<T>;

	fn sandwich_product(&self, rhs: SeaMotorPGA3D<T>) -> Self::Output {
		let a_mv: &SeaCoreMultivector4D<T> = (&**self).into();
		let b_mv: SeaCoreMultivector4D<T> = rhs.0.into();
		let out_mv = a_mv.sandwich_product(&b_mv, &Self::METRIC);
		SeaMotorPGA3D::<T>(SeaMultivectorPGA3D::from(out_mv))
	}
}
impl<T:SeaCoefficient> SandwichProduct<&SeaMotorPGA3D<T>> for SeaFlectorPGA3D<T> {
	type Output = SeaMotorPGA3D<T>;

	fn sandwich_product(&self, rhs: &SeaMotorPGA3D<T>) -> Self::Output {
		let a_mv: &SeaCoreMultivector4D<T> = (&**self).into();
		let b_mv: &SeaCoreMultivector4D<T> = (&**rhs).into();
		let out_mv = a_mv.sandwich_product(&b_mv, &Self::METRIC);
		SeaMotorPGA3D::<T>(SeaMultivectorPGA3D::from(out_mv))
	}
}
impl<T:SeaCoefficient> Div<SeaMotorPGA3D<T>> for SeaFlectorPGA3D<T> {
	type Output = SeaMotorPGA3D<T>;
	fn div(self, rhs: SeaMotorPGA3D<T>) -> Self::Output {
		self.sandwich_product(rhs)
	}
}
impl<T:SeaCoefficient> Div<&SeaMotorPGA3D<T>> for SeaFlectorPGA3D<T> {
	type Output = SeaMotorPGA3D<T>;
	fn div(self, rhs: &SeaMotorPGA3D<T>) -> Self::Output {
		self.sandwich_product(rhs)
	}
}
impl<T:SeaCoefficient> Div<SeaMotorPGA3D<T>> for &SeaFlectorPGA3D<T> {
	type Output = SeaMotorPGA3D<T>;
	fn div(self, rhs: SeaMotorPGA3D<T>) -> Self::Output {
		self.sandwich_product(rhs)
	}
}
impl<T:SeaCoefficient> Div<&SeaMotorPGA3D<T>> for &SeaFlectorPGA3D<T> {
	type Output = SeaMotorPGA3D<T>;
	fn div(self, rhs: &SeaMotorPGA3D<T>) -> Self::Output {
		self.sandwich_product(rhs)
	}
}

//Versor
impl<T:SeaCoefficient> SandwichProduct<SeaVersorPGA3D<T>> for SeaFlectorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;

	fn sandwich_product(&self, rhs: SeaVersorPGA3D<T>) -> Self::Output {
		let a_mv: &SeaCoreMultivector4D<T> = (&**self).into();
		let b_mv: SeaCoreMultivector4D<T> = rhs.0.into();
		let out_mv = a_mv.sandwich_product(&b_mv, &Self::METRIC);
		SeaVersorPGA3D::<T>(SeaMultivectorPGA3D::from(out_mv))
	}
}
impl<T:SeaCoefficient> SandwichProduct<&SeaVersorPGA3D<T>> for SeaFlectorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;

	fn sandwich_product(&self, rhs: &SeaVersorPGA3D<T>) -> Self::Output {
		let a_mv: &SeaCoreMultivector4D<T> = (&**self).into();
		let b_mv: &SeaCoreMultivector4D<T> = (&**rhs).into();
		let out_mv = a_mv.sandwich_product(&b_mv, &Self::METRIC);
		SeaVersorPGA3D::<T>(SeaMultivectorPGA3D::from(out_mv))
	}
}
impl<T:SeaCoefficient> Div<SeaVersorPGA3D<T>> for SeaFlectorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;
	fn div(self, rhs: SeaVersorPGA3D<T>) -> Self::Output {
		self.sandwich_product(rhs)
	}
}
impl<T:SeaCoefficient> Div<&SeaVersorPGA3D<T>> for SeaFlectorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;
	fn div(self, rhs: &SeaVersorPGA3D<T>) -> Self::Output {
		self.sandwich_product(rhs)
	}
}
impl<T:SeaCoefficient> Div<SeaVersorPGA3D<T>> for &SeaFlectorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;
	fn div(self, rhs: SeaVersorPGA3D<T>) -> Self::Output {
		self.sandwich_product(rhs)
	}
}
impl<T:SeaCoefficient> Div<&SeaVersorPGA3D<T>> for &SeaFlectorPGA3D<T> {
	type Output = SeaVersorPGA3D<T>;
	fn div(self, rhs: &SeaVersorPGA3D<T>) -> Self::Output {
		self.sandwich_product(rhs)
	}
}