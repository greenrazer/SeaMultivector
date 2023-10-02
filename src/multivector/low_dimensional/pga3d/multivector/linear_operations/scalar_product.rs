use super::*;

impl<T:SeaCoefficient> ScalarProduct<Self> for SeaMultivectorPGA3D<T> {
	type Output = T;

	fn scalar_product(&self, rhs: Self) -> Self::Output {
		self.0.scalar_product(&rhs.0, &Self::METRIC)
	}
}
impl<T:SeaCoefficient> ScalarProduct<&Self> for SeaMultivectorPGA3D<T> {
	type Output = T;

	fn scalar_product(&self, rhs: &Self) -> Self::Output {
		self.0.scalar_product(&rhs.0, &Self::METRIC)
	}
}

// Motor
impl<T:SeaCoefficient> ScalarProduct<SeaMotorPGA3D<T>> for SeaMultivectorPGA3D<T> {
	type Output = T;

	fn scalar_product(&self, rhs: SeaMotorPGA3D<T>) -> Self::Output {
		self.scalar_product(&*rhs)
	}
}
impl<T:SeaCoefficient> ScalarProduct<&SeaMotorPGA3D<T>> for SeaMultivectorPGA3D<T> {
	type Output = T;

	fn scalar_product(&self, rhs: &SeaMotorPGA3D<T>) -> Self::Output {
		self.scalar_product(&**rhs)
	}
}

// Flector
impl<T:SeaCoefficient> ScalarProduct<SeaFlectorPGA3D<T>> for SeaMultivectorPGA3D<T> {
	type Output = T;

	fn scalar_product(&self, rhs: SeaFlectorPGA3D<T>) -> Self::Output {
		self.scalar_product(&*rhs)
	}
}
impl<T:SeaCoefficient> ScalarProduct<&SeaFlectorPGA3D<T>> for SeaMultivectorPGA3D<T> {
	type Output = T;

	fn scalar_product(&self, rhs: &SeaFlectorPGA3D<T>) -> Self::Output {
		self.scalar_product(&**rhs)
	}
}

// Versor
impl<T:SeaCoefficient> ScalarProduct<SeaVersorPGA3D<T>> for SeaMultivectorPGA3D<T> {
	type Output = T;

	fn scalar_product(&self, rhs: SeaVersorPGA3D<T>) -> Self::Output {
		self.scalar_product(&*rhs)
	}
}
impl<T:SeaCoefficient> ScalarProduct<&SeaVersorPGA3D<T>> for SeaMultivectorPGA3D<T> {
	type Output = T;

	fn scalar_product(&self, rhs: &SeaVersorPGA3D<T>) -> Self::Output {
		self.scalar_product(&**rhs)
	}
}
