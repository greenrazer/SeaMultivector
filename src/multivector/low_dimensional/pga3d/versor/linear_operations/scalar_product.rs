use super::*;

impl<T:SeaCoefficient> ScalarProduct<Self> for SeaVersorPGA3D<T> {
	type Output = T;

	fn scalar_product(&self, rhs: Self) -> Self::Output {
		(**self).scalar_product(&*rhs)
	}
}
impl<T:SeaCoefficient> ScalarProduct<&Self> for SeaVersorPGA3D<T> {
	type Output = T;

	fn scalar_product(&self, rhs: &Self) -> Self::Output {
		(**self).scalar_product(&**rhs)
	}
}

// Multivector
impl<T:SeaCoefficient> ScalarProduct<SeaMultivectorPGA3D<T>> for SeaVersorPGA3D<T> {
	type Output = T;

	fn scalar_product(&self, rhs: SeaMultivectorPGA3D<T>) -> Self::Output {
		(**self).scalar_product(rhs)
	}
}
impl<T:SeaCoefficient> ScalarProduct<&SeaMultivectorPGA3D<T>> for SeaVersorPGA3D<T> {
	type Output = T;

	fn scalar_product(&self, rhs: &SeaMultivectorPGA3D<T>) -> Self::Output {
		(**self).scalar_product(rhs)
	}
}

// Motor
impl<T:SeaCoefficient> ScalarProduct<SeaMotorPGA3D<T>> for SeaVersorPGA3D<T> {
	type Output = T;

	fn scalar_product(&self, rhs: SeaMotorPGA3D<T>) -> Self::Output {
		(**self).scalar_product(rhs)
	}
}
impl<T:SeaCoefficient> ScalarProduct<&SeaMotorPGA3D<T>> for SeaVersorPGA3D<T> {
	type Output = T;

	fn scalar_product(&self, rhs: &SeaMotorPGA3D<T>) -> Self::Output {
		(**self).scalar_product(rhs)
	}
}

// Flector
impl<T:SeaCoefficient> ScalarProduct<SeaFlectorPGA3D<T>> for SeaVersorPGA3D<T> {
	type Output = T;

	fn scalar_product(&self, rhs: SeaFlectorPGA3D<T>) -> Self::Output {
		(**self).scalar_product(rhs)
	}
}
impl<T:SeaCoefficient> ScalarProduct<&SeaFlectorPGA3D<T>> for SeaVersorPGA3D<T> {
	type Output = T;

	fn scalar_product(&self, rhs: &SeaFlectorPGA3D<T>) -> Self::Output {
		(**self).scalar_product(rhs)
	}
}

