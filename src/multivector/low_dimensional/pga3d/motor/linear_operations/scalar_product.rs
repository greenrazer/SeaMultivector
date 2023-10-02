use super::*;

impl<T:SeaCoefficient> ScalarProduct<Self> for SeaMotorPGA3D<T> {
	type Output = T;

	fn scalar_product(&self, rhs: Self) -> Self::Output {
		(**self).scalar_product(&*rhs)
	}
}
impl<T:SeaCoefficient> ScalarProduct<&Self> for SeaMotorPGA3D<T> {
	type Output = T;

	fn scalar_product(&self, rhs: &Self) -> Self::Output {
		(**self).scalar_product(&**rhs)
	}
}

// Multivector
impl<T:SeaCoefficient> ScalarProduct<SeaMultivectorPGA3D<T>> for SeaMotorPGA3D<T> {
	type Output = T;

	fn scalar_product(&self, rhs: SeaMultivectorPGA3D<T>) -> Self::Output {
		(**self).scalar_product(rhs)
	}
}
impl<T:SeaCoefficient> ScalarProduct<&SeaMultivectorPGA3D<T>> for SeaMotorPGA3D<T> {
	type Output = T;

	fn scalar_product(&self, rhs: &SeaMultivectorPGA3D<T>) -> Self::Output {
		(**self).scalar_product(rhs)
	}
}

// Flector
impl<T:SeaCoefficient> ScalarProduct<SeaFlectorPGA3D<T>> for SeaMotorPGA3D<T> {
	type Output = T;

	fn scalar_product(&self, rhs: SeaFlectorPGA3D<T>) -> Self::Output {
		(**self).scalar_product(rhs)
	}
}
impl<T:SeaCoefficient> ScalarProduct<&SeaFlectorPGA3D<T>> for SeaMotorPGA3D<T> {
	type Output = T;

	fn scalar_product(&self, rhs: &SeaFlectorPGA3D<T>) -> Self::Output {
		(**self).scalar_product(rhs)
	}
}

// Versor
impl<T:SeaCoefficient> ScalarProduct<SeaVersorPGA3D<T>> for SeaMotorPGA3D<T> {
	type Output = T;

	fn scalar_product(&self, rhs: SeaVersorPGA3D<T>) -> Self::Output {
		(**self).scalar_product(rhs)
	}
}
impl<T:SeaCoefficient> ScalarProduct<&SeaVersorPGA3D<T>> for SeaMotorPGA3D<T> {
	type Output = T;

	fn scalar_product(&self, rhs: &SeaVersorPGA3D<T>) -> Self::Output {
		(**self).scalar_product(rhs)
	}
}



