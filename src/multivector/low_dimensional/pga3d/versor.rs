use super::*;

mod scalar_operations;
mod linear_operations;
mod non_linear_operations;
mod non_geometric_linear_operations;

mod iterators;

#[derive(Debug, Clone, PartialEq)]
pub struct SeaVersorPGA3D<T:SeaCoefficient>(pub SeaMultivectorPGA3D<T>);

impl<T: SeaCoefficient> MultivectorConstMetric<T> for SeaVersorPGA3D<T> {
	const METRIC: SeaMetric = SeaMultivectorPGA3D::<T>::METRIC;
}

impl<T:SeaCoefficient> Multivector<T> for SeaVersorPGA3D<T> {
	type MultivectorIndexType = SeaBasisBladeIndexLow;
	type GeneralMultivectorOutput = SeaMultivectorPGA3D<T>;

	fn basis_blade_coefficient(&self, index: &Self::MultivectorIndexType) -> T {
		(**self).basis_blade_coefficient(index)
	}

	fn scalar(&self) -> T {
		(**self).scalar()
	}

	fn pseudoscalar(&self) -> Self::GeneralMultivectorOutput {
		(**self).pseudoscalar()
	}

	fn grade(&self, grade: &SeaGrade) -> Self::GeneralMultivectorOutput {
		(**self).grade(grade)
	}
}
impl<T:SeaCoefficient> MultivectorMut<T> for SeaVersorPGA3D<T> {}

impl<T:SeaCoefficient> Versor<T> for SeaVersorPGA3D<T> {
	fn inverse(&self) -> Self {
		Self((**self).reversion_conjugated())
	}
}

impl<T: SeaCoefficient> Deref for SeaVersorPGA3D<T> {
	type Target = SeaMultivectorPGA3D<T>;
	fn deref(&self) -> &Self::Target {
		&self.0
	}
}
impl<T: SeaCoefficient> DerefMut for SeaVersorPGA3D<T> {
	fn deref_mut(&mut self) -> &mut Self::Target {
		&mut self.0
	}
}

impl<T:SeaCoefficient> Default for SeaVersorPGA3D<T> {
	fn default() -> Self {
		Self(SeaMultivectorPGA3D::<T>::default())
	}
}
