use super::*;

mod scalar_operations;
mod linear_operations;
mod non_linear_operations;
mod non_geometric_linear_operations;

mod iterators;

mod matrix;

#[derive(Debug, Clone, PartialEq)]
pub struct SeaFlectorPGA3D<T:SeaCoefficient>(pub SeaMultivectorPGA3D<T>);

impl<T: SeaCoefficient> MultivectorConstMetric<T> for SeaFlectorPGA3D<T> {
	const METRIC: SeaMetric = SeaMultivectorPGA3D::<T>::METRIC;
}

impl<T:SeaCoefficient> Multivector<T> for SeaFlectorPGA3D<T> {
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
impl<T:SeaCoefficient> MultivectorMut<T> for SeaFlectorPGA3D<T> {}

impl<T:SeaCoefficient> Versor<T> for SeaFlectorPGA3D<T> {
	fn inverse(&self) -> Self {
		Self((**self).reversion_conjugated())
	}
}

impl<T:SeaCoefficient> Flector<T> for SeaFlectorPGA3D<T> {}

impl<T: SeaCoefficient> Deref for SeaFlectorPGA3D<T> {
	type Target = SeaMultivectorPGA3D<T>;
	fn deref(&self) -> &Self::Target {
		&self.0
	}
}
impl<T: SeaCoefficient> DerefMut for SeaFlectorPGA3D<T> {
	fn deref_mut(&mut self) -> &mut Self::Target {
		&mut self.0
	}
}

impl<T:SeaCoefficient> Default for SeaFlectorPGA3D<T> {
	fn default() -> Self {
		Self(SeaMultivectorPGA3D::<T>::default())
	}
}
