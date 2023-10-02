use super::*;

impl<T:SeaCoefficient> From<SeaCoreMultivector4D<T>> for SeaMultivectorPGA3D<T> {
	fn from(value: SeaCoreMultivector4D<T>) -> Self {
		Self(value)
	}
}
impl<T:SeaCoefficient> From<&SeaCoreMultivector4D<T>> for SeaMultivectorPGA3D<T> {
	fn from(value: &SeaCoreMultivector4D<T>) -> Self {
		Self(value.clone())
	}
}