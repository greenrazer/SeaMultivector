use super::*;

impl<T:SeaCoefficient> From<SeaFlectorPGA3D<T>> for SeaMultivectorPGA3D<T> {
	fn from(value: SeaFlectorPGA3D<T>) -> Self {
		value.0	
	}
}
impl<T:SeaCoefficient> From<&SeaFlectorPGA3D<T>> for SeaMultivectorPGA3D<T> {
	fn from(value: &SeaFlectorPGA3D<T>) -> Self {
		value.0.clone()
	}
}
