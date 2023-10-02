use super::*;

impl<T:SeaCoefficient> From<SeaVersorPGA3D<T>> for SeaMultivectorPGA3D<T> {
	fn from(value: SeaVersorPGA3D<T>) -> Self {
		value.0	
	}
}
impl<T:SeaCoefficient> From<&SeaVersorPGA3D<T>> for SeaMultivectorPGA3D<T> {
	fn from(value: &SeaVersorPGA3D<T>) -> Self {
		value.0.clone()
	}
}
