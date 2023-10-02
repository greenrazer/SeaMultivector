use super::*;

impl<T:SeaCoefficient> From<SeaMotorPGA3D<T>> for SeaMultivectorPGA3D<T> {
	fn from(value: SeaMotorPGA3D<T>) -> Self {
		value.0	
	}
}
impl<T:SeaCoefficient> From<&SeaMotorPGA3D<T>> for SeaMultivectorPGA3D<T> {
	fn from(value: &SeaMotorPGA3D<T>) -> Self {
		value.0.clone()
	}
}
