use super::*;

impl<T:SeaCoefficient> From<T> for SeaMotorPGA3D<T> {
	fn from(value: T) -> Self {
		SeaMotorPGA3D(SeaMultivectorPGA3D::<T>::new_lifted_scalar(value))	
	}
}