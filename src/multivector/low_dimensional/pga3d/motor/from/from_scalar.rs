use super::*;

impl<T:SeaCoefficient> From<T> for SeaMotorPGA3D<T> {
	fn from(value: T) -> Self {
		SeaMultivectorPGA3D::<T>::new_scalar(value)
	}
}