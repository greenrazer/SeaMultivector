use super::*;

impl<T:SeaCoefficient> From<T> for SeaMultivectorPGA3D<T> {
	fn from(value: T) -> Self {
		Self::new_scalar(value).into()	
	}
}

