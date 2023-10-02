use super::*;

impl<T: SeaCoefficient> Into<SeaCoreMultivector4D<T>> for SeaMultivectorPGA3D<T> {
	fn into(self) -> SeaCoreMultivector4D<T> {
		self.0
	}
}
impl<'a, T: SeaCoefficient> Into<&'a SeaCoreMultivector4D<T>> for &'a SeaMultivectorPGA3D<T> {
	fn into(self) -> &'a SeaCoreMultivector4D<T> {
		&self.0
	}
}
