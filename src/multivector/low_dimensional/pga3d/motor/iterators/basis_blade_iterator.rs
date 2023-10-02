use super::*;

impl<T:SeaCoefficient> MultivectorIter for SeaMotorPGA3D<T> {
	fn iter(&self) -> <&Self as IntoIterator>::IntoIter {
		self.into_iter()
	}
}

impl<T:SeaCoefficient> MultivectorIterMut for SeaMotorPGA3D<T> {
	fn iter_mut(&mut self) -> <&mut Self as IntoIterator>::IntoIter{
		self.into_iter()
	}
}

impl<'a, T:SeaCoefficient> IntoIterator for &'a SeaMotorPGA3D<T> {
	type IntoIter =  SeaCoreSparseMultivector4DIter<'a, T>;
	type Item = SeaBasisBladeLowRef<'a, T>;

	fn into_iter(self) -> Self::IntoIter {
		self.0.iter()
	}

}
impl<'a, T:SeaCoefficient> IntoIterator for &'a mut SeaMotorPGA3D<T> {
	type IntoIter =  SeaCoreSparseMultivector4DIterMut<'a, T>;
	type Item = SeaBasisBladeLowRefMut<'a, T>;

	fn into_iter(self) -> Self::IntoIter {
		self.0.iter_mut()
	}

}