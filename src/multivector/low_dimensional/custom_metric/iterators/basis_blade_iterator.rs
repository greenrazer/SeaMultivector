use super::*;

impl<T:SeaCoefficient, I:SeaBasisBladeSet> MultivectorIter for SeaLowDimensionalMultivectorCustomMetric<T, I> {
	fn iter(&self) -> <&Self as IntoIterator>::IntoIter {
		self.into_iter()
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeSet> MultivectorIterMut for SeaLowDimensionalMultivectorCustomMetric<T, I> {
	fn iter_mut(&mut self) -> <&mut Self as IntoIterator>::IntoIter {
		self.into_iter()
	}
}

impl<'a, T:SeaCoefficient, I:SeaBasisBladeSet> IntoIterator for &'a SeaLowDimensionalMultivectorCustomMetric<T, I> {
	type Item = SeaBasisBladeLowRef<'a, T>;
	type IntoIter = SparseMultivectorIter<'a, T, I>;

	fn into_iter(self) -> Self::IntoIter {
		self.multivector.iter()
	}
}
impl<'a, T:SeaCoefficient, I:SeaBasisBladeSet> IntoIterator for &'a mut SeaLowDimensionalMultivectorCustomMetric<T, I> {
	type Item = SeaBasisBladeLowRefMut<'a, T>;
	type IntoIter = SparseMultivectorIterMut<'a, T, I>;

	fn into_iter(self) -> Self::IntoIter {
		self.multivector.iter_mut()
	}
}