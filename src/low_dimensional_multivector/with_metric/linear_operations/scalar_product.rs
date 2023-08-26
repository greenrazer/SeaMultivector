use super::*;

impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> ScalarProduct<Self> for SeaMultivectorLowSparseCoefficentsWithMetric<T, I> {
	type Output = T;
	fn scalar_product(&self, rhs: Self) -> Self::Output {
		self.scalar_product(&rhs)
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> ScalarProduct<&Self> for SeaMultivectorLowSparseCoefficentsWithMetric<T, I> {
	type Output = T;
	fn scalar_product(&self, rhs: &Self) -> Self::Output {
		assert_eq!(self.metric, rhs.metric);
		self.multivector.scalar_product(&rhs.multivector, &self.metric)
	}
}
