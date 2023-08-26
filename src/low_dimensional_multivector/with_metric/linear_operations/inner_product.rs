use super::*;

impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> InnerProduct<Self> for SeaMultivectorLowSparseCoefficentsWithMetric<T, I> {
	type Output = Self;
	fn inner_product(&self, rhs: Self) -> Self::Output {
		self.inner_product(&rhs)
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> InnerProduct<&Self> for SeaMultivectorLowSparseCoefficentsWithMetric<T, I> {
	type Output = Self;
	fn inner_product(&self, rhs: &Self) -> Self::Output {
		assert_eq!(self.metric, rhs.metric);
		Self {
			metric: self.metric.clone(),
			multivector: self.multivector.inner_product(&rhs.multivector, &self.metric)
		}
	}
}