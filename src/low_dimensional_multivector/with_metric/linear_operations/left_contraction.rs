use super::*;

impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> LeftContraction<Self> for SeaMultivectorLowSparseCoefficentsWithMetric<T, I> {
	type Output = Self;
	fn left_contraction(&self, rhs: Self) -> Self::Output {
		self.left_contraction(&rhs)
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> LeftContraction<&Self> for SeaMultivectorLowSparseCoefficentsWithMetric<T, I> {
	type Output = Self;
	fn left_contraction(&self, rhs: &Self) -> Self::Output {
		assert_eq!(self.metric, rhs.metric);
		Self {
			metric: self.metric.clone(),
			multivector: self.multivector.left_contraction(&rhs.multivector, &self.metric)
		}
	}
}