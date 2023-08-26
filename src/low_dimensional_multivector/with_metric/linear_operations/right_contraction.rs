use super::*;

impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> RightContraction<Self> for SeaMultivectorLowSparseCoefficentsWithMetric<T, I> {
	type Output = Self;
	fn right_contraction(&self, rhs: Self) -> Self::Output {
		self.right_contraction(&rhs)
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> RightContraction<&Self> for SeaMultivectorLowSparseCoefficentsWithMetric<T, I> {
	type Output = Self;
	fn right_contraction(&self, rhs: &Self) -> Self::Output {
		assert_eq!(self.metric, rhs.metric);
		Self {
			metric: self.metric.clone(),
			multivector: self.multivector.right_contraction(&rhs.multivector, &self.metric)
		}
	}
}