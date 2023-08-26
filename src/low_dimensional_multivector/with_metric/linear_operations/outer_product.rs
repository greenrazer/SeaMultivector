use super::*;

impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> OuterProduct<Self> for SeaMultivectorLowSparseCoefficentsWithMetric<T, I> {
	type Output = Self;
	fn outer_product(&self, rhs: Self) -> Self::Output {
		self.outer_product(&rhs)
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> OuterProduct<&Self> for SeaMultivectorLowSparseCoefficentsWithMetric<T, I> {
	type Output = Self;
	fn outer_product(&self, rhs: &Self) -> Self::Output {
		assert_eq!(self.metric.dimensions, rhs.metric.dimensions);
		Self {
			metric: self.metric.clone(),
			multivector: self.multivector.outer_product(&rhs.multivector)
		}
	}
}