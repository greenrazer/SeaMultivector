use super::*;

impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> RegressiveProduct<Self> for SeaMultivectorLowSparseCoefficentsWithMetric<T, I> {
	type Output = Self;
	
	fn regressive_product(&self, rhs: Self) -> Self::Output {
		self.regressive_product(&rhs)
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> RegressiveProduct<&Self> for SeaMultivectorLowSparseCoefficentsWithMetric<T, I> {
	type Output = Self;
	
	fn regressive_product(&self, rhs: &Self) -> Self::Output {
		assert_eq!(self.metric, rhs.metric);
		Self {
			metric: self.metric.clone(),
			multivector: self.multivector.regressive_product(&rhs.multivector, &self.metric)
		}
	}
}