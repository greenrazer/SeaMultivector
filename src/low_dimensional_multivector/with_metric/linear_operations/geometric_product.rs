use super::*;

impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> GeometricProduct<Self> for SeaMultivectorLowSparseCoefficentsWithMetric<T, I> {
	type Output = Self;
	fn geometric_product(&self, rhs: Self) -> Self::Output {
		self.geometric_product(&rhs)
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> GeometricProduct<&Self> for SeaMultivectorLowSparseCoefficentsWithMetric<T, I> {
	type Output = Self;
	fn geometric_product(&self, rhs: &Self) -> Self::Output {
		assert_eq!(self.metric, rhs.metric);
		Self {
			metric: self.metric.clone(),
			multivector: self.multivector.geometric_product(&rhs.multivector, &self.metric)
		}
	}
}