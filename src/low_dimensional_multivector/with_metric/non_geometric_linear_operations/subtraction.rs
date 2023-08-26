use super::*;

impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> SubAssign<Self> for SeaMultivectorLowSparseCoefficentsWithMetric<T, I> {
	fn sub_assign(&mut self, rhs: Self) {
		*self -= &rhs
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> SubAssign<&Self> for SeaMultivectorLowSparseCoefficentsWithMetric<T, I> {
	fn sub_assign(&mut self, rhs: &Self) {
		assert_eq!(self.metric, rhs.metric);
		self.multivector -= &rhs.multivector
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> Sub<Self> for SeaMultivectorLowSparseCoefficentsWithMetric<T, I> {
    type Output = SeaMultivectorLowSparseCoefficentsWithMetric<T, I>;

    fn sub(self, rhs: Self) -> Self::Output {
		&self - &rhs
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> Sub<&Self> for SeaMultivectorLowSparseCoefficentsWithMetric<T, I> {
    type Output = SeaMultivectorLowSparseCoefficentsWithMetric<T, I>;

    fn sub(self, rhs: &Self) -> Self::Output {
		&self - rhs
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> Sub<SeaMultivectorLowSparseCoefficentsWithMetric<T, I>> for &SeaMultivectorLowSparseCoefficentsWithMetric<T, I> {
    type Output = SeaMultivectorLowSparseCoefficentsWithMetric<T, I>;

    fn sub(self, rhs: SeaMultivectorLowSparseCoefficentsWithMetric<T, I>) -> Self::Output {
		self - &rhs
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> Sub<Self> for &SeaMultivectorLowSparseCoefficentsWithMetric<T, I> {
    type Output = SeaMultivectorLowSparseCoefficentsWithMetric<T, I>;

    fn sub(self, rhs: Self) -> Self::Output {
		assert_eq!(self.metric, rhs.metric);
		SeaMultivectorLowSparseCoefficentsWithMetric {
			metric: self.metric.clone(),
			multivector: &self.multivector - &rhs.multivector
		}
	}
}