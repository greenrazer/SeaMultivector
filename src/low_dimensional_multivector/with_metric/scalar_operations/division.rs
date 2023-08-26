use super::*;

impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> DivAssign<T> for SeaMultivectorLowSparseCoefficentsWithMetric<T, I> {
    fn div_assign(&mut self, rhs: T) {
		self.multivector /= rhs
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> Div<T> for SeaMultivectorLowSparseCoefficentsWithMetric<T, I> {
    type Output = SeaMultivectorLowSparseCoefficentsWithMetric<T, I>;

    fn div(self, rhs: T) -> Self::Output {
		&self / rhs
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> Div<T> for &SeaMultivectorLowSparseCoefficentsWithMetric<T, I> {
    type Output = SeaMultivectorLowSparseCoefficentsWithMetric<T, I>;

    fn div(self, rhs: T) -> Self::Output {
		SeaMultivectorLowSparseCoefficentsWithMetric {
			metric: self.metric.clone(),
			multivector: &self.multivector / rhs
		}
	}
}
