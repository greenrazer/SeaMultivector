use super::*;

impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> MulAssign<T> for SeaMultivectorLowSparseCoefficentsWithMetric<T, I> {
    fn mul_assign(&mut self, rhs: T){
		self.multivector *= rhs
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> Mul<T> for SeaMultivectorLowSparseCoefficentsWithMetric<T, I> {
    type Output = SeaMultivectorLowSparseCoefficentsWithMetric<T, I>;

    fn mul(self, rhs: T) -> Self::Output {
		&self * rhs
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> Mul<T> for &SeaMultivectorLowSparseCoefficentsWithMetric<T, I> {
    type Output = SeaMultivectorLowSparseCoefficentsWithMetric<T, I>;

    fn mul(self, rhs: T) -> Self::Output {
		SeaMultivectorLowSparseCoefficentsWithMetric {
			metric: self.metric.clone(),
			multivector: &self.multivector * rhs
		}
	}
}
