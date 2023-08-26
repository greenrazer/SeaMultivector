use super::*;

impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> SubAssign<T> for SeaMultivectorLowSparseCoefficentsWithMetric<T, I> {
    fn sub_assign(&mut self, rhs: T){
		self.multivector -= rhs
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> Sub<T> for SeaMultivectorLowSparseCoefficentsWithMetric<T, I> {
    type Output = SeaMultivectorLowSparseCoefficentsWithMetric<T, I>;

    fn sub(self, rhs: T) -> Self::Output {
		&self - rhs
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> Sub<T> for &SeaMultivectorLowSparseCoefficentsWithMetric<T, I> {
    type Output = SeaMultivectorLowSparseCoefficentsWithMetric<T, I>;

    fn sub(self, rhs: T) -> Self::Output {
		SeaMultivectorLowSparseCoefficentsWithMetric {
			metric: self.metric.clone(),
			multivector: &self.multivector - rhs
		}
	}
}
