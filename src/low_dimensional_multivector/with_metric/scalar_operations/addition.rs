use super::*;

impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> AddAssign<T> for SeaMultivectorLowSparseCoefficentsWithMetric<T, I> {
    fn add_assign(&mut self, rhs: T){
		self.multivector += rhs
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> Add<T> for SeaMultivectorLowSparseCoefficentsWithMetric<T, I> {
    type Output = SeaMultivectorLowSparseCoefficentsWithMetric<T, I>;

    fn add(self, rhs: T) -> Self::Output {
		&self + rhs
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> Add<T> for &SeaMultivectorLowSparseCoefficentsWithMetric<T, I> {
    type Output = SeaMultivectorLowSparseCoefficentsWithMetric<T, I>;

    fn add(self, rhs: T) -> Self::Output {
		SeaMultivectorLowSparseCoefficentsWithMetric {
			metric: self.metric.clone(),
			multivector: &self.multivector + rhs
		}
	}
}
