use super::*;

impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> AddAssign<Self> for SeaMultivectorLowSparseCoefficentsWithMetric<T, I> {
	fn add_assign(&mut self, rhs: Self) {
		*self += &rhs
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> AddAssign<&Self> for SeaMultivectorLowSparseCoefficentsWithMetric<T, I> {
	fn add_assign(&mut self, rhs: &Self) {
		assert_eq!(self.metric, rhs.metric);
		self.multivector += &rhs.multivector
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> Add<Self> for SeaMultivectorLowSparseCoefficentsWithMetric<T, I> {
    type Output = SeaMultivectorLowSparseCoefficentsWithMetric<T, I>;

    fn add(self, rhs: Self) -> Self::Output {
		&self + &rhs
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> Add<&Self> for SeaMultivectorLowSparseCoefficentsWithMetric<T, I> {
    type Output = SeaMultivectorLowSparseCoefficentsWithMetric<T, I>;

    fn add(self, rhs: &Self) -> Self::Output {
		&self + rhs
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> Add<SeaMultivectorLowSparseCoefficentsWithMetric<T, I>> for &SeaMultivectorLowSparseCoefficentsWithMetric<T, I> {
    type Output = SeaMultivectorLowSparseCoefficentsWithMetric<T, I>;

    fn add(self, rhs: SeaMultivectorLowSparseCoefficentsWithMetric<T, I>) -> Self::Output {
		self + &rhs
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> Add<Self> for &SeaMultivectorLowSparseCoefficentsWithMetric<T, I> {
    type Output = SeaMultivectorLowSparseCoefficentsWithMetric<T, I>;

    fn add(self, rhs: Self) -> Self::Output {
		assert_eq!(self.metric, rhs.metric);
		SeaMultivectorLowSparseCoefficentsWithMetric {
			metric: self.metric.clone(),
			multivector: &self.multivector + &rhs.multivector
		}
	}
}