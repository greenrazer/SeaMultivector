use super::*;

impl<T:SeaCoefficient, I:SeaBasisBladeSet> SubAssign<Self> for SeaSparseMultivectorLowDimensionsWithMetric<T, I> {
	fn sub_assign(&mut self, rhs: Self) {
		*self -= &rhs
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeSet> SubAssign<&Self> for SeaSparseMultivectorLowDimensionsWithMetric<T, I> {
	fn sub_assign(&mut self, rhs: &Self) {
		assert_eq!(self.metric, rhs.metric);
		self.multivector -= &rhs.multivector
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeSet> Sub<Self> for SeaSparseMultivectorLowDimensionsWithMetric<T, I> {
    type Output = SeaSparseMultivectorLowDimensionsWithMetric<T, I>;

    fn sub(self, rhs: Self) -> Self::Output {
		&self - &rhs
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeSet> Sub<&Self> for SeaSparseMultivectorLowDimensionsWithMetric<T, I> {
    type Output = SeaSparseMultivectorLowDimensionsWithMetric<T, I>;

    fn sub(self, rhs: &Self) -> Self::Output {
		&self - rhs
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeSet> Sub<SeaSparseMultivectorLowDimensionsWithMetric<T, I>> for &SeaSparseMultivectorLowDimensionsWithMetric<T, I> {
    type Output = SeaSparseMultivectorLowDimensionsWithMetric<T, I>;

    fn sub(self, rhs: SeaSparseMultivectorLowDimensionsWithMetric<T, I>) -> Self::Output {
		self - &rhs
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeSet> Sub<Self> for &SeaSparseMultivectorLowDimensionsWithMetric<T, I> {
    type Output = SeaSparseMultivectorLowDimensionsWithMetric<T, I>;

    fn sub(self, rhs: Self) -> Self::Output {
		assert_eq!(self.metric, rhs.metric);
		SeaSparseMultivectorLowDimensionsWithMetric {
			metric: self.metric.clone(),
			multivector: &self.multivector - &rhs.multivector
		}
	}
}