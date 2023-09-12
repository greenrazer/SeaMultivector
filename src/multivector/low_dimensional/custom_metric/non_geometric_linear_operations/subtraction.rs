use super::*;

impl<T:SeaCoefficient, I:SeaBasisBladeSet> SubAssign<Self> for SeaLowDimensionalMultivectorCustomMetric<T, I> {
	fn sub_assign(&mut self, rhs: Self) {
		*self -= &rhs
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeSet> SubAssign<&Self> for SeaLowDimensionalMultivectorCustomMetric<T, I> {
	fn sub_assign(&mut self, rhs: &Self) {
		assert_eq!(self.metric, rhs.metric);
		self.multivector -= &rhs.multivector
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeSet> Sub<Self> for SeaLowDimensionalMultivectorCustomMetric<T, I> {
    type Output = SeaLowDimensionalMultivectorCustomMetric<T, I>;

    fn sub(self, rhs: Self) -> Self::Output {
		&self - &rhs
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeSet> Sub<&Self> for SeaLowDimensionalMultivectorCustomMetric<T, I> {
    type Output = SeaLowDimensionalMultivectorCustomMetric<T, I>;

    fn sub(self, rhs: &Self) -> Self::Output {
		&self - rhs
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeSet> Sub<SeaLowDimensionalMultivectorCustomMetric<T, I>> for &SeaLowDimensionalMultivectorCustomMetric<T, I> {
    type Output = SeaLowDimensionalMultivectorCustomMetric<T, I>;

    fn sub(self, rhs: SeaLowDimensionalMultivectorCustomMetric<T, I>) -> Self::Output {
		self - &rhs
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeSet> Sub<Self> for &SeaLowDimensionalMultivectorCustomMetric<T, I> {
    type Output = SeaLowDimensionalMultivectorCustomMetric<T, I>;

    fn sub(self, rhs: Self) -> Self::Output {
		assert_eq!(self.metric, rhs.metric);
		SeaLowDimensionalMultivectorCustomMetric {
			metric: self.metric.clone(),
			multivector: &self.multivector - &rhs.multivector
		}
	}
}