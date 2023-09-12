use super::*;

impl<T:SeaCoefficient, I:SeaBasisBladeSet> SubAssign<T> for SeaSparseMultivectorLowDimensionsWithMetric<T, I> {
    fn sub_assign(&mut self, rhs: T){
		self.multivector -= rhs
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeSet> Sub<T> for SeaSparseMultivectorLowDimensionsWithMetric<T, I> {
    type Output = SeaSparseMultivectorLowDimensionsWithMetric<T, I>;

    fn sub(self, rhs: T) -> Self::Output {
		&self - rhs
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeSet> Sub<T> for &SeaSparseMultivectorLowDimensionsWithMetric<T, I> {
    type Output = SeaSparseMultivectorLowDimensionsWithMetric<T, I>;

    fn sub(self, rhs: T) -> Self::Output {
		SeaSparseMultivectorLowDimensionsWithMetric {
			metric: self.metric.clone(),
			multivector: &self.multivector - rhs
		}
	}
}

#[cfg(test)]
mod tests {
	use crate::tests::catch_unwind_silent;

	use super::*;
	
	#[test]
	fn sub() {
		let mut a = SeaSparseMultivectorLowDimensionsWithMetric::<f32, u64>::from_str("3,3> + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345 +e6 + 12e126 + 8e46").unwrap();

		let out = SeaSparseMultivectorLowDimensionsWithMetric::<f32, u64>::from_str("3,3> -1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345 +e6 + 12e126 + 8e46").unwrap();
		assert_eq!(
			&a - 2.0,
			out
		);
		a -= 2.0;
		assert_eq!(
			a,
			out
		);
	}

}