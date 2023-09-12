use super::*;

impl<T:SeaCoefficient, I:SeaBasisBladeSet> MulAssign<T> for SeaSparseMultivectorLowDimensionsWithMetric<T, I> {
    fn mul_assign(&mut self, rhs: T){
		self.multivector *= rhs
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeSet> Mul<T> for SeaSparseMultivectorLowDimensionsWithMetric<T, I> {
    type Output = SeaSparseMultivectorLowDimensionsWithMetric<T, I>;

    fn mul(self, rhs: T) -> Self::Output {
		&self * rhs
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeSet> Mul<T> for &SeaSparseMultivectorLowDimensionsWithMetric<T, I> {
    type Output = SeaSparseMultivectorLowDimensionsWithMetric<T, I>;

    fn mul(self, rhs: T) -> Self::Output {
		SeaSparseMultivectorLowDimensionsWithMetric {
			metric: self.metric.clone(),
			multivector: &self.multivector * rhs
		}
	}
}

#[cfg(test)]
mod tests {
	use crate::tests::catch_unwind_silent;

	use super::*;

	#[test]
	fn mul() {
		let mut a = SeaSparseMultivectorLowDimensionsWithMetric::<f32, u64>::from_str("3,3> + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345 +e6 + 12e126 + 8e46").unwrap();

		let out = SeaSparseMultivectorLowDimensionsWithMetric::<f32, u64>::from_str("3,3> +2 + 2e1 + 4e2 + 6e12 + 8e3 + 10e13 + 12e23 + 14e123 + 16e4 + 18e14 + 20e24 + 22e124 + 24e345 + 26e1345 + 28e2345 + 30e12345 +2e6 + 24e126 + 16e46").unwrap();
		assert_eq!(
			&a * 2.0,
			out
		);
		a *= 2.0;
		assert_eq!(
			a,
			out
		);
	}

}
