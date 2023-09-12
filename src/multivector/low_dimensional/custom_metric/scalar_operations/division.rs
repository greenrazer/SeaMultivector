use super::*;

impl<T:SeaCoefficient, I:SeaBasisBladeSet> DivAssign<T> for SeaMultivectorLowSparseCoefficentsWithMetric<T, I> {
    fn div_assign(&mut self, rhs: T) {
		self.multivector /= rhs
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeSet> Div<T> for SeaMultivectorLowSparseCoefficentsWithMetric<T, I> {
    type Output = SeaMultivectorLowSparseCoefficentsWithMetric<T, I>;

    fn div(self, rhs: T) -> Self::Output {
		&self / rhs
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeSet> Div<T> for &SeaMultivectorLowSparseCoefficentsWithMetric<T, I> {
    type Output = SeaMultivectorLowSparseCoefficentsWithMetric<T, I>;

    fn div(self, rhs: T) -> Self::Output {
		SeaMultivectorLowSparseCoefficentsWithMetric {
			metric: self.metric.clone(),
			multivector: &self.multivector / rhs
		}
	}
}

#[cfg(test)]
mod tests {
	use crate::tests::catch_unwind_silent;

	use super::*;

	#[test]
	fn div() {
		let mut a = SeaMultivectorLowSparseCoefficentsWithMetric::<f32, u64>::from_str("3,3> + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345 +e6 + 12e126 + 8e46").unwrap();

		let out = SeaMultivectorLowSparseCoefficentsWithMetric::<f32, u64>::from_str("3,3> +0.5 + 0.5e1 + 1e2 + 1.5e12 + 2e3 + 2.5e13 + 3e23 + 3.5e123 + 4e4 + 4.5e14 + 5e24 + 5.5e124 + 6e345 + 6.5e1345 + 7e2345 + 7.5e12345 +0.5e6 + 6e126 + 4e46").unwrap();
		assert_eq!(
			&a / 2.0,
			out
		);
		a /= 2.0;
		assert_eq!(
			a,
			out
		);
	}

}
