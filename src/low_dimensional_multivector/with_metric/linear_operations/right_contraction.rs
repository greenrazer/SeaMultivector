use super::*;

impl<T:SeaCoefficient, I:SeaBasisBladeSet> RightContraction<Self> for SeaMultivectorLowSparseCoefficentsWithMetric<T, I> {
	fn right_contraction(&self, rhs: Self) -> Self {
		self.right_contraction(&rhs)
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeSet> RightContraction<&Self> for SeaMultivectorLowSparseCoefficentsWithMetric<T, I> {
	fn right_contraction(&self, rhs: &Self) -> Self {
		assert_eq!(self.metric, rhs.metric);
		Self {
			metric: self.metric.clone(),
			multivector: self.multivector.right_contraction(&rhs.multivector, &self.metric)
		}
	}
}

#[cfg(test)]
mod tests {
	use crate::tests::catch_unwind_silent;

	use super::*;

	#[test]
	fn works() {
		let a = SeaMultivectorLowSparseCoefficentsWithMetric::<f32, u64>::from_str("3,3> + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345 +e6 + 12e126 + 8e46").unwrap();
		let b = SeaMultivectorLowSparseCoefficentsWithMetric::<f32, u64>::from_str("3,3>  + 3e13 + 1e4 + 2e124  +4e345 +83e1235 + 70e1245 + 4e1345 + 8e2345 + 12e123456").unwrap();

		let out = SeaMultivectorLowSparseCoefficentsWithMetric::<f32, u64>::from_str("3,3>115+59e1-105e2+1050e3+1245e4+8e6-71e12-18e35-39e45+13e135+14e235+45e245+15e1235").unwrap();
		assert_eq!(
			a.right_contraction(&b),
			out
		);

		let a = SeaMultivectorLowSparseCoefficentsWithMetric::<f32, u64>::from_str("3,3> + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345 +e6 + 12e126 + 8e46").unwrap();
		let b = SeaMultivectorLowSparseCoefficentsWithMetric::<f32, u64>::from_str("3,3>  + 3e13 + 1e4 + 2e124  +4e345 +83e1235 + 70e1245 + 4e1345 + 8e2345 + 12e123456").unwrap();
		assert_eq!(
			&a >> &b,
			out
		);
		assert_eq!(
			&a >> b,
			out
		);

		let b = SeaMultivectorLowSparseCoefficentsWithMetric::<f32, u64>::from_str("3,3>  + 3e13 + 1e4 + 2e124  +4e345 +83e1235 + 70e1245 + 4e1345 + 8e2345 + 12e123456").unwrap();
		assert_eq!(
			a >> &b,
			out
		);

		let a = SeaMultivectorLowSparseCoefficentsWithMetric::<f32, u64>::from_str("3,3> + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345 +e6 + 12e126 + 8e46").unwrap();
		assert_eq!(
			a >> b,
			out
		);
	}

	#[test]
	fn panics_on_mismatched_metrics() {
		let a = SeaMultivectorLowSparseCoefficentsWithMetric::<f32, u64>::from_str("1,2,3> + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345 +e6 + 12e126 + 8e46").unwrap();
		let b = SeaMultivectorLowSparseCoefficentsWithMetric::<f32, u64>::from_str("3,2,1>  + 3e13 + 1e4 + 2e124  +4e345 +83e1235 + 70e1245 + 4e1345 + 8e2345 + 12e123456").unwrap();

		let result = catch_unwind_silent(|| a.right_contraction(&b));
		assert!(result.is_err());
	}
}
