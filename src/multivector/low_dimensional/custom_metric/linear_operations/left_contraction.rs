use super::*;

impl<T:SeaCoefficient, I:SeaBasisBladeSet> LeftContraction<Self> for SeaSparseMultivectorLowDimensionsWithMetric<T, I> {
	fn left_contraction(&self, rhs: Self) -> Self {
		self.left_contraction(&rhs)
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeSet> LeftContraction<&Self> for SeaSparseMultivectorLowDimensionsWithMetric<T, I> {
	fn left_contraction(&self, rhs: &Self) -> Self {
		assert_eq!(self.metric, rhs.metric);
		Self {
			metric: self.metric.clone(),
			multivector: self.multivector.left_contraction(&rhs.multivector, &self.metric)
		}
	}
}

#[cfg(test)]
mod tests {
	use crate::tests::catch_unwind_silent;

    use super::*;


	#[test]
	fn works() {
		let a = SeaSparseMultivectorLowDimensionsWithMetric::<f32, u64>::from_str("3,3> + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345 +e6 + 12e126 + 8e46").unwrap();
		let b = SeaSparseMultivectorLowDimensionsWithMetric::<f32, u64>::from_str("3,3>  + 3e13 + 1e4 + 2e124  +4e345 +83e1235 + 70e1245 + 4e1345 + 8e2345 + 12e123456").unwrap();

		let out = SeaSparseMultivectorLowDimensionsWithMetric::<f32, u64>::from_str("3,3>115+56e1+78e2+3e3-5e4+189e5+180e6-16e12+3e13-4e14+202e15+168e16+2e24-215e25-156e26-333e35-262e45+2e124-228e125-144e126-198e135-156e145+19e235+38e245-120e345-132e356-84e456+179e1235+70e1245+4e1345-120e1356-72e1456+8e2345+108e2356+60e2456-36e3456+12e12345+96e12356+48e12456-24e13456+12e23456+12e123456").unwrap();
		assert_eq!(
			a.left_contraction(&b),
			out
		);

		let a = SeaSparseMultivectorLowDimensionsWithMetric::<f32, u64>::from_str("3,3> + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345 +e6 + 12e126 + 8e46").unwrap();
		let b = SeaSparseMultivectorLowDimensionsWithMetric::<f32, u64>::from_str("3,3>  + 3e13 + 1e4 + 2e124  +4e345 +83e1235 + 70e1245 + 4e1345 + 8e2345 + 12e123456").unwrap();

		let out = SeaSparseMultivectorLowDimensionsWithMetric::<f32, u64>::from_str("3,3>115+56e1+78e2+3e3-5e4+189e5+180e6-16e12+3e13-4e14+202e15+168e16+2e24-215e25-156e26-333e35-262e45+2e124-228e125-144e126-198e135-156e145+19e235+38e245-120e345-132e356-84e456+179e1235+70e1245+4e1345-120e1356-72e1456+8e2345+108e2356+60e2456-36e3456+12e12345+96e12356+48e12456-24e13456+12e23456+12e123456").unwrap();
		assert_eq!(
			&a << &b,
			out
		);
		assert_eq!(
			a << &b,
			out
		);

		let a = SeaSparseMultivectorLowDimensionsWithMetric::<f32, u64>::from_str("3,3> + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345 +e6 + 12e126 + 8e46").unwrap();
		assert_eq!(
			&a << b,
			out
		);

		let b = SeaSparseMultivectorLowDimensionsWithMetric::<f32, u64>::from_str("3,3>  + 3e13 + 1e4 + 2e124  +4e345 +83e1235 + 70e1245 + 4e1345 + 8e2345 + 12e123456").unwrap();
		assert_eq!(
			a << b,
			out
		);
	}

	#[test]
	fn panics_on_mismatched_metrics() {
		let a = SeaSparseMultivectorLowDimensionsWithMetric::<f32, u64>::from_str("1,2,3> + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345 +e6 + 12e126 + 8e46").unwrap();
		let b = SeaSparseMultivectorLowDimensionsWithMetric::<f32, u64>::from_str("3,2,1>  + 3e13 + 1e4 + 2e124  +4e345 +83e1235 + 70e1245 + 4e1345 + 8e2345 + 12e123456").unwrap();

		let result = catch_unwind_silent(|| a.left_contraction(&b));
		assert!(result.is_err());
	}
}
