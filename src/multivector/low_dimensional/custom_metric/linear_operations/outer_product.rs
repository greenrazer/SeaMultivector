use super::*;

impl<T:SeaCoefficient, I:SeaBasisBladeSet> OuterProduct<Self> for SeaLowDimensionalMultivectorCustomMetric<T, I> {
	type Output = Self;

	fn outer_product(&self, rhs: Self) -> Self::Output {
		self.outer_product(&rhs)
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeSet> OuterProduct<&Self> for SeaLowDimensionalMultivectorCustomMetric<T, I> {
	type Output = Self;

	fn outer_product(&self, rhs: &Self) -> Self::Output {
		assert_eq!(self.metric.dimensions, rhs.metric.dimensions);
		Self {
			metric: self.metric.clone(),
			multivector: self.multivector.outer_product(&rhs.multivector)
		}
	}
}

#[cfg(test)]
mod tests {
	use crate::tests::catch_unwind_silent;

    use super::*;

	#[test]
	fn works() {
		let a = SeaLowDimensionalMultivectorCustomMetric::<f32, u64>::from_str("3,3> + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345 +e6 + 12e126 + 8e46").unwrap();
		let b = SeaLowDimensionalMultivectorCustomMetric::<f32, u64>::from_str("3,3>  + 3e13 + 1e4 + 2e124  +4e345 +83e1235 + 70e1245 + 4e1345 + 8e2345 + 12e123456").unwrap();

		let out = SeaLowDimensionalMultivectorCustomMetric::<f32, u64>::from_str("3,3>e4 +3e13 +e14 +2e24 +4e34 -e46 -6e123 +5e124 +29e134 +3e136 +6e234 +4e345 -15e1234 +83e1235 +70e1245 -14e1246 +8e1345 +24e1346 +16e2345 -4e3456 -372e12345 +83e12356 +70e12456 +4e13456 +8e23456 -700e123456").unwrap();
		assert_eq!(
			a.outer_product(&b),
			out
		);

		let a = SeaLowDimensionalMultivectorCustomMetric::<f32, u64>::from_str("3,3> + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345 +e6 + 12e126 + 8e46").unwrap();
		let b = SeaLowDimensionalMultivectorCustomMetric::<f32, u64>::from_str("3,3>  + 3e13 + 1e4 + 2e124  +4e345 +83e1235 + 70e1245 + 4e1345 + 8e2345 + 12e123456").unwrap();
		assert_eq!(
			&a ^&b,
			out
		);

		assert_eq!(
			a ^&b,
			out
		);

		let a = SeaLowDimensionalMultivectorCustomMetric::<f32, u64>::from_str("3,3> + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345 +e6 + 12e126 + 8e46").unwrap();
		assert_eq!(
			&a ^b,
			out
		);

		let b = SeaLowDimensionalMultivectorCustomMetric::<f32, u64>::from_str("3,3>  + 3e13 + 1e4 + 2e124  +4e345 +83e1235 + 70e1245 + 4e1345 + 8e2345 + 12e123456").unwrap();
		assert_eq!(
			a^b,
			out
		);
	}
}
