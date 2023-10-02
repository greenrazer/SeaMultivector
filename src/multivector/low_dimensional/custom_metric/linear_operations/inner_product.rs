use super::*;

impl<T:SeaCoefficient, I:SeaBasisBladeSet> InnerProduct<Self> for SeaLowDimensionalMultivectorCustomMetric<T, I> {
	type Output = Self;

	fn inner_product(&self, rhs: Self) -> Self::Output {
		self.inner_product(&rhs)
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeSet> InnerProduct<&Self> for SeaLowDimensionalMultivectorCustomMetric<T, I> {
	type Output = Self;

	fn inner_product(&self, rhs: &Self) -> Self::Output {
		assert_eq!(self.metric, rhs.metric);
		Self {
			metric: self.metric.clone(),
			multivector: self.multivector.inner_product(&rhs.multivector, &self.metric)
		}
	}
}

impl<T:SeaCoefficient, I:SeaBasisBladeSet> BitOr<Self> for SeaLowDimensionalMultivectorCustomMetric<T, I> {
	type Output = SeaLowDimensionalMultivectorCustomMetric<T, I>;
	fn bitor(self, rhs: Self) -> Self::Output {
		self.inner_product(rhs)
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeSet> BitOr<&Self> for SeaLowDimensionalMultivectorCustomMetric<T, I> {
	type Output = SeaLowDimensionalMultivectorCustomMetric<T, I>;
	fn bitor(self, rhs: &Self) -> Self::Output {
		self.inner_product(rhs)
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeSet> BitOr<SeaLowDimensionalMultivectorCustomMetric<T, I>> for &SeaLowDimensionalMultivectorCustomMetric<T, I> {
	type Output = SeaLowDimensionalMultivectorCustomMetric<T, I>;
	fn bitor(self, rhs: SeaLowDimensionalMultivectorCustomMetric<T, I>) -> Self::Output {
		self.inner_product(rhs)
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeSet> BitOr<Self> for &SeaLowDimensionalMultivectorCustomMetric<T, I> {
	type Output = SeaLowDimensionalMultivectorCustomMetric<T, I>;
	fn bitor(self, rhs: Self) -> Self::Output {
		self.inner_product(rhs)
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

		let out = SeaLowDimensionalMultivectorCustomMetric::<f32, u64>::from_str("3,3>115+115e1 -27e2 +1053e3 +1240e4 +189e5 +188e6 -87e12 +3e13 -4e14 +202e15 +168e16 +2e24 -215e25 -156e26 -351e35 -301e45 +2e124 -228e125 -144e126 -185e135 -156e145 +33e235 +83e245 -120e345 -132e356 -84e456 +194e1235 +70e1245 +4e1345 -120e1356 -72e1456 +8e2345 +108e2356 +60e2456 -36e3456 +12e12345 +96e12356 +48e12456 -24e13456 +12e23456 +12e123456\u{200B}").unwrap();
		assert_eq!(
			a.inner_product(&b),
			out
		);

		assert_eq!(
			&a | &b ,
			out
		);

		assert_eq!(
			a | &b ,
			out
		);

		let a = SeaLowDimensionalMultivectorCustomMetric::<f32, u64>::from_str("3,3> + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345 +e6 + 12e126 + 8e46").unwrap();
		assert_eq!(
			&a | b ,
			out
		);

		let b = SeaLowDimensionalMultivectorCustomMetric::<f32, u64>::from_str("3,3>  + 3e13 + 1e4 + 2e124  +4e345 +83e1235 + 70e1245 + 4e1345 + 8e2345 + 12e123456").unwrap();
		assert_eq!(
			a | b ,
			out
		);
	}

	#[test]
	fn panics_on_mismatched_metrics() {
		let a = SeaLowDimensionalMultivectorCustomMetric::<f32, u64>::from_str("1,2,3> + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345 +e6 + 12e126 + 8e46").unwrap();
		let b = SeaLowDimensionalMultivectorCustomMetric::<f32, u64>::from_str("3,2,1>  + 3e13 + 1e4 + 2e124  +4e345 +83e1235 + 70e1245 + 4e1345 + 8e2345 + 12e123456").unwrap();

		let result = catch_unwind_silent(|| a.inner_product(&b));
		assert!(result.is_err());
	}
}
