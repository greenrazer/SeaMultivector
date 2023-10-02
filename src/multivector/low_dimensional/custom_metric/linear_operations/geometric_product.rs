use super::*;

impl<T:SeaCoefficient, I:SeaBasisBladeSet> GeometricProduct<Self> for SeaLowDimensionalMultivectorCustomMetric<T, I> {
	type Output = Self;

	fn geometric_product(&self, rhs: Self) -> Self::Output {
		self.geometric_product(&rhs)
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeSet> GeometricProduct<&Self> for SeaLowDimensionalMultivectorCustomMetric<T, I> {
	type Output = Self;

	fn geometric_product(&self, rhs: &Self) -> Self::Output {
		assert_eq!(self.metric, rhs.metric);
		Self {
			metric: self.metric.clone(),
			multivector: self.multivector.geometric_product(&rhs.multivector, &self.metric)
		}
	}
}

impl<T:SeaCoefficient, I:SeaBasisBladeSet> Mul<Self> for SeaLowDimensionalMultivectorCustomMetric<T, I> {
	type Output = SeaLowDimensionalMultivectorCustomMetric<T, I>;
	fn mul(self, rhs: Self) -> Self::Output {
		self.geometric_product(rhs)
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeSet> Mul<&Self> for SeaLowDimensionalMultivectorCustomMetric<T, I> {
	type Output = SeaLowDimensionalMultivectorCustomMetric<T, I>;
	fn mul(self, rhs: &Self) -> Self::Output {
		self.geometric_product(rhs)
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeSet> Mul<SeaLowDimensionalMultivectorCustomMetric<T, I>> for &SeaLowDimensionalMultivectorCustomMetric<T, I> {
	type Output = SeaLowDimensionalMultivectorCustomMetric<T, I>;
	fn mul(self, rhs: SeaLowDimensionalMultivectorCustomMetric<T, I>) -> Self::Output {
		self.geometric_product(rhs)
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeSet> Mul<Self> for &SeaLowDimensionalMultivectorCustomMetric<T, I> {
	type Output = SeaLowDimensionalMultivectorCustomMetric<T, I>;
	fn mul(self, rhs: Self) -> Self::Output {
		self.geometric_product(rhs)
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

		let out = SeaLowDimensionalMultivectorCustomMetric::<f32, u64>::from_str("3,3>115+115e1-27e2+1053e3+1240e4+189e5+188e6-21e12+983e13+1159e14+202e15+168e16-919e23-1075e24-215e25-156e26+17e34-351e35-301e45+23e46-846e123-991e124-228e125-128e126+17e134-265e135+3e136-228e145-17e234+143e235-36e236+135e245+303e345-1160e356-924e456-15e1234+230e1235+124e1245-14e1246-560e1256+442e1345+24e1346-152e1356-72e1456-393e2345+44e2356+60e2456-40e3456-360e12345+179e12356+118e12456+76e13456-28e23456-700e123456").unwrap();
		assert_eq!(
			a.geometric_product(&b),
			out
		);

		assert_eq!(
			&a * &b,
			out
		);

		assert_eq!(
			a * &b,
			out
		);

		let a = SeaLowDimensionalMultivectorCustomMetric::<f32, u64>::from_str("3,3> + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345 +e6 + 12e126 + 8e46").unwrap();
		assert_eq!(
			&a * b,
			out
		);

		let b = SeaLowDimensionalMultivectorCustomMetric::<f32, u64>::from_str("3,3>  + 3e13 + 1e4 + 2e124  +4e345 +83e1235 + 70e1245 + 4e1345 + 8e2345 + 12e123456").unwrap();
		assert_eq!(
			a * b,
			out
		);
	}

	#[test]
	fn panics_on_mismatched_metrics() {
		let a = SeaLowDimensionalMultivectorCustomMetric::<f32, u64>::from_str("1,2,3> + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345 +e6 + 12e126 + 8e46").unwrap();
		let b = SeaLowDimensionalMultivectorCustomMetric::<f32, u64>::from_str("3,2,1>  + 3e13 + 1e4 + 2e124  +4e345 +83e1235 + 70e1245 + 4e1345 + 8e2345 + 12e123456").unwrap();

		let result = catch_unwind_silent(|| a.geometric_product(&b));
		assert!(result.is_err());
	}
}
