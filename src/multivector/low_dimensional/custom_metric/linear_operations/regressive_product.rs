use super::*;

impl<T:SeaCoefficient, I:SeaBasisBladeSet> RegressiveProduct<Self> for SeaLowDimensionalMultivectorCustomMetric<T, I> {
	type Output = Self;

	fn regressive_product(&self, rhs: Self) -> Self::Output {
		self.regressive_product(&rhs)
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeSet> RegressiveProduct<&Self> for SeaLowDimensionalMultivectorCustomMetric<T, I> {
	type Output = Self;

	fn regressive_product(&self, rhs: &Self) -> Self::Output {
		assert_eq!(self.metric, rhs.metric);
		Self {
			metric: self.metric.clone(),
			multivector: self.multivector.regressive_product(&rhs.multivector, &self.metric)
		}
	}
}

impl<T:SeaCoefficient, I:SeaBasisBladeSet> BitAnd<Self> for SeaLowDimensionalMultivectorCustomMetric<T, I> {
	type Output = SeaLowDimensionalMultivectorCustomMetric<T, I>;
	
	fn bitand(self, rhs: Self) -> Self::Output {
		self.regressive_product(rhs)
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeSet> BitAnd<&Self> for SeaLowDimensionalMultivectorCustomMetric<T, I> {
	type Output = SeaLowDimensionalMultivectorCustomMetric<T, I>;
	
	fn bitand(self, rhs: &Self) -> Self::Output {
		self.regressive_product(rhs)
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeSet> BitAnd<SeaLowDimensionalMultivectorCustomMetric<T, I>> for &SeaLowDimensionalMultivectorCustomMetric<T, I> {
	type Output = SeaLowDimensionalMultivectorCustomMetric<T, I>;
	
	fn bitand(self, rhs: SeaLowDimensionalMultivectorCustomMetric<T, I>) -> Self::Output {
		self.regressive_product(rhs)
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeSet> BitAnd<Self> for &SeaLowDimensionalMultivectorCustomMetric<T, I> {
	type Output = SeaLowDimensionalMultivectorCustomMetric<T, I>;
	
	fn bitand(self, rhs: Self) -> Self::Output {
		self.regressive_product(rhs)
	}
}

#[cfg(test)]
mod tests {
	use crate::tests::catch_unwind_silent;

    use super::*;

	#[test]
	fn works() {
		let a = SeaLowDimensionalMultivectorCustomMetric::<f32, u32>::from_str("3,2> + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345").unwrap();
		let b = SeaLowDimensionalMultivectorCustomMetric::<f32, u32>::from_str("3,2>  + 3e13 + 1e4 + 2e124  +4e345").unwrap();

		let out = SeaLowDimensionalMultivectorCustomMetric::<f32, u32>::from_str("3,2> 12 +70e3 +83e4 +45e13 +26e14 +28e24 +30e124 +60e345").unwrap();
		assert_eq!(
			a.regressive_product(&b),
			out
		);

		assert_eq!(
			&a & &b,
			out
		);

		assert_eq!(
			a & &b,
			out
		);

		let a = SeaLowDimensionalMultivectorCustomMetric::<f32, u32>::from_str("3,2> + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345").unwrap();
		assert_eq!(
			&a & b,
			out
		);

		let b = SeaLowDimensionalMultivectorCustomMetric::<f32, u32>::from_str("3,2>  + 3e13 + 1e4 + 2e124  +4e345").unwrap();
		assert_eq!(
			a & b,
			out
		);
	}
}
