use super::*;

impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> ScalarProduct<Self> for SeaMultivectorLowSparseCoefficentsWithMetric<T, I> {
	type Output = T;
	fn scalar_product(&self, rhs: Self) -> Self::Output {
		self.scalar_product(&rhs)
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> ScalarProduct<&Self> for SeaMultivectorLowSparseCoefficentsWithMetric<T, I> {
	type Output = T;
	fn scalar_product(&self, rhs: &Self) -> Self::Output {
		assert_eq!(self.metric, rhs.metric);
		self.multivector.scalar_product(&rhs.multivector, &self.metric)
	}
}

#[cfg(test)]
mod tests {
	use crate::tests::catch_unwind_silent;

    use super::*;


	#[test]
	fn works() {	
		let a = SeaMultivectorLowSparseCoefficentsWithMetric::<f32, u32>::from_str("3,1,1> + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345").unwrap();
		let b = SeaMultivectorLowSparseCoefficentsWithMetric::<f32, u32>::from_str("3,1,1> + 3e13 + 1e4 + 2e124  +4e345").unwrap();

		assert_eq!(
			a.scalar_product(&b),
			-1.0
		);
	}
}

