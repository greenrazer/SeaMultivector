use super::*;

impl<T:SeaCoefficient, I:SeaBasisBladeBits> SubAssign<T> for SeaMultivectorLowSparseCoefficentsWithMetric<T, I> {
    fn sub_assign(&mut self, rhs: T){
		self.multivector -= rhs
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeBits> Sub<T> for SeaMultivectorLowSparseCoefficentsWithMetric<T, I> {
    type Output = SeaMultivectorLowSparseCoefficentsWithMetric<T, I>;

    fn sub(self, rhs: T) -> Self::Output {
		&self - rhs
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeBits> Sub<T> for &SeaMultivectorLowSparseCoefficentsWithMetric<T, I> {
    type Output = SeaMultivectorLowSparseCoefficentsWithMetric<T, I>;

    fn sub(self, rhs: T) -> Self::Output {
		SeaMultivectorLowSparseCoefficentsWithMetric {
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
		let mut a = SeaMultivectorLowSparseCoefficentsWithMetric::<f32, u64>::from_str("3,3> + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345 +e6 + 12e126 + 8e46").unwrap();

		let out = SeaMultivectorLowSparseCoefficentsWithMetric::<f32, u64>::from_str("3,3> -1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345 +e6 + 12e126 + 8e46").unwrap();
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