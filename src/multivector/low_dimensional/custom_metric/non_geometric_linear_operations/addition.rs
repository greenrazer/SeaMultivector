use super::*;

impl<T:SeaCoefficient, I:SeaBasisBladeSet> AddAssign<&SeaBasisBladeLow<T>> for SeaSparseMultivectorLowDimensionsWithMetric<T, I> {
	fn add_assign(&mut self, rhs: &SeaBasisBladeLow<T>) {
		self.multivector += rhs
	}
}

impl<T:SeaCoefficient, I:SeaBasisBladeSet> AddAssign<Self> for SeaSparseMultivectorLowDimensionsWithMetric<T, I> {
	fn add_assign(&mut self, rhs: Self) {
		*self += &rhs
	}
}

impl<T:SeaCoefficient, I:SeaBasisBladeSet> AddAssign<&Self> for SeaSparseMultivectorLowDimensionsWithMetric<T, I> {
	fn add_assign(&mut self, rhs: &Self) {
		assert_eq!(self.metric, rhs.metric);
		self.multivector += &rhs.multivector
	}
}

impl<T:SeaCoefficient, I:SeaBasisBladeSet> Add<Self> for SeaSparseMultivectorLowDimensionsWithMetric<T, I> {
    type Output = SeaSparseMultivectorLowDimensionsWithMetric<T, I>;

    fn add(self, rhs: Self) -> Self::Output {
		&self + &rhs
	}
}

impl<T:SeaCoefficient, I:SeaBasisBladeSet> Add<&Self> for SeaSparseMultivectorLowDimensionsWithMetric<T, I> {
    type Output = SeaSparseMultivectorLowDimensionsWithMetric<T, I>;

    fn add(self, rhs: &Self) -> Self::Output {
		&self + rhs
	}
}

impl<T:SeaCoefficient, I:SeaBasisBladeSet> Add<SeaSparseMultivectorLowDimensionsWithMetric<T, I>> for &SeaSparseMultivectorLowDimensionsWithMetric<T, I> {
    type Output = SeaSparseMultivectorLowDimensionsWithMetric<T, I>;

    fn add(self, rhs: SeaSparseMultivectorLowDimensionsWithMetric<T, I>) -> Self::Output {
		self + &rhs
	}
}

impl<T:SeaCoefficient, I:SeaBasisBladeSet> Add<Self> for &SeaSparseMultivectorLowDimensionsWithMetric<T, I> {
    type Output = SeaSparseMultivectorLowDimensionsWithMetric<T, I>;

    fn add(self, rhs: Self) -> Self::Output {
		assert_eq!(self.metric, rhs.metric);
		SeaSparseMultivectorLowDimensionsWithMetric {
			metric: self.metric.clone(),
			multivector: &self.multivector + &rhs.multivector
		}
	}
}

#[cfg(test)]
mod tests {
	use crate::tests::catch_unwind_silent;

    use super::*;

	#[test]
	fn add_basis_blades() {
		let mut a = SeaSparseMultivectorLowDimensionsWithMetric::<f32, u8>::from_str("3>-1 + e1 + 2e2 + 3e12 + 5e13 + 6e23 + 7e123").unwrap();
		a += &SeaBasisBladeLow::from_str("4e3").unwrap();
		assert_eq!(
			a,
			SeaSparseMultivectorLowDimensionsWithMetric::<f32, u8>::from_str("3>-1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123").unwrap()
		);

		let mut a = SeaSparseMultivectorLowDimensionsWithMetric::<f32, u8>::from_str("3>-1 + e1 + 2e2 + 3e12 + 5e13").unwrap();
		a += &SeaBasisBladeLow::from_str("4e3").unwrap();
		a += &SeaBasisBladeLow::from_str("7e123").unwrap();
		a += &SeaBasisBladeLow::from_str("6e23").unwrap();
		assert_eq!(
			a,
			SeaSparseMultivectorLowDimensionsWithMetric::<f32, u8>::from_str("3>-1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123").unwrap()
		);
	}

	#[test]
	fn add_mv() {
		let mut a = SeaSparseMultivectorLowDimensionsWithMetric::<f32, u8>::from_str("3>-1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123").unwrap();
		let b = SeaSparseMultivectorLowDimensionsWithMetric::<f32, u8>::from_str("3>-1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123").unwrap();
		assert_eq!(
			&a + &b,
			SeaSparseMultivectorLowDimensionsWithMetric::<f32, u8>::from_str("3>-2 + 2e1 + 4e2 + 6e12 + 8e3 + 10e13 + 12e23 + 14e123").unwrap()
		);

		a += &b;
		assert_eq!(
			a,
			SeaSparseMultivectorLowDimensionsWithMetric::<f32, u8>::from_str("3>-2 + 2e1 + 4e2 + 6e12 + 8e3 + 10e13 + 12e23 + 14e123").unwrap()
		);
	}
}
