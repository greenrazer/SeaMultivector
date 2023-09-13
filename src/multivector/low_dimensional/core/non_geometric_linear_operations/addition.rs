use super::*;

impl<T:SeaCoefficient, I:SeaBasisBladeSet> AddAssign<&Self> for SeaLowDimensionalSparseMultivector<T, I> {
	fn add_assign(&mut self, rhs: &Self) {
		let bb_out = self.basis_blade_set.union(&rhs.basis_blade_set);

		let size_diff = bb_out.count_basis_blades() as usize - self.coefficients.len();
		self.coefficients.splice(0..0, repeat(T::zero()).take(size_diff));

		let mut self_idx = 0;
		let mut rhs_idx = 0;
		for (i, idx) in bb_out.iter().enumerate() {
			if self.basis_blade_set.has(&idx) {
				self.coefficients.swap(i, size_diff+self_idx);
				self_idx += 1;
			}
			
			if rhs.basis_blade_set.has(&idx) {
				self.coefficients[i] += rhs.coefficients[rhs_idx];
				rhs_idx += 1;
			}
		}

		self.basis_blade_set = bb_out;
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeSet> Add<Self> for &SeaLowDimensionalSparseMultivector<T, I> {
    type Output = SeaLowDimensionalSparseMultivector<T, I>;

    fn add(self, rhs: Self) -> Self::Output {
		let mut out = self.clone();
		out += rhs;
		out
	}
}

#[cfg(test)]
mod tests{
	use std::time::Instant;

	use crate::tests::catch_unwind_silent;

	use super::*;

	#[test]
	fn add() {
		let mut a = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("-1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123").unwrap();
		let b = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("-1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123").unwrap();
		assert_eq!(
			&a + &b,
			SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("-2 + 2e1 + 4e2 + 6e12 + 8e3 + 10e13 + 12e23 + 14e123").unwrap()
		);

		a += &b;
		assert_eq!(
			a,
			SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("-2 + 2e1 + 4e2 + 6e12 + 8e3 + 10e13 + 12e23 + 14e123").unwrap()
		);

		let mut a = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("-1 + e1 + 2e2 + 3e12 + 5e13 + 6e23 + 7e123").unwrap();
		a += &SeaBasisBladeLow::from_str("4e3").unwrap();
		assert_eq!(
			a,
			SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("-1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123").unwrap()
		);

		let mut a = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("-1 + e1 + 2e2 + 3e12 + 5e13").unwrap();
		a += &SeaBasisBladeLow::from_str("4e3").unwrap();
		a += &SeaBasisBladeLow::from_str("7e123").unwrap();
		a += &SeaBasisBladeLow::from_str("6e23").unwrap();
		assert_eq!(
			a,
			SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("-1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123").unwrap()
		);
	}
}
