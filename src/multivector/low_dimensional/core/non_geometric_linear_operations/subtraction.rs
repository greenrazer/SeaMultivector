use super::*;

impl<T:SeaCoefficient, I:SeaBasisBladeSet> SubAssign<&Self> for SeaLowDimensionalSparseMultivector<T, I> {
	fn sub_assign(&mut self, rhs: &Self) {
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
				self.coefficients[i] -= rhs.coefficients[rhs_idx];
				rhs_idx += 1;
			}
		}

		self.basis_blade_set = bb_out;
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeSet> Sub<Self> for &SeaLowDimensionalSparseMultivector<T, I> {
    type Output = SeaLowDimensionalSparseMultivector<T, I>;

    fn sub(self, rhs: Self) -> Self::Output {
		let mut out = self.clone();
		out -= rhs;
		out
	}
}
