use super::*;

impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> AddAssign<&Self> for SeaMultivectorLowSparseCoefficents<T, I> {
	fn add_assign(&mut self, rhs: &Self) {
		let bb_out = &self.basis_blade_bits | &rhs.basis_blade_bits;

		let size_diff = bb_out.count_ones() as usize - self.coefficients.len();
		self.coefficients.splice(0..0, repeat(T::zero()).take(size_diff));

		let mut self_idx = 0;
		let mut rhs_idx = 0;
		for (i, idx) in bb_out.into_iter().enumerate() {
			if self.basis_blade_bits.has(&idx) {
				self.coefficients.swap(i, size_diff+self_idx);
				self_idx += 1;
			}
			
			if rhs.basis_blade_bits.has(&idx) {
				self.coefficients[i] += rhs.coefficients[rhs_idx];
				rhs_idx += 1;
			}
		}

		self.basis_blade_bits = bb_out;
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> Add<Self> for &SeaMultivectorLowSparseCoefficents<T, I> {
    type Output = SeaMultivectorLowSparseCoefficents<T, I>;

    fn add(self, rhs: Self) -> Self::Output {
		let mut out = self.clone();
		out += rhs;
		out
	}
}
