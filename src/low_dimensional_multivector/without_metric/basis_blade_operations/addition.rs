use super::*;

impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> AddAssign<&SeaBasisBladeLow<T>> for SeaMultivectorLowSparseCoefficents<T, I> {
	#[allow(clippy::suspicious_op_assign_impl)]
	fn add_assign(&mut self, rhs: &SeaBasisBladeLow<T>) {
		let new_idx = self.basis_blade_bits.count_basis_blades_less_than(&rhs.index);
		if !self.basis_blade_bits.has(&rhs.index) {
			self.coefficients.insert(new_idx as usize, T::zero());
		}
		self.basis_blade_bits |= &rhs.index;
		self.coefficients[new_idx as usize] += rhs.coefficient;
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> Add<&SeaBasisBladeLow<T>> for &SeaMultivectorLowSparseCoefficents<T, I> {
	type Output = SeaMultivectorLowSparseCoefficents<T, I>;
	#[allow(clippy::suspicious_op_assign_impl)]
	fn add(self, rhs: &SeaBasisBladeLow<T>) -> Self::Output{
		let mut out = self.clone();
		out += rhs;
		out
	}
}
