use super::*;

impl<T:SeaCoefficient, I:SeaBasisBladeBits> DivAssign<T> for SeaMultivectorLowSparseCoefficents<T, I> {
	fn div_assign(&mut self, rhs: T) {
		for i in self.coefficients.iter_mut() {
            *i /= rhs;
        }
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeBits> Div<T> for SeaMultivectorLowSparseCoefficents<T, I> {
	type Output = SeaMultivectorLowSparseCoefficents<T, I>;
	fn div(self, rhs: T) -> Self::Output {
		&self / rhs
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeBits> Div<T> for &SeaMultivectorLowSparseCoefficents<T, I> {
	type Output = SeaMultivectorLowSparseCoefficents<T, I>;
	fn div(self, rhs: T) -> Self::Output {
		SeaMultivectorLowSparseCoefficents {
			basis_blade_bits: self.basis_blade_bits.clone(),
			coefficients: self.coefficients.iter().map(|x| *x/rhs).collect()
		}
	}
}
