use super::*;

impl<T:SeaCoefficient, I:SeaBasisBladeSet> MulAssign<T> for SeaSparseMultivectorLowDimensions<T, I> {
	fn mul_assign(&mut self, rhs: T) {
		for i in self.coefficients.iter_mut() {
            *i *= rhs;
        }
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeSet> Mul<T> for SeaSparseMultivectorLowDimensions<T, I> {
	type Output = SeaSparseMultivectorLowDimensions<T, I>;
	fn mul(self, rhs: T) -> Self::Output {
		&self * rhs
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeSet> Mul<T> for &SeaSparseMultivectorLowDimensions<T, I> {
	type Output = SeaSparseMultivectorLowDimensions<T, I>;
	fn mul(self, rhs: T) -> Self::Output {
		SeaSparseMultivectorLowDimensions {
			basis_blade_bits: self.basis_blade_bits.clone(),
			coefficients: self.coefficients.iter().map(|x| *x*rhs).collect()
		}
	}
}

