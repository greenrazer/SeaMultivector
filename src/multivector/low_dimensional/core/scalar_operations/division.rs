use super::*;

impl<T:SeaCoefficient, I:SeaBasisBladeSet> DivAssign<T> for SeaLowDimensionalSparseMultivector<T, I> {
	fn div_assign(&mut self, rhs: T) {
		for i in self.coefficients.iter_mut() {
            *i /= rhs;
        }
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeSet> Div<T> for SeaLowDimensionalSparseMultivector<T, I> {
	type Output = SeaLowDimensionalSparseMultivector<T, I>;
	fn div(self, rhs: T) -> Self::Output {
		&self / rhs
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeSet> Div<T> for &SeaLowDimensionalSparseMultivector<T, I> {
	type Output = SeaLowDimensionalSparseMultivector<T, I>;
	fn div(self, rhs: T) -> Self::Output {
		SeaLowDimensionalSparseMultivector {
			basis_blade_bits: self.basis_blade_bits.clone(),
			coefficients: self.coefficients.iter().map(|x| *x/rhs).collect()
		}
	}
}
