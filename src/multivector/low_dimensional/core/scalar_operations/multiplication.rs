use super::*;

impl<T:SeaCoefficient, I:SeaBasisBladeSet> MulAssign<T> for SeaLowDimensionalSparseMultivector<T, I> {
	fn mul_assign(&mut self, rhs: T) {
		for i in self.coefficients.iter_mut() {
            *i *= rhs;
        }
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeSet> Mul<T> for SeaLowDimensionalSparseMultivector<T, I> {
	type Output = SeaLowDimensionalSparseMultivector<T, I>;
	fn mul(self, rhs: T) -> Self::Output {
		&self * rhs
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeSet> Mul<T> for &SeaLowDimensionalSparseMultivector<T, I> {
	type Output = SeaLowDimensionalSparseMultivector<T, I>;
	fn mul(self, rhs: T) -> Self::Output {
		SeaLowDimensionalSparseMultivector {
			basis_blade_set: self.basis_blade_set.clone(),
			coefficients: self.coefficients.iter().map(|x| *x*rhs).collect()
		}
	}
}

