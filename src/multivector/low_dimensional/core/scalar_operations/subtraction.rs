use super::*;

impl<T:SeaCoefficient, I:SeaBasisBladeSet> SubAssign<T> for SeaLowDimensionalSparseMultivector<T, I> {
    fn sub_assign(&mut self, rhs: T) {
		*self += &SeaBasisBladeLow {
			index: SeaBasisBladeIndexLow::SCALAR,
			coefficient: -rhs
		};
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeSet> Sub<T> for SeaLowDimensionalSparseMultivector<T, I> {
    type Output = SeaLowDimensionalSparseMultivector<T, I>;

    fn sub(self, rhs: T) -> Self::Output {
		&self - rhs
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeSet> Sub<T> for &SeaLowDimensionalSparseMultivector<T, I> {
    type Output = SeaLowDimensionalSparseMultivector<T, I>;

    fn sub(self, rhs: T) -> Self::Output {
		let mut out = self.clone();
		out += &SeaBasisBladeLow {
			index: SeaBasisBladeIndexLow::SCALAR,
			coefficient: -rhs
		};
		out 
	}
}
