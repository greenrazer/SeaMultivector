use super::*;

impl<T:SeaCoefficient, I:SeaBasisBladeSet> AddAssign<T> for SeaSparseMultivectorLowDimensions<T, I> {
    fn add_assign(&mut self, rhs: T) {
		*self += &SeaBasisBladeLow {
			index: SeaBasisBladeIndexLow::SCALAR,
			coefficient: rhs
		};
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeSet> Add<T> for SeaSparseMultivectorLowDimensions<T, I> {
    type Output = SeaSparseMultivectorLowDimensions<T, I>;

    fn add(self, rhs: T) -> Self::Output {
		&self + rhs
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeSet> Add<T> for &SeaSparseMultivectorLowDimensions<T, I> {
    type Output = SeaSparseMultivectorLowDimensions<T, I>;

    fn add(self, rhs: T) -> Self::Output {
		let mut out = self.clone();
		out += &SeaBasisBladeLow {
			index: SeaBasisBladeIndexLow::SCALAR,
			coefficient: rhs
		};
		out 
	}
}
