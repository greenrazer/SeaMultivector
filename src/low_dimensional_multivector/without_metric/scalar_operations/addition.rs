use super::*;

impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> AddAssign<T> for SeaMultivectorLowSparseCoefficents<T, I> {
    fn add_assign(&mut self, rhs: T) {
		*self += &SeaBasisBladeLow {
			index: SeaBasisBladeIndexLow::ZERO,
			coefficient: rhs
		};
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> Add<T> for SeaMultivectorLowSparseCoefficents<T, I> {
    type Output = SeaMultivectorLowSparseCoefficents<T, I>;

    fn add(self, rhs: T) -> Self::Output {
		&self + rhs
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> Add<T> for &SeaMultivectorLowSparseCoefficents<T, I> {
    type Output = SeaMultivectorLowSparseCoefficents<T, I>;

    fn add(self, rhs: T) -> Self::Output {
		let mut out = self.clone();
		out += &SeaBasisBladeLow {
			index: SeaBasisBladeIndexLow::ZERO,
			coefficient: rhs
		};
		out 
	}
}
