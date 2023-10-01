use crate::coefficient::SeaCoefficient;

pub trait MultivectorConjugation {
	fn reversion_conjugated(&self) -> Self;
	fn involution_conjugated(&self) -> Self;
	fn clifford_conjugated(&self) -> Self;
}

pub trait MultivectorConjugationMut {
	fn reversion_conjugation(&mut self);
	fn involution_conjugation(&mut self);
	fn clifford_conjugation(&mut self);
}

pub trait MultivectorDual {
	fn dual(&self) -> Self;
	fn undual(&self) -> Self;
	fn poincare_dual(&self) -> Self;
	fn right_complimented(&self) -> Self;
	fn left_complimented(&self) -> Self;
}

pub trait MultivectorDualMut {
	fn dualization(&mut self);
	fn undualization(&mut self);
	fn poincare_dualization(&mut self);
	fn right_compliment(&mut self);
	fn left_compliment(&mut self);
}

pub trait MultivectorNorm<T: SeaCoefficient> {
	fn scalar_norm(&self) -> T;
	fn scalar_normalized(&self) -> Self;
}

pub trait MultivectorNormMut {
	fn scalar_normalize(&mut self);
}