use crate::coefficient::SeaCoefficient;


pub trait Multivector{
	fn reversion_conjugated(&self) -> Self;
	fn involution_conjugated(&self) -> Self;
	fn clifford_conjugated(&self) -> Self;
	
	fn poincare_dual(&self) -> Self;
	fn right_complimented(&self) -> Self;
	fn left_complimented(&self) -> Self;
}
pub trait MultivectorMut: Multivector{
	fn reversion_conjugation(&mut self);
	fn involution_conjugation(&mut self);
	fn clifford_conjugation(&mut self);
	
	fn poincare_dualization(&mut self);
	fn right_compliment(&mut self);
	fn left_compliment(&mut self);
}

pub trait Versor: Multivector {}
pub trait VersorMut: Versor + MultivectorMut {}

pub trait Motor: Multivector {}
pub trait MotorMut: Motor + MultivectorMut {}

pub trait Glide: Multivector {}
pub trait GlideMut: Glide + MultivectorMut {}