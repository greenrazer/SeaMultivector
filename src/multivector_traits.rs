use crate::coefficient::SeaCoefficient;

pub trait Multivector{
	fn new() -> Self;

	fn reversed(&self) -> Self;
	fn involuted(&self) -> Self;
	fn conjugated(&self) -> Self;
}

impl<T: SeaCoefficient> Multivector for T {
	fn new() -> Self {
		T::zero()
	}

	fn reversed(&self) -> Self {
		*self
	}

	fn involuted(&self) -> Self {
		*self
	}

	fn conjugated(&self) -> Self {
		*self
	}
}

pub trait MultivectorMut{
	fn reversion(&mut self);
	fn involution(&mut self);
	fn conjugation(&mut self);
}

impl<T: SeaCoefficient> MultivectorMut for T {
	fn reversion(&mut self){}
	fn involution(&mut self){}
	fn conjugation(&mut self){}
}