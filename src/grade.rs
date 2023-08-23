use std::ops::{Deref, DerefMut};

use crate::coefficient::SeaCoefficient;

#[derive(Debug, PartialEq, PartialOrd, Clone, Copy)]
pub struct SeaGrade(pub u8);
impl SeaGrade {
	pub fn is_reversion_negitive(&self) -> bool {
		if self.0 == 0{
			return false;
		}
		((self.0*(self.0-1)) >> 1)% 2 == 1
	}
	pub fn is_involution_negitive(&self) -> bool {
		self.0 % 2 == 1
	}
	pub fn is_conjugation_negitive(&self) -> bool {
		((self.0*(self.0+1)) >> 1) % 2 == 1
	}

	pub fn sign_reversion<T:SeaCoefficient>(&self) -> T {
		if self.is_reversion_negitive() {
			return -T::one();
		}
		T::one()
	}

	pub fn sign_involution<T:SeaCoefficient>(&self) -> T {
		if self.is_involution_negitive() {
			return -T::one();
		}
		T::one()
	}

	pub fn sign_conjugation<T:SeaCoefficient>(&self) -> T {
		if self.is_conjugation_negitive() {
			return -T::one();
		}
		T::one()
	}
}

impl Deref for SeaGrade {
	type Target = u8;
	fn deref(&self) -> &Self::Target {
		&self.0
	}
}
impl DerefMut for SeaGrade {
	fn deref_mut(&mut self) -> &mut Self::Target {
		&mut self.0
	}
}