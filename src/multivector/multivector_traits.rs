use std::ops::{Mul, Add, Sub, Div, Neg, AddAssign, MulAssign, SubAssign, DivAssign, Shl, Shr, BitOr, BitXor, BitAnd, Not};

use crate::coefficient::SeaCoefficient;
use crate::multivector::multivector_product_traits::{GeometricProduct, OuterProduct, InnerProduct, LeftContraction, RightContraction, ScalarProduct, RegressiveProduct};


pub trait Multivector<T:SeaCoefficient>: 
		Neg + Not
		+ Add<T> + Mul<T> + Sub<T> + Div<T>
		+ Add<Self> + for<'a> Add<&'a Self>
		+ Sub<Self> + for<'a> Sub<&'a Self>
		+ LeftContraction<Self> + for<'a> LeftContraction<&'a Self>
		+ Shl<Self> + for<'a> Shl<&'a Self> 
		+ RightContraction<Self> + for<'a> RightContraction<&'a Self>
		+ Shr<Self> + for<'a> Shr<&'a Self>
		+ InnerProduct<Self> + for<'a> InnerProduct<&'a Self>
		+ BitOr<Self> + for<'a> BitOr<&'a Self>
		+ OuterProduct<Self> + for<'a> OuterProduct<&'a Self>
		+ BitXor<Self> + for<'a> BitXor<&'a Self>
		+ GeometricProduct<Self> + for<'a> GeometricProduct<&'a Self>
		+ Mul<Self> + for<'a> Mul<&'a Self>
		+ RegressiveProduct<Self> + for<'a> RegressiveProduct<&'a Self>
		+ BitAnd<Self> + for<'a> BitAnd<&'a Self> 
		+ ScalarProduct<Self> + for<'a> ScalarProduct<&'a Self>
			where
				Self: Sized,
				for<'a> &'a Self: Neg + Not
								+ Add<T> + Mul<T> + Sub<T> + Div<T>
								+ Add<Self> + Add<&'a Self> 
								+ Sub<Self> + Sub<&'a Self>
								+ Shl<Self> + Shl<&'a Self>
								+ Shr<Self> + Shr<&'a Self>
								+ BitOr<Self> + BitOr<&'a Self>
								+ BitXor<Self> + BitXor<&'a Self>
								+ Mul<Self> + Mul<&'a Self>
								+ BitAnd<Self> + BitAnd<&'a Self>
{
	type MultivectorIndexType;
	fn basis_blade_coefficient(&self, index: &Self::MultivectorIndexType) -> T;
	fn reversion_conjugated(&self) -> Self;
	fn involution_conjugated(&self) -> Self;
	fn clifford_conjugated(&self) -> Self;
	
	fn dual(&self) -> Self;
	fn undual(&self) -> Self;
	fn poincare_dual(&self) -> Self;
	fn right_complimented(&self) -> Self;
	fn left_complimented(&self) -> Self;

	fn scalar_norm(&self) -> T;
	fn scalar_normalized(&self) -> Self;

	fn exp(&self) -> Self;
}
pub trait MultivectorMut<T:SeaCoefficient>: 
	Multivector<T>
	+ AddAssign<T> + MulAssign<T> + SubAssign<T> + DivAssign<T>
	+ AddAssign<Self> + for<'a> AddAssign<&'a Self>
	+ SubAssign<Self> + for<'a> SubAssign<&'a Self>
		where 
			for<'a> &'a Self: Neg + Not
							+ Add<T> + Mul<T> + Sub<T> + Div<T>
							+ Add<Self> + Add<&'a Self> 
							+ Sub<Self> + Sub<&'a Self>
							+ Shl<Self> + Shl<&'a Self>
							+ Shr<Self> + Shr<&'a Self>
							+ BitOr<Self> + BitOr<&'a Self>
							+ BitXor<Self> + BitXor<&'a Self>
							+ Mul<Self> + Mul<&'a Self>
							+ BitAnd<Self> + BitAnd<&'a Self>
{
	fn reversion_conjugation(&mut self);
	fn involution_conjugation(&mut self);
	fn clifford_conjugation(&mut self);
	
	fn dualization(&mut self);
	fn undualization(&mut self);
	fn poincare_dualization(&mut self);
	fn right_compliment(&mut self);
	fn left_compliment(&mut self);

	fn scalar_normalize(&mut self);
}

pub trait Versor<T:SeaCoefficient>: 
	Multivector<T> 
		where 
			for<'a> &'a Self: Neg + Not
				+ Add<T> + Mul<T> + Sub<T> + Div<T>
				+ Add<Self> + Add<&'a Self> 
				+ Sub<Self> + Sub<&'a Self>
				+ Shl<Self> + Shl<&'a Self>
				+ Shr<Self> + Shr<&'a Self>
				+ BitOr<Self> + BitOr<&'a Self>
				+ BitXor<Self> + BitXor<&'a Self>
				+ Mul<Self> + Mul<&'a Self>
				+ BitAnd<Self> + BitAnd<&'a Self>
{
	fn inverse(&self) -> Self;
}

pub trait Motor<T:SeaCoefficient>: 
	Versor<T>
		where 
			for<'a> &'a Self: Neg + Not
				+ Add<T> + Mul<T> + Sub<T> + Div<T>
				+ Add<Self> + Add<&'a Self> 
				+ Sub<Self> + Sub<&'a Self>
				+ Shl<Self> + Shl<&'a Self>
				+ Shr<Self> + Shr<&'a Self>
				+ BitOr<Self> + BitOr<&'a Self>
				+ BitXor<Self> + BitXor<&'a Self>
				+ Mul<Self> + Mul<&'a Self>
				+ BitAnd<Self> + BitAnd<&'a Self>
{}

pub trait Flector<T:SeaCoefficient>: 
	Versor<T>
		where 
			for<'a> &'a Self: Neg + Not
				+ Add<T> + Mul<T> + Sub<T> + Div<T>
				+ Add<Self> + Add<&'a Self> 
				+ Sub<Self> + Sub<&'a Self>
				+ Shl<Self> + Shl<&'a Self>
				+ Shr<Self> + Shr<&'a Self>
				+ BitOr<Self> + BitOr<&'a Self>
				+ BitXor<Self> + BitXor<&'a Self>
				+ Mul<Self> + Mul<&'a Self>
				+ BitAnd<Self> + BitAnd<&'a Self>
{}