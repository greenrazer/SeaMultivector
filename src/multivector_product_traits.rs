use crate::coefficient::SeaCoefficient;

pub trait LeftContraction<Rhs> {
	fn left_contraction(&self, rhs: Rhs) -> Self;
}

pub trait RightContraction<Rhs> {
	fn right_contraction(&self, rhs: Rhs) -> Self;
}

pub trait InnerProduct<Rhs> {
	fn inner_product(&self, rhs: Rhs) -> Self;
}

pub trait ScalarProduct<Rhs> {
	type Output;
	
	fn scalar_product(&self, rhs: Rhs) -> Self::Output;
}

pub trait OuterProduct<Rhs> {
	fn outer_product(&self, rhs: Rhs) -> Self;
}

pub trait RegressiveProduct<Rhs> {
	fn regressive_product(&self, rhs: Rhs) -> Self;
}

pub trait GeometricProduct<Rhs> {
	fn geometric_product(&self, rhs: Rhs) -> Self;
}
