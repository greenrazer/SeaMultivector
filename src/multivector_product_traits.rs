use crate::coefficient::SeaCoefficient;

pub trait LeftContraction<Rhs=Self> {
	type Output;

	fn left_contraction(&self, rhs: Rhs) -> Self::Output;
}

pub trait RightContraction<Rhs=Self> {
	type Output;

	fn right_contraction(&self, rhs: Rhs) -> Self::Output;
}

pub trait InnerProduct<Rhs=Self> {
	type Output;

	fn inner_product(&self, rhs: Rhs) -> Self::Output;
}

pub trait ScalarProduct<Rhs=Self> {
	type Output;

	fn scalar_product(&self, rhs: Rhs) -> Self::Output;
}

pub trait OuterProduct<Rhs=Self> {
	type Output;

	fn outer_product(&self, rhs: Rhs) -> Self::Output;
}

pub trait RegressiveProduct<Rhs=Self> {
	type Output;

	fn regressive_product(&self, rhs: Rhs) -> Self::Output;
}

pub trait GeometricProduct<Rhs=Self> {
	type Output;

	fn geometric_product(&self, rhs: Rhs) -> Self::Output;
}
