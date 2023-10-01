pub trait LeftContraction<Rhs> {
	type Output;

	fn left_contraction(&self, rhs: Rhs) -> Self::Output;
}

pub trait RightContraction<Rhs> {
	type Output;

	fn right_contraction(&self, rhs: Rhs) -> Self::Output;
}

pub trait InnerProduct<Rhs> {
	type Output;

	fn inner_product(&self, rhs: Rhs) -> Self::Output;
}

pub trait ScalarProduct<Rhs> {
	type Output;
	
	fn scalar_product(&self, rhs: Rhs) -> Self::Output;
}

pub trait OuterProduct<Rhs> {
	type Output;

	fn outer_product(&self, rhs: Rhs) -> Self::Output;
}

pub trait RegressiveProduct<Rhs> {
	type Output;

	fn regressive_product(&self, rhs: Rhs) -> Self::Output;
}

pub trait GeometricProduct<Rhs> {
	type Output;

	fn geometric_product(&self, rhs: Rhs) -> Self::Output;
}

pub trait SandwichProduct<Rhs> {
	type Output;

	fn sandwich_product(&self, rhs: Rhs) -> Self::Output;
}
