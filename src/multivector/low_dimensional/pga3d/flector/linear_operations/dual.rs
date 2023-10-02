use super::*;

impl<T:SeaCoefficient> MultivectorDual for SeaFlectorPGA3D<T> {
	fn poincare_dual(&self) -> Self {
		Self((**self).poincare_dual())
	}

	fn right_complimented(&self) -> Self {
		Self((**self).right_complimented())
	}

	fn left_complimented(&self) -> Self {
		Self((**self).left_complimented())
	}

	fn dual(&self) -> Self {
		Self((**self).dual())
	}

	fn undual(&self) -> Self {
		Self((**self).undual())
	}
}

impl<T:SeaCoefficient> MultivectorDualMut for SeaFlectorPGA3D<T> {
	fn poincare_dualization(&mut self) {
		(**self).poincare_dualization()
	}

	fn right_compliment(&mut self) {
		(**self).right_compliment()
	}

	fn left_compliment(&mut self) {
		(**self).left_compliment()
	}

	fn dualization(&mut self) {
		(**self).dualization()
	}

	fn undualization(&mut self) {
		(**self).undualization()
	}
}

impl<T:SeaCoefficient> Not for SeaFlectorPGA3D<T> {
	type Output = SeaFlectorPGA3D<T>;

	fn not(self) -> Self::Output {
		self.dual()
	}
}
impl<T:SeaCoefficient> Not for &SeaFlectorPGA3D<T> {
	type Output = SeaFlectorPGA3D<T>;

	fn not(self) -> Self::Output {
		self.dual()
	}
}