use super::*;

impl<T:SeaCoefficient> MultivectorDual for SeaMultivectorPGA3D<T> {
	fn poincare_dual(&self) -> Self {
		Self(self.0.poincare_dual())
	}

	fn right_complimented(&self) -> Self {
		Self(self.0.right_complimented())
	}

	fn left_complimented(&self) -> Self {
		Self(self.0.left_complimented())
	}

	fn dual(&self) -> Self {
		Self(self.0.dual(&Self::METRIC))
	}

	fn undual(&self) -> Self {
		Self(self.0.undual(&Self::METRIC))
	}
}

impl<T:SeaCoefficient> MultivectorDualMut for SeaMultivectorPGA3D<T> {
	fn poincare_dualization(&mut self) {
		self.0.poincare_dualization()
	}

	fn right_compliment(&mut self) {
		self.0.right_compliment()
	}

	fn left_compliment(&mut self) {
		self.0.left_compliment()
	}

	fn dualization(&mut self) {
		self.0.dualization(&Self::METRIC)
	}

	fn undualization(&mut self) {
		self.0.undualization(&Self::METRIC)
	}
}