use super::*;

impl<T:SeaCoefficient, I:SeaBasisBladeSet> MultivectorDual for SeaLowDimensionalMultivectorCustomMetric<T, I> {
	fn dual(&self) -> Self {
		Self {
			metric: self.metric.clone(),
			multivector: self.multivector.dual(&self.metric)
		}
	}

	fn undual(&self) -> Self {
		Self {
			metric: self.metric.clone(),
			multivector: self.multivector.undual(&self.metric)
		}
	}

	fn poincare_dual(&self) -> Self {
		Self {
			metric: self.metric.clone(),
			multivector: self.multivector.poincare_dual()
		}
	}

	fn right_complimented(&self) -> Self {
		Self {
			metric: self.metric.clone(),
			multivector: self.multivector.right_complimented()
		}
	}

	fn left_complimented(&self) -> Self {
		Self {
			metric: self.metric.clone(),
			multivector: self.multivector.left_complimented()
		}
	}
}

impl<T:SeaCoefficient, I:SeaBasisBladeSet> MultivectorDualMut for SeaLowDimensionalMultivectorCustomMetric<T, I> {
	fn dualization(&mut self) {
		self.multivector.dualization(&self.metric);
	}
	
	fn undualization(&mut self) {
		self.multivector.undualization(&self.metric);
	}
	
	fn poincare_dualization(&mut self) {
		self.multivector.poincare_dualization()
	}

	fn right_compliment(&mut self) {
		self.multivector.right_compliment()
	}

	fn left_compliment(&mut self) {
		self.multivector.left_compliment()
	}
}