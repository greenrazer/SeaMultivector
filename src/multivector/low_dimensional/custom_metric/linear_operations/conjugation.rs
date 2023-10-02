use super::*;

impl<T:SeaCoefficient, I:SeaBasisBladeSet> MultivectorConjugation for SeaLowDimensionalMultivectorCustomMetric<T, I> {
	fn reversion_conjugated(&self) -> Self {
		Self {
			metric: self.metric.clone(),
			multivector: self.multivector.reversion_conjugated()
		}
	}
	
	fn involution_conjugated(&self) -> Self {
		Self {
			metric: self.metric.clone(),
			multivector: self.multivector.involution_conjugated()
		}
	}
	
	fn clifford_conjugated(&self) -> Self {
		Self {
			metric: self.metric.clone(),
			multivector: self.multivector.clifford_conjugated()
		}
	}

	fn negated(&self) -> Self {
		Self {
			metric: self.metric.clone(),
			multivector: self.multivector.negated()
		}
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeSet> MultivectorConjugationMut for SeaLowDimensionalMultivectorCustomMetric<T, I> {
	fn reversion_conjugation(&mut self) {
		self.multivector.reversion_conjugation()
	}

	fn involution_conjugation(&mut self) {
		self.multivector.involution_conjugation()
	}

	fn clifford_conjugation(&mut self) {
		self.multivector.clifford_conjugation()
	}

	fn negate(&mut self) {
		self.multivector.negate()
	}
}

impl<T:SeaCoefficient, I:SeaBasisBladeSet> Neg for &SeaLowDimensionalMultivectorCustomMetric<T, I> {
	type Output = SeaLowDimensionalMultivectorCustomMetric<T, I>;

	fn neg(self) -> Self::Output {
		SeaLowDimensionalMultivectorCustomMetric {
			metric: self.metric.clone(),
			multivector: -&self.multivector
		}
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeSet> Neg for SeaLowDimensionalMultivectorCustomMetric<T, I> {
	type Output = SeaLowDimensionalMultivectorCustomMetric<T, I>;

	fn neg(self) -> Self::Output {
		-&self
	}
}