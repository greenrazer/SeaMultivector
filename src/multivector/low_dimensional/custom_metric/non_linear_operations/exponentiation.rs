use super::*;


impl<T:SeaCoefficient, I:SeaBasisBladeSet> MultivectorExponentiation for SeaLowDimensionalMultivectorCustomMetric<T, I> {
	fn exp(&self) -> Self {
		Self {
			metric: self.metric.clone(),
			multivector: self.multivector.exp(&self.metric)
		}
	}
}
