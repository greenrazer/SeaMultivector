use super::*;

impl<T:SeaCoefficient, I:SeaBasisBladeSet> MultivectorNorm<T> for SeaLowDimensionalMultivectorCustomMetric<T, I> {
	fn scalar_norm(&self) -> T {
		self.multivector.scalar_norm(&self.metric)
	}

	fn scalar_normalized(&self) -> Self {
		Self { 
			metric: self.metric.clone(), 
			multivector: self.multivector.scalar_normalized(&self.metric) 
		}
	}
}
impl<T:SeaCoefficient, I:SeaBasisBladeSet> MultivectorNormMut for SeaLowDimensionalMultivectorCustomMetric<T, I> {
	fn scalar_normalize(&mut self) {
		self.multivector.scalar_normalize(&self.metric);
	}
}