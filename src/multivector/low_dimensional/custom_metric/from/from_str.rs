use super::*;

impl<T:SeaCoefficient, I:SeaBasisBladeSet> FromStr for SeaLowDimensionalMultivectorCustomMetric<T, I> {
	type Err = ParseMultivectorError;
	
	fn from_str(s: &str) -> Result<Self, Self::Err> {
		let (metric_str, mv_str) = match s.split_once('>') {
			Some(r) => r,
			None => return Err(ParseMultivectorError::MissingMetric)
		};
		let metric = match  SeaMetric::from_str(metric_str){
			Ok(v) => v,
			Err(e) => return Err(ParseMultivectorError::InvalidMetric(e))
		};
		if metric.dimensions != I::DIMENSIONS {
			return Err(ParseMultivectorError::InvalidMetricDimensions);
		}
		let mv = SeaLowDimensionalSparseMultivector::from_str(mv_str)?;

		Ok(Self {
			metric,
			multivector: mv 
		})
	}
}