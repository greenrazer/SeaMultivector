use std::str::FromStr;

#[derive(Debug, PartialEq)]
enum SeaParseMetricError {
	MissingNumPositive,
	InvalidNumPositive,
	InvalidNumNegitive,
	InvalidNumZero,
	InvalidMetricLength
}

#[derive(Debug, PartialEq)]
struct SeaMetric {
	positive: u8,
	negitive: u8,
	zero: u8,
	dimension: u8
}
impl SeaMetric {
	fn new(p: u8, n: u8, z: u8) -> Self {
		Self { 
			positive: p, 
			negitive: n, 
			zero: z, 
			dimension: p+n+z
		}
	}
}
impl FromStr for SeaMetric {
	type Err = SeaParseMetricError;
	fn from_str(s: &str) -> Result<Self, Self::Err> {
		if s.trim().is_empty() {
			return Err(SeaParseMetricError::MissingNumPositive);
		}
		
		let mut split_iter = s.split(',');
		
		let p = match split_iter.next().unwrap().trim().parse::<u8>() {
			Ok(p) => p,
			Err(_) => return Err(SeaParseMetricError::InvalidNumPositive)
		};
		
		let n_raw = match split_iter.next() {
			Some(r) => r.trim(),
			None => return Ok(Self::new(p, 0, 0))
		};
		let n = match n_raw.parse::<u8>() {
			Ok(n) => n,
			Err(_) => return Err(SeaParseMetricError::InvalidNumNegitive)
		};

		let z_raw = match split_iter.next() {
			Some(r) => r.trim(),
			None => return Ok(Self::new(p, n, 0))
		};
		let z = match z_raw.parse::<u8>() {
			Ok(z) => z,
			Err(_) => return Err(SeaParseMetricError::InvalidNumZero)
		};

		match split_iter.next() {
			Some(_) => return Err(SeaParseMetricError::InvalidMetricLength),
			None => {}
		};

		Ok(Self::new(p, n, z))
	}
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn metric_from_str() {
		assert_eq!(SeaMetric::from_str("1").unwrap(), SeaMetric::new(1, 0, 0));
		assert_eq!(SeaMetric::from_str("1, 2").unwrap(), SeaMetric::new(1, 2, 0));
		assert_eq!(SeaMetric::from_str("1, 2, 3").unwrap(), SeaMetric::new(1, 2, 3));

		assert_eq!(
			SeaMetric::from_str("").unwrap_err(), 
			SeaParseMetricError::MissingNumPositive
		);
		assert_eq!(
			SeaMetric::from_str("a").unwrap_err(), 
			SeaParseMetricError::InvalidNumPositive
		);
		assert_eq!(
			SeaMetric::from_str("1,a").unwrap_err(), 
			SeaParseMetricError::InvalidNumNegitive
		);
		assert_eq!(
			SeaMetric::from_str("1,2,c").unwrap_err(), 
			SeaParseMetricError::InvalidNumZero
		);
		assert_eq!(
			SeaMetric::from_str("1,2,3,4").unwrap_err(), 
			SeaParseMetricError::InvalidMetricLength
		);

	}
}