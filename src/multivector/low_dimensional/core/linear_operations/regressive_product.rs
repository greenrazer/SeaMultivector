use super::*;

// Regressive Product
impl<T:SeaCoefficient, I:SeaBasisBladeSet> SeaLowDimensionalSparseMultivector<T, I> {
	pub fn regressive_product(&self, rhs: &Self, metric: &SeaMetric) -> Self {
		let (a_dual, b_dual) = (self.dual(metric), rhs.dual(metric));

		let mut out_dual = &a_dual ^ &b_dual;
		
		out_dual.undualization(metric);
		
		out_dual
	}
}

#[cfg(test)]
mod tests{
	use std::time::Instant;

	use crate::tests::catch_unwind_silent;

	use super::*;

	#[test]
	fn euclidean() {
		let m = SeaMetric::new(3, 0, 0);
		
		let a = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("1").unwrap();
		let b = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("3e12").unwrap();

		let out = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("").unwrap();
		assert_eq!(
			a.regressive_product(&b, &m),
			out
		);

		let a = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("e1 + e2").unwrap();
		let b = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("e2 + e3").unwrap();

		let out = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("").unwrap();
		assert_eq!(
			a.regressive_product(&b, &m),
			out
		);

		let a = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("1 + e1 + e3").unwrap();
		let b = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("     e1  + 3e12").unwrap();

		let out = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("3").unwrap();
		assert_eq!(
			a.regressive_product(&b, &m),
			out
		);

		let a = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("1 + 2e1 + 1e12 + 2e3").unwrap();
		let b = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("  + 1e2 + 2e13 + 1e23").unwrap();

		let out = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("2 -2e1 -e2").unwrap();
		assert_eq!(
			a.regressive_product(&b, &m),
			out
		);

		let m = SeaMetric::new(4, 0, 0);

		let a = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str(" + 1 + e1 + 2e12 + e13 + e234 ").unwrap();
		let b = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("  + 1e4 + 2e124 + 3e13 +4e34").unwrap();

		let out = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("8 -2e1 +3e3 +2e24").unwrap();
		assert_eq!(
			a.regressive_product(&b, &m),
			out
		);
	}

	#[test]
	fn projective() {
		let m = SeaMetric::new(3, 0, 1);

		let a = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("1 + 2e1 + 1e12 + 2e3").unwrap();
		let b = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("  + 1e2 + 2e13 + 1e23").unwrap();

		let out = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("").unwrap();
		assert_eq!(
			a.regressive_product(&b, &m),
			out
		);

		let a = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str(" + 2e3 + 2e14 + 1e124 ").unwrap();
		let b = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("  + 1e24 + 2e134").unwrap();

		let out = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("2e14").unwrap();
		assert_eq!(
			a.regressive_product(&b, &m),
			out
		);

		let a = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str(" + 1 + e1 + 2e12 + e13 + e234 ").unwrap();
		let b = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("  + 3e13 + 1e4 + 2e124  +4e34").unwrap();

		let out = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("+8 -2e1 +3e3 -2e24").unwrap();
		assert_eq!(
			a.regressive_product(&b, &m),
			out
		);

		let a = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("  + e1 + 2e2 - e234 + e1234 ").unwrap();
		let b = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("  - 3e13 + 1e4 - 2e124  +4e34").unwrap();

		let out = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("+e4 +3e3 -2e24 +4e34 -3e13 -2e124").unwrap();
		assert_eq!(
			a.regressive_product(&b, &m),
			out
		);

		let m = SeaMetric::new(3, 0, 2);

		let a = SeaLowDimensionalSparseMultivector::<f32, u32>::from_str(" + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345").unwrap();
		let b = SeaLowDimensionalSparseMultivector::<f32, u32>::from_str("  + 3e12 + 1e4 + 2e124  +4e345").unwrap();
		
		let out = SeaLowDimensionalSparseMultivector::<f32, u32>::from_str("48 +83e4 -39e1 -42e2 +28e3 +26e14 +28e24 +45e12 +60e345 +30e124").unwrap();
		assert_eq!(
			a.regressive_product(&b, &m),
			out
		);
	}

	#[test]
	fn hyperbolic() {
		let m = SeaMetric::new(3, 1, 0);

		let a = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("1 + 2e1 + 1e12 + 2e3").unwrap();
		let b = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("  + 1e2 + 2e13 + 1e23").unwrap();

		let out = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("").unwrap();
		assert_eq!(
			a.regressive_product(&b, &m),
			out
		);

		let a = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str(" + e1 + 2e2 - e234 + e1234  ").unwrap();
		let b = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("  - 3e13 + 1e4 - 2e124  +4e34").unwrap();

		let out = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("3e3 +e4 -3e13 +2e24 +4e34 -2e124").unwrap();
		assert_eq!(
			a.regressive_product(&b, &m),
			out
		);

		let a = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str(" + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e34 + 13e134 + 14e234 + 15e1234").unwrap();
		let b = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("  + 3e13 + 1e4 + 2e124  +4e34").unwrap();

		let out = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str(" -33 -43e1 -12e2 +70e3 +83e4 -14e12 +45e13 +26e14 +28e24 +60e34 +30e124").unwrap();
		assert_eq!(
			a.regressive_product(&b, &m),
			out
		);
	}

	#[test]
	fn other() {
		let m = SeaMetric::new(3, 2, 0);

		let a = SeaLowDimensionalSparseMultivector::<f32, u32>::from_str(" + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345").unwrap();
		let b = SeaLowDimensionalSparseMultivector::<f32, u32>::from_str("  + 3e13 + 1e4 + 2e124  +4e345").unwrap();

		let out = SeaLowDimensionalSparseMultivector::<f32, u32>::from_str("12 +70e3 +83e4 +45e13 +26e14 +28e24 +30e124 +60e345").unwrap();
		assert_eq!(
			a.regressive_product(&b, &m),
			out
		);

		let m = SeaMetric::new(3, 1, 1);
		
		let a = SeaLowDimensionalSparseMultivector::<f32, u32>::from_str(" + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345").unwrap();
		let b = SeaLowDimensionalSparseMultivector::<f32, u32>::from_str("  + 3e13 + 1e4 + 2e124  +4e345").unwrap();

		let out = SeaLowDimensionalSparseMultivector::<f32, u32>::from_str("12 -14e3 +83e4 +45e13 +26e14 +28e24 +60e345 +30e124").unwrap();
		assert_eq!(
			a.regressive_product(&b, &m),
			out
		);
	}

	#[test]
	fn dimensions_6() {
		let m = SeaMetric::new(3, 3, 0);

		let a = SeaLowDimensionalSparseMultivector::<f32, u64>::from_str(" + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345 +e6 + 12e126 + 8e46").unwrap();
		let b = SeaLowDimensionalSparseMultivector::<f32, u64>::from_str("  + 3e13 + 1e4 + 2e124  +4e345 +83e1235 + 70e1245 + 4e1345 + 8e2345 + 12e123456").unwrap();

		let out = SeaLowDimensionalSparseMultivector::<f32, u64>::from_str(" -604 -36e1 -72e2 +48e3 +96e4 +12e6 +36e12 +60e13 +108e14 +72e23 +120e24 +96e46 +84e123 +132e124 +144e126 +144e345 +156e1345 +168e2345 +180e12345").unwrap();
		assert_eq!(
			a.regressive_product(&b, &m),
			out
		);

		let a = SeaLowDimensionalSparseMultivector::<f32, u64>::from_str("1 + 2e1 + 1e12 + 2e3").unwrap();
		let b = SeaLowDimensionalSparseMultivector::<f32, u64>::from_str("  + 1e2 + 2e13 + 1e23").unwrap();

		let out = SeaLowDimensionalSparseMultivector::<f32, u64>::from_str("").unwrap();
		assert_eq!(
			a.regressive_product(&b, &m),
			out
		);

		// let now = Instant::now();
		// let num_iters = 10_000;
		// let mut t = a.clone();
		// for _ in 0..num_iters {
		// 	t = t.geometric_product(&b, &m);
		// }
		// let elapsed_time = now.elapsed();
		// println!("time time time {:?}", elapsed_time/num_iters);
	}
}