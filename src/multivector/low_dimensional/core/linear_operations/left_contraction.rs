use super::*;

impl<T:SeaCoefficient, I:SeaBasisBladeSet> SeaLowDimensionalSparseMultivector<T, I> {
	pub fn left_contraction(&self, rhs: &Self, metric: &SeaMetric) -> Self {
		let (num_metric_non_zeros, metric_non_zeros, metric_negs) = Self::metric_product_constants(metric);
		
		let output_basis_blade_set = {
			self.basis_blade_set.left_contraction(
				(&rhs.basis_blade_set, num_metric_non_zeros)
			)
		};

		let mut output_coefficients = unsafe {
			Self::create_coefficients_array(&output_basis_blade_set)
		};

		for (out_idx, o) in output_basis_blade_set.iter().enumerate() {
			
			let curr_rhs_bbb = rhs.basis_blade_set.left_contraction_with_basis_blade(&o,  &metric_non_zeros);
			
			let dot_mask = self.basis_blade_set.intersection(&curr_rhs_bbb);
			
			let total_neg_mask = I::total_negations(&o, &metric_negs);

			output_coefficients[out_idx] = Self::special_dot_product(
				self,
				rhs,
				&dot_mask,
				&total_neg_mask,
				&o
			);
		}

		Self {
			basis_blade_set: output_basis_blade_set,
			coefficients: output_coefficients
		}
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

		let out = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("3e12").unwrap();
		assert_eq!(
			a.left_contraction(&b, &m),
			out
		);

		let a = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("e1 + e2").unwrap();
		let b = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("e2 + e3").unwrap();

		let out = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("1").unwrap();
		assert_eq!(
			a.left_contraction(&b, &m),
			out
		);

		let a = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("1 + e1 + e3").unwrap();
		let b = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("     e1  + 3e12").unwrap();

		let out = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("1 + e1 + 3e2 +3e12").unwrap();
		assert_eq!(
			a.left_contraction(&b, &m),
			out
		);

		let a = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("1 + 2e1 + 1e12 + 2e3").unwrap();
		let b = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("  + 1e2 + 2e13 + 1e23").unwrap();

		let out = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("-4e1 - e2 + 4e3 + 2e13 +1e23").unwrap();
		assert_eq!(
			a.left_contraction(&b, &m),
			out
		);

		let m = SeaMetric::new(4, 0, 0);

		let a = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str(" + 1 + e1 + 2e12 + e13 + e234 ").unwrap();
		let b = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("  + 1e4 + 2e124 + 3e13 +4e34").unwrap();

		let out = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("-3 + 3e3 -3e4 + 3e13 + 2e24 + 4e34 + 2e124").unwrap();
		assert_eq!(
			a.left_contraction(&b, &m),
			out
		);
	}

	#[test]
	fn projective() {
		let m = SeaMetric::new(3, 0, 1);
		let a = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("  + e1 + 2e2 - e234 + e1234 ").unwrap();
		let b = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("  - 3e13 + 1e4 - 2e124  +4e34").unwrap();

		let out = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("-3e3 +4e14 -2e24").unwrap();
		assert_eq!(
			a.left_contraction(&b, &m),
			out
		);

		let a = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str(" + 1 + e1 + 2e12 + e13 + e234 ").unwrap();
		let b = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("  + 3e13 + 1e4 + 2e124  +4e34").unwrap();

		let out = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("-3 -3e4 + 3e3 +2e24 + 4e34 + 3e13 + 2e124").unwrap();
		assert_eq!(
			a.left_contraction(&b, &m),
			out
		);
	}

	#[test]
	fn hyperbolic() {
		let m = SeaMetric::new(3, 1, 0);

		let a = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("1 + 2e1 + 1e12 + 2e3").unwrap();
		let b = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("  + 1e2 + 2e13 + 1e23").unwrap();

		let out = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("-4e1 -e2 +4e3 +2e13 +e23").unwrap();
		assert_eq!(
			a.left_contraction(&b, &m),
			out
		);

		let a = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str(" + e1 + 2e2 - e234 + e1234  ").unwrap();
		let b = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("  - 3e13 + 1e4 - 2e124  +4e34").unwrap();

		let out = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("-3e3 + 4e14 -2e24").unwrap();
		assert_eq!(
			a.left_contraction(&b, &m),
			out
		);

		let a = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str(" + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e34 + 13e134 + 14e234 + 15e1234").unwrap();
		let b = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("  + 3e13 + 1e4 + 2e124  +4e34").unwrap();

		let out = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("47 +8e1 -18e2 +35e3 + 11e4 -16e12 + 3e13 -4e14 +2e24 +4e34 + 2e124").unwrap();
		//4.371µs
		assert_eq!(
			a.left_contraction(&b, &m),
			out
		);
	}

	#[test]
	fn other() {
		let m = SeaMetric::new(3, 2, 0);

		let a = SeaLowDimensionalSparseMultivector::<f32, u32>::from_str(" + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345").unwrap();
		let b = SeaLowDimensionalSparseMultivector::<f32, u32>::from_str("  + 3e13 + 1e4 + 2e124  +4e345").unwrap();

		let out = SeaLowDimensionalSparseMultivector::<f32, u32>::from_str("-49 +8e1 -18e2 +3e3 -5e4 -16e12 +3e13 -4e14 +2e24 +32e35 +16e45 +2e124 + 4e345").unwrap();
		assert_eq!(
			a.left_contraction(&b, &m),
			out
		);

		let m = SeaMetric::new(3, 1, 1);
		
		let a = SeaLowDimensionalSparseMultivector::<f32, u32>::from_str(" + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345").unwrap();
		let b = SeaLowDimensionalSparseMultivector::<f32, u32>::from_str("  + 3e13 + 1e4 + 2e124  +4e345").unwrap();

		let out = SeaLowDimensionalSparseMultivector::<f32, u32>::from_str("-1 +8e1 -18e2 + 3e3 -5e4 +32e35 +16e45 -16e12 + 3e13 -4e14 +2e24 +4e345 +2e124").unwrap();
		//3.867µs
		assert_eq!(
			a.left_contraction(&b, &m),
			out
		);
	}

	#[test]
	fn dimensions_6() {
		let m = SeaMetric::new(3, 3, 0);

		let a = SeaLowDimensionalSparseMultivector::<f32, u64>::from_str(" + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345 +e6 + 12e126 + 8e46").unwrap();
		let b = SeaLowDimensionalSparseMultivector::<f32, u64>::from_str("  + 3e13 + 1e4 + 2e124  +4e345 +83e1235 + 70e1245 + 4e1345 + 8e2345 + 12e123456").unwrap();

		let out = SeaLowDimensionalSparseMultivector::<f32, u64>::from_str("115+56e1+78e2+3e3-5e4+189e5+180e6-16e12+3e13-4e14+202e15+168e16+2e24-215e25-156e26-333e35-262e45+2e124-228e125-144e126-198e135-156e145+19e235+38e245-120e345-132e356-84e456+179e1235+70e1245+4e1345-120e1356-72e1456+8e2345+108e2356+60e2456-36e3456+12e12345+96e12356+48e12456-24e13456+12e23456+12e123456").unwrap();
		//25.833µs
		//23.094µs
		assert_eq!(
			a.left_contraction(&b, &m),
			out
		);

		let a = SeaLowDimensionalSparseMultivector::<f32, u64>::from_str("1 + 2e1 + 1e12 + 2e3").unwrap();
		let b = SeaLowDimensionalSparseMultivector::<f32, u64>::from_str("  + 1e2 + 2e13 + 1e23").unwrap();

		let out = SeaLowDimensionalSparseMultivector::<f32, u64>::from_str("-4e1 -e2 +4e3 +2e13 +e23").unwrap();
		//2.081µs
		assert_eq!(
			a.left_contraction(&b, &m),
			out
		);

		// let now = Instant::now();
		// let num_iters = 10_000;
		// let mut t = a.clone();
		// for _ in 0..num_iters {
		// 	t = t.left_contraction(&b, &m);
		// }
		// let elapsed_time = now.elapsed();
		// println!("time time time {:?}", elapsed_time/num_iters);
	}
}