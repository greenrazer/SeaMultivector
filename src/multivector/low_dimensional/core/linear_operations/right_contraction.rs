use super::*;

// Right Contraction
impl<T:SeaCoefficient, I:SeaBasisBladeSet> SeaLowDimensionalSparseMultivector<T, I> {
	pub fn right_contraction(&self, rhs: &Self, metric: &SeaMetric) -> Self {
		let (num_metric_non_zeros, metric_non_zeros, metric_negs) = Self::metric_product_constants(metric);

		let output_basis_blade_set = {
			self.basis_blade_set.right_contraction(
				(&rhs.basis_blade_set, num_metric_non_zeros)
			)
		};

		let mut output_coefficients = unsafe {
			Self::create_coefficients_array(&output_basis_blade_set)
		};

		for (out_idx, o) in output_basis_blade_set.iter().enumerate() {
			let curr_rhs_bbb = rhs.basis_blade_set.right_contraction_with_basis_blade(&o,  &metric_non_zeros);
			
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

		let out = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("").unwrap();
		assert_eq!(
			a.right_contraction(&b, &m),
			out
		);

		let a = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("e1 + e2").unwrap();
		let b = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("e2 + e3").unwrap();

		let out = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("1").unwrap();
		assert_eq!(
			a.right_contraction(&b, &m),
			out
		);

		let a = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("1 + e1 + e3").unwrap();
		let b = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("     e1  + 3e12").unwrap();

		let out = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("1").unwrap();
		assert_eq!(
			a.right_contraction(&b, &m),
			out
		);

		let a = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("e12").unwrap();
		let b = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("e2").unwrap();

		let out = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("e1").unwrap();
		assert_eq!(
			a.right_contraction(&b, &m),
			out
		);

		let a = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("1 + 2e1 + 1e12 + 2e3").unwrap();
		let b = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("  + 1e2 + 2e13 + 1e23").unwrap();

		let out = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("e1").unwrap();
		assert_eq!(
			a.right_contraction(&b, &m),
			out
		);

		let m = SeaMetric::new(4, 0, 0);

		let a = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str(" + 1 + e1 + 2e12 + e13 + e234 ").unwrap();
		let b = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("  + 1e4 + 2e124 + 3e13 +4e34").unwrap();

		let out = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("-3 -4e2 + e23").unwrap();
		assert_eq!(
			a.right_contraction(&b, &m),
			out
		);
	}

	#[test]
	fn projective() {
		let m = SeaMetric::new(3, 0, 1);
		let a = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("  + e1 + 2e2 - e234 + e1234 ").unwrap();
		let b = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("  - 3e13 + 1e4 - 2e124  +4e34").unwrap();

		let out = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("-3e24").unwrap();
		assert_eq!(
			a.right_contraction(&b, &m),
			out
		);

		let a = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str(" + 1 + e1 + 2e12 + e13 + e234 ").unwrap();
		let b = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("  + 3e13 + 1e4 + 2e124  +4e34").unwrap();

		let out = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("-3").unwrap();
		assert_eq!(
			a.right_contraction(&b, &m),
			out
		);
	}

	#[test]
	fn hyperbolic() {
		let m = SeaMetric::new(3, 1, 0);

		let a = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("1 + 2e1 + 1e12 + 2e3").unwrap();
		let b = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("  + 1e2 + 2e13 + 1e23").unwrap();

		let out = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("e1").unwrap();
		assert_eq!(
			a.right_contraction(&b, &m),
			out
		);

		let a = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str(" + e1 + 2e2 - e234 + e1234  ").unwrap();
		let b = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("  - 3e13 + 1e4 - 2e124  +4e34").unwrap();

		let out = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("-4e2 - 2e3 +4e12 + e23 - 3e24 - e123").unwrap();
		assert_eq!(
			a.right_contraction(&b, &m),
			out
		);

		let a = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str(" + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e34 + 13e134 + 14e234 + 15e1234").unwrap();
		let b = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("  + 3e13 + 1e4 + 2e124  +4e34").unwrap();

		let out = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("47+43e1+67e2+18e3-39e4+49e12-13e13-14e23+45e24-15e123").unwrap();
		assert_eq!(
			a.right_contraction(&b, &m),
			out
		);
	}

	#[test]
	fn other() {
		let m = SeaMetric::new(3, 2, 0);

		let a = SeaLowDimensionalSparseMultivector::<f32, u32>::from_str(" + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345").unwrap();
		let b = SeaLowDimensionalSparseMultivector::<f32, u32>::from_str("  + 3e13 + 1e4 + 2e124  +4e345").unwrap();

		let out = SeaLowDimensionalSparseMultivector::<f32, u32>::from_str("-49-61e1-45e2-71e12-18e35-39e45+13e135+14e235+45e245+15e1235").unwrap();
		assert_eq!(
			a.right_contraction(&b, &m),
			out
		);

		let m = SeaMetric::new(3, 1, 1);
		
		let a = SeaLowDimensionalSparseMultivector::<f32, u32>::from_str(" + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345").unwrap();
		let b = SeaLowDimensionalSparseMultivector::<f32, u32>::from_str("  + 3e13 + 1e4 + 2e124  +4e345").unwrap();

		let out = SeaLowDimensionalSparseMultivector::<f32, u32>::from_str("-1 -9e1 +11e2 -18e35 -39e45 -11e12 +13e135 +14e235 +45e245 +15e1235").unwrap();
		assert_eq!(
			a.right_contraction(&b, &m),
			out
		);
	}

	#[test]
	fn dimensions_6() {
		let m = SeaMetric::new(3, 3, 0);

		let a = SeaLowDimensionalSparseMultivector::<f32, u64>::from_str(" + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345 +e6 + 12e126 + 8e46").unwrap();
		let b = SeaLowDimensionalSparseMultivector::<f32, u64>::from_str("  + 3e13 + 1e4 + 2e124  +4e345 +83e1235 + 70e1245 + 4e1345 + 8e2345 + 12e123456").unwrap();

		let out = SeaLowDimensionalSparseMultivector::<f32, u64>::from_str("115+59e1-105e2+1050e3+1245e4+8e6-71e12-18e35-39e45+13e135+14e235+45e245+15e1235\u{200B}").unwrap();
		assert_eq!(
			a.right_contraction(&b, &m),
			out
		);

		let a = SeaLowDimensionalSparseMultivector::<f32, u64>::from_str("1 + 2e1 + 1e12 + 2e3").unwrap();
		let b = SeaLowDimensionalSparseMultivector::<f32, u64>::from_str("  + 1e2 + 2e13 + 1e23").unwrap();

		let out = SeaLowDimensionalSparseMultivector::<f32, u64>::from_str("e1").unwrap();
		assert_eq!(
			a.right_contraction(&b, &m),
			out
		);
	}
}