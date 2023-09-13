use super::*;

// Scalar Product
impl<T:SeaCoefficient, I:SeaBasisBladeSet> SeaLowDimensionalSparseMultivector<T, I> {
	pub fn scalar_product(&self, rhs: &Self, metric: &SeaMetric) -> T {
		let (_, metric_non_zeros, metric_negs) = Self::metric_product_constants(metric);
			
		let mut dot_mask = self.basis_blade_set.intersection(&rhs.basis_blade_set);
		dot_mask.intersection_mut(&I::outer_product_zeros(&metric_non_zeros));
		
		let total_neg_mask = I::total_negations(&SeaBasisBladeIndexLow::SCALAR, &metric_negs);

		let mut out = T::zero();
		for i in dot_mask.iter() {
			let a_idx = self.basis_blade_set.count_basis_blades_less_than(&i) as usize;
			let b_idx = rhs.basis_blade_set.count_basis_blades_less_than(&i) as usize;
			let c_a = self.coefficients[a_idx];
			let c_b = rhs.coefficients[b_idx];
			let c = c_a*c_b;
			
			out += if total_neg_mask.has(&i) {-c} else {c};
		}

		out
	}
}

#[cfg(test)]
mod tests{
	use std::time::Instant;

	use crate::tests::catch_unwind_silent;

	use super::*;

	use num_traits::One;

	#[test]
	fn euclidean() {
		let m = SeaMetric::new(3, 0, 0);
		
		let a = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("1").unwrap();
		let b = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("3e12").unwrap();
		assert_eq!(
			a.scalar_product(&b, &m),
			0.0
		);

		let a = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("e1 + e2").unwrap();
		let b = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("e2 + e3").unwrap();
		assert_eq!(
			a.scalar_product(&b, &m),
			1.0
		);

		let a = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("1 + e1 + e3").unwrap();
		let b = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("     e1  + 3e12").unwrap();
		assert_eq!(
			a.scalar_product(&b, &m),
			1.0
		);

		let a = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("1 + 2e1 + 1e12 + 2e3").unwrap();
		let b = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("  + 1e2 + 2e13 + 1e23").unwrap();
		assert_eq!(
			a.scalar_product(&b, &m),
			0.0
		);

		let m = SeaMetric::new(4, 0, 0);

		let a = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str(" + 1 + e1 + 2e12 + e13 + e234 ").unwrap();
		let b = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str(" + 3e13 + 1e4 + 2e124 +4e34").unwrap();
		assert_eq!(
			a.scalar_product(&b, &m),
			-3.0
		);
	}

	#[test]
	fn projective() {
		let m = SeaMetric::new(3, 0, 1);

		let a = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("1 + 2e1 + 1e12 + 2e3").unwrap();
		let b = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("  + 1e2 + 2e13 + 1e23").unwrap();
		assert_eq!(
			a.scalar_product(&b, &m),
			0.0
		);

		let a = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str(" + 2e3 + 2e14 + 1e124 ").unwrap();
		let b = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("  + 1e24 + 2e134").unwrap();
		assert_eq!(
			a.scalar_product(&b, &m),
			0.0
		);

		let a = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str(" + 1 + e1 + 2e12 + e13 + e234 ").unwrap();
		let b = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("  + 3e13 + 1e4 + 2e124  +4e34").unwrap();
		assert_eq!(
			a.scalar_product(&b, &m),
			-3.0
		);

		let a = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("  + e1 + 2e2 - e234 + e1234 ").unwrap();
		let b = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("  - 3e13 + 1e4 - 2e124  +4e34").unwrap();
		assert_eq!(
			a.scalar_product(&b, &m),
			0.0
		);
	}

	#[test]
	fn hyperbolic() {
		let m = SeaMetric::new(3, 1, 0);

		let a = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("1 + 2e1 + 1e12 + 2e3").unwrap();
		let b = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("  + 1e2 + 2e13 + 1e23").unwrap();
		assert_eq!(
			a.scalar_product(&b, &m),
			0.0
		);

		let a = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str(" + e1 + 2e2 - e234 + e1234  ").unwrap();
		let b = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("  - 3e13 + 1e4 - 2e124  +4e34").unwrap();
		assert_eq!(
			a.scalar_product(&b, &m),
			0.0
		);

		let a = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str(" + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e34 + 13e134 + 14e234 + 15e1234").unwrap();
		let b = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("  + 3e13 + 1e4 + 2e124  +4e34").unwrap();
		assert_eq!(
			a.scalar_product(&b, &m),
			47.0
		);
	}

	#[test]
	fn other() {
		let m = SeaMetric::new(3, 2, 0);

		let a = SeaLowDimensionalSparseMultivector::<f32, u32>::from_str(" + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345").unwrap();
		let b = SeaLowDimensionalSparseMultivector::<f32, u32>::from_str("  + 3e13 + 1e4 + 2e124  +4e345").unwrap();

		assert_eq!(
			a.scalar_product(&b, &m),
			-49.0
		);

		let m = SeaMetric::new(3, 1, 1);
		
		let a = SeaLowDimensionalSparseMultivector::<f32, u32>::from_str(" + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345").unwrap();
		let b = SeaLowDimensionalSparseMultivector::<f32, u32>::from_str("  + 3e13 + 1e4 + 2e124  +4e345").unwrap();

		assert_eq!(
			a.scalar_product(&b, &m),
			-1.0
		);
	}

	#[test]
	fn dimensions_6() {
		let m = SeaMetric::new(3, 3, 0);

		let a = SeaLowDimensionalSparseMultivector::<f32, u64>::from_str(" + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345 +e6 + 12e126 + 8e46").unwrap();
		let b = SeaLowDimensionalSparseMultivector::<f32, u64>::from_str("  + 3e13 + 1e4 + 2e124  +4e345 +83e1235 + 70e1245 + 4e1345 + 8e2345 + 12e123456").unwrap();

		assert_eq!(
			a.scalar_product(&b, &m),
			115.0
		);

		let a = SeaLowDimensionalSparseMultivector::<f32, u64>::from_str("1 + 2e1 + 1e12 + 2e3").unwrap();
		let b = SeaLowDimensionalSparseMultivector::<f32, u64>::from_str("  + 1e2 + 2e13 + 1e23").unwrap();
		assert_eq!(
			a.scalar_product(&b, &m),
			0.0
		);
	}
}
