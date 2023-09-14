use super::*;

impl<T:SeaCoefficient, I:SeaBasisBladeSet> SeaLowDimensionalSparseMultivector<T, I>{
	pub fn sandwich_product_versor(&self, rhs: &Self, metric: &SeaMetric) -> Self {
		let (num_metric_non_zeros, metric_non_zeros, metric_negs) = Self::metric_product_constants(metric);

		let stage_1_output_basis_blade_set = {
			self.basis_blade_set.geometric_product(
				(&rhs.basis_blade_set, num_metric_non_zeros)
			)
		};

		let stage_2_output_basis_blade_set = {
			stage_1_output_basis_blade_set.geometric_product(
				(&self.basis_blade_set, num_metric_non_zeros)
			)
		};

		let mut stage_1_output_coefficients = unsafe {
			Self::create_coefficients_array(&stage_1_output_basis_blade_set)
		};

		let mut stage_2_output_coefficients = unsafe {
			Self::create_coefficients_array(&stage_2_output_basis_blade_set)
		};

		for (out_idx, o) in stage_1_output_basis_blade_set.iter().enumerate() {
			
			let curr_rhs_bbb = rhs.basis_blade_set.geometric_product_with_basis_blade(&o, &metric_non_zeros);

			let dot_mask = self.basis_blade_set.intersection(&curr_rhs_bbb);
			
			let total_neg_mask = I::total_negations(&o, &metric_negs);
			
			stage_1_output_coefficients[out_idx] = Self::special_dot_product_raw(
				&self.basis_blade_set,
				&self.coefficients,
				&rhs.basis_blade_set,
				&rhs.coefficients,
				&dot_mask,
				&total_neg_mask,
				&o
			);
		}

		for (out_idx, o) in stage_2_output_basis_blade_set.iter().enumerate() {
			
			let curr_rhs_bbb = self.basis_blade_set.geometric_product_with_basis_blade(&o, &metric_non_zeros);

			let dot_mask = stage_1_output_basis_blade_set.intersection(&curr_rhs_bbb);
			
			let mut total_neg_mask = I::total_negations(&o, &metric_negs);
			let bbbs_reversion_neg_swap = I::reversion_negations().geometric_product_with_basis_blade(&o, &metric_non_zeros);
			total_neg_mask.symmetric_difference_mut(&bbbs_reversion_neg_swap);
			
			stage_2_output_coefficients[out_idx] = Self::special_dot_product_raw(
				&stage_1_output_basis_blade_set,
				&stage_1_output_coefficients,
				&self.basis_blade_set,
				&self.coefficients,
				&dot_mask,
				&total_neg_mask,
				&o
			);
		}

		let mut output = Self {
			basis_blade_set: stage_2_output_basis_blade_set,
			coefficients: stage_2_output_coefficients
		};

		output.remove_zeros();

		output
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
		//774ns
		assert_eq!(
			a.sandwich_product_versor(&b, &m),
			out
		);

		let a = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("e1 + e2").unwrap();
		let b = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("e2 + e3").unwrap();

		let out = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("2e1 - 2e3").unwrap();
		//3.382µs
		assert_eq!(
			a.sandwich_product_versor(&b, &m),
			out
		);

		let a = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("1 + e1 + e3").unwrap();
		let b = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("     e1  + 3e12").unwrap();

		let out = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("2 + e1 + 2e3 + 3e12 + 6e23 + 6e123").unwrap();
		assert_eq!(
			a.sandwich_product_versor(&b, &m),
			out
		);

		let a = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("1 + 2e1 + 1e12 + 2e3").unwrap();
		let b = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("  + 1e2 + 2e13 + 1e23").unwrap();

		let out = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("4 +2e1  -8e2 +8e12 -14e13 -4e23 -4e123").unwrap();
		//7.542µs
		assert_eq!(
			a.sandwich_product_versor(&b, &m),
			out
		);

		let m = SeaMetric::new(4, 0, 0);

		let a = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str(" + 1 + e1 + 2e12 + e13 + e234 ").unwrap();
		let b = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("  + 1e4 + 2e124 + 3e13 +4e34").unwrap();

		let out = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("6e4 + 20e12 - 16e13 + 16e14 - 4e23 - 6e24 + 24e34 - 20e123 + 14e124 + 16e134 - 6e234").unwrap();
		//21.11µs
		//10.495µs
		assert_eq!(
			a.sandwich_product_versor(&b, &m),
			out
		);

		// let _m = SeaMetric::new(4, 0, 0);

		// let _a = SeaLowDimensionalSparseMultivector::<f32, u64>::from_str("1 + 1e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123").unwrap();
		// let _b = SeaLowDimensionalSparseMultivector::<f32, u64>::from_str("1 + 1e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123").unwrap();

		// let now = Instant::now();
		// let num_iters = 10_000;
		// let mut t = a.clone();
		// for _ in 0..num_iters {
		// 	t = t.sandwich_product_versor(&b, &m);
		// }
		// let elapsed_time = now.elapsed();
		//u8 8.612µs
		//u16 8.625µs
		//u32 5.222µs
		//u64 6.029µs
		//might as well only use 32 bit integers
		// println!("time time time {:?}", elapsed_time/num_iters);
	}

	#[test]
	fn projective() {
		let m = SeaMetric::new(3, 0, 1);

		let a = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("1 + 2e1 + 1e12 + 2e3").unwrap();
		let b = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("  + 1e2 + 2e13 + 1e23").unwrap();

		let out = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("4 +2e1  -8e2 +8e12 -14e13 -4e23 -4e123").unwrap();
		assert_eq!(
			a.sandwich_product_versor(&b, &m),
			out
		);

		let a = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str(" + 2e3 + 2e14 + 1e124 ").unwrap();
		let b = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("  + 1e24 + 2e134").unwrap();

		let out = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("+4e24 + 8e134").unwrap();
		assert_eq!(
			a.sandwich_product_versor(&b, &m),
			out
		);

		let a = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str(" + 1 + e1 + 2e12 + e13 + e234 ").unwrap();
		let b = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("  + 3e13 + 1e4 + 2e124  +4e34").unwrap();

		let out = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("5e4 + 16e14 - 6e24 + 20e34 + 12e12 - 9e13 - 12e23 + 16e124 + 16e134 - 6e234 - 12e123").unwrap();
		//12.67µs
		assert_eq!(
			a.sandwich_product_versor(&b, &m),
			out
		);

		let a = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("  + e1 + 2e2 - e234 + e1234 ").unwrap();
		let b = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("  - 3e13 + 1e4 - 2e124  +4e34").unwrap();

		let out = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("- 5e4 - 12e14 + 6e24 + 20e34 - 9e13 + 12e23 - 16e124 ").unwrap();
		//8.403µs
		assert_eq!(
			a.sandwich_product_versor(&b, &m),
			out
		);
	}

	#[test]
	fn hyperbolic() {
		let m = SeaMetric::new(3, 1, 0);

		let a = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("1 + 2e1 + 1e12 + 2e3").unwrap();
		let b = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("  + 1e2 + 2e13 + 1e23").unwrap();

		let out = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("4 + 2e1 - 8e2 + 8e12 - 14e13 - 4e23 - 4e123").unwrap();
		assert_eq!(
			a.sandwich_product_versor(&b, &m),
			out
		);

		let a = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str(" + e1 + 2e2 - e234 + e1234  ").unwrap();
		let b = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("  - 3e13 + 1e4 - 2e124  +4e34").unwrap();

		let out = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str(" - 4e3 - 5e4 + 8e12 - 5e13 - 12e14 + 20e23 + 10e24 + 12e34 - 8e123 - 20e124 + 8e134").unwrap();
		assert_eq!(
			a.sandwich_product_versor(&b, &m),
			out
		);

		let a = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str(" + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e34 + 13e134 + 14e234 + 15e1234").unwrap();
		let b = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("  + 3e13 + 1e4 + 2e124  +4e34").unwrap();

		let out = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str(" - 744 - 738e1 + 100e2 - 40e3 + e4 - 164e12 - 109e13 + 320e14 - 2454e23 + 2078e24 - 2610e34 - 2452e123 + 2102e124 - 2608e134 + 146e234 + 14e1234").unwrap();
		assert_eq!(
			a.sandwich_product_versor(&b, &m),
			out
		);
	}
}