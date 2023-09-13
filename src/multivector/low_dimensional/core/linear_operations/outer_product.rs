use super::*;

impl<T:SeaCoefficient, I:SeaBasisBladeSet> OuterProduct<&Self> for SeaLowDimensionalSparseMultivector<T, I> {
	type Output = Self;

	fn outer_product(&self, rhs: &Self) -> Self::Output {
		let output_basis_blade_set = {
			self.basis_blade_set.outer_product(&rhs.basis_blade_set)
		};

		let mut output_coefficients = unsafe {
			Self::create_coefficients_array(&output_basis_blade_set)
		};

		let reverse_rhs_bbbs = rhs.basis_blade_set.dual();

		for (out_idx, o) in output_basis_blade_set.iter().enumerate() {
			
			let curr_rhs_bbb = &reverse_rhs_bbbs.outer_product_with_basis_blade(&o);
			
			let dot_mask = self.basis_blade_set.intersection(curr_rhs_bbb);
			
			let cannonical_neg_mask = I::cannonical_reordering_negations(&o);

			output_coefficients[out_idx] = Self::special_dot_product(
				self,
				rhs,
				&dot_mask,
				&cannonical_neg_mask,
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

	// #[test]
	// fn thing() {

	// 	let a = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("-1 + 0.8682e14 - 0.4924e24 -1e34").unwrap();
	// 	let b = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("e3").unwrap();

	// 	println!("{:?}", a.outer_product(&b));
	// }

	#[test]
	fn works() {
		let a = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("1").unwrap();
		let b = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("3e12").unwrap();

		let out = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("3e12").unwrap();
		assert_eq!(
			a.outer_product(&b),
			out
		);

		let a = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("e1 + e2").unwrap();
		let b = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("e2 + e3").unwrap();

		let out = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("e12 + e13 + e23").unwrap();
		assert_eq!(
			a.outer_product(&b),
			out
		);

		let a = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("1 + e1 + e3").unwrap();
		let b = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("     e1  + 3e12").unwrap();

		let out = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("e1 + 3e12 - e13 + 3e123").unwrap();
		assert_eq!(
			a.outer_product(&b),
			out
		);

		let a = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("e12").unwrap();
		let b = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("e2").unwrap();

		let out = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("").unwrap();
		assert_eq!(
			a.outer_product(&b),
			out
		);

		let a = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("1 + 2e1 + 1e12 + 2e3").unwrap();
		let b = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("  + 1e2 + 2e13 + 1e23").unwrap();

		let out = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("e2 +2e12 +2e13 -e23 +2e123").unwrap();
		assert_eq!(
			a.outer_product(&b),
			out
		);

		let a = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str(" + 1 + e1 + 2e12 + e13 + e234 ").unwrap();
		let b = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("  + 1e4 + 2e124 + 3e13 +4e34").unwrap();

		let out = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("e4 +3e13 +e14 +4e34 +4e124 +5e134 +8e1234").unwrap();
		assert_eq!(
			a.outer_product(&b),
			out
		);

		let a = SeaLowDimensionalSparseMultivector::<f32, u32>::from_str(" + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345").unwrap();
		let b = SeaLowDimensionalSparseMultivector::<f32, u32>::from_str("  + 3e13 + 1e4 + 2e124  +4e345").unwrap();

		let out = SeaLowDimensionalSparseMultivector::<f32, u32>::from_str("e4 +3e13 +e14 +2e24 +4e34 -6e123 +5e124 +29e134 +6e234 +4e345 -15e1234 +4e1345 +8e2345 +12e12345").unwrap();
		assert_eq!(
			a.outer_product(&b),
			out
		);

		let a = SeaLowDimensionalSparseMultivector::<f32, u32>::from_str("15 +14e1 -13e2 +7e45 -11e35 +12e12 +6e145 -5e245 +3e345 -10e135 +9e235 +4e1245 -2e1345 +e2345 -8e1235 +e12345").unwrap();
		let b = SeaLowDimensionalSparseMultivector::<f32, u32>::from_str("-2e35 +4e12 +3e345 -e1235").unwrap();

		let out = SeaLowDimensionalSparseMultivector::<f32, u32>::from_str("-30e35 +60e12 +45e345 -28e135 +26e235 +28e1245 +42e1345 -39e2345 -83e1235 +48e12345").unwrap();
		assert_eq!(
			a.outer_product(&b),
			out
		);

		let a = SeaLowDimensionalSparseMultivector::<f32, u64>::from_str(" + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345 +e6 + 12e126 + 8e46").unwrap();
		let b = SeaLowDimensionalSparseMultivector::<f32, u64>::from_str("  + 3e13 + 1e4 + 2e124  +4e345 +83e1235 + 70e1245 + 4e1345 + 8e2345 + 12e123456").unwrap();

		let out = SeaLowDimensionalSparseMultivector::<f32, u64>::from_str("e4 +3e13 +e14 +2e24 +4e34 -e46 -6e123 +5e124 +29e134 +3e136 +6e234 +4e345 -15e1234 +83e1235 +70e1245 -14e1246 +8e1345 +24e1346 +16e2345 -4e3456 -372e12345 +83e12356 +70e12456 +4e13456 +8e23456 -700e123456").unwrap();
		//6.69µs
		assert_eq!(
			a.outer_product(&b),
			out
		);

		// let now = Instant::now();
		// let num_iters = 10_000;
		// let mut t = a.clone();
		// for _ in 0..num_iters {
		// 	t = a.outer_product(&b);
		// }
		// let elapsed_time = now.elapsed();
		// println!("time time time {:?}", elapsed_time/num_iters);
	}
}