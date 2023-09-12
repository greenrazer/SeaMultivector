use super::*;

// Inner Product
impl<T:SeaCoefficient, I:SeaBasisBladeSet> SeaSparseMultivectorLowDimensions<T, I> {
	pub fn inner_product(&self, rhs: &Self, metric: &SeaMetric) -> Self {
		let (num_metric_non_zeros, metric_non_zeros, metric_negs) = Self::metric_product_constants(metric);

		let output_basis_blade_bits = {
			self.basis_blade_bits.inner_product(
				(&rhs.basis_blade_bits, num_metric_non_zeros)
			)
		};

		let mut output_coefficients = unsafe {
			Self::create_coefficients_array(&output_basis_blade_bits)
		};

		for (out_idx, o) in output_basis_blade_bits.iter().enumerate() {

			let left_contraction_rhs_bbb = rhs.basis_blade_bits.left_contraction_with_basis_blade(&o, &metric_non_zeros);
			
			let mut dot_mask_left = self.basis_blade_bits.intersection(&left_contraction_rhs_bbb);
			
			let total_neg_mask = I::total_negations(&o, &metric_negs);

			output_coefficients[out_idx] = Self::special_dot_product(
				self,
				rhs,
				&dot_mask_left,
				&total_neg_mask,
				&o
			);

			if !o.is_zero() {
				let mut right_contraction_rhs_bbb = rhs.basis_blade_bits.right_contraction_with_basis_blade(&o, &metric_non_zeros);

				let mut dot_mask_right = self.basis_blade_bits.intersection(&right_contraction_rhs_bbb);

				output_coefficients[out_idx] += Self::special_dot_product(
					self,
					rhs,
					&dot_mask_right,
					&total_neg_mask,
					&o
				);
			}
		}

		Self {
			basis_blade_bits: output_basis_blade_bits,
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
		
		let a = SeaSparseMultivectorLowDimensions::<f32, u8>::from_str("1").unwrap();
		let b = SeaSparseMultivectorLowDimensions::<f32, u8>::from_str("3e12").unwrap();

		let out = SeaSparseMultivectorLowDimensions::<f32, u8>::from_str("3e12").unwrap();
		assert_eq!(
			a.inner_product(&b, &m),
			out
		);

		let a = SeaSparseMultivectorLowDimensions::<f32, u8>::from_str("e1 + e2").unwrap();
		let b = SeaSparseMultivectorLowDimensions::<f32, u8>::from_str("e2 + e3").unwrap();

		let out = SeaSparseMultivectorLowDimensions::<f32, u8>::from_str("1").unwrap();
		assert_eq!(
			a.inner_product(&b, &m),
			out
		);

		let a = SeaSparseMultivectorLowDimensions::<f32, u8>::from_str("1 + e1 + e3").unwrap();
		let b = SeaSparseMultivectorLowDimensions::<f32, u8>::from_str("     e1  + 3e12").unwrap();

		let out = SeaSparseMultivectorLowDimensions::<f32, u8>::from_str("1 + e1 + 3e2 + 3e12").unwrap();
		assert_eq!(
			a.inner_product(&b, &m),
			out
		);

		let a = SeaSparseMultivectorLowDimensions::<f32, u8>::from_str("e12").unwrap();
		let b = SeaSparseMultivectorLowDimensions::<f32, u8>::from_str("e2").unwrap();

		let out = SeaSparseMultivectorLowDimensions::<f32, u8>::from_str("e1").unwrap();
		assert_eq!(
			a.inner_product(&b, &m),
			out
		);

		let a = SeaSparseMultivectorLowDimensions::<f32, u8>::from_str("1 + 2e1 + 1e12 + 2e3").unwrap();
		let b = SeaSparseMultivectorLowDimensions::<f32, u8>::from_str("  + 1e2 + 2e13 + 1e23").unwrap();

		let out = SeaSparseMultivectorLowDimensions::<f32, u8>::from_str("-3e1 -e2 +4e3 +2e13 + e23").unwrap();
		assert_eq!(
			a.inner_product(&b, &m),
			out
		);

		let m = SeaMetric::new(4, 0, 0);

		let a = SeaSparseMultivectorLowDimensions::<f32, u16>::from_str(" + 1 + e1 + 2e12 + e13 + e234 ").unwrap();
		let b = SeaSparseMultivectorLowDimensions::<f32, u16>::from_str("  + 1e4 + 2e124 + 3e13 +4e34").unwrap();

		let out = SeaSparseMultivectorLowDimensions::<f32, u16>::from_str("-3 -4e2 +3e3 -3e4 +3e13 +e23 +2e24 +4e34 +2e124​").unwrap();
		assert_eq!(
			a.inner_product(&b, &m),
			out
		);
	}

	#[test]
	fn projective() {
		let m = SeaMetric::new(3, 0, 1);
		let a = SeaSparseMultivectorLowDimensions::<f32, u16>::from_str("  + e1 + 2e2 - e234 + e1234 ").unwrap();
		let b = SeaSparseMultivectorLowDimensions::<f32, u16>::from_str("  - 3e13 + 1e4 - 2e124  +4e34").unwrap();

		let out = SeaSparseMultivectorLowDimensions::<f32, u16>::from_str("-3e3 +4e14 -5e24").unwrap();
		assert_eq!(
			a.inner_product(&b, &m),
			out
		);

		let a = SeaSparseMultivectorLowDimensions::<f32, u16>::from_str(" + 1 + e1 + 2e12 + e13 + e234 ").unwrap();
		let b = SeaSparseMultivectorLowDimensions::<f32, u16>::from_str("  + 3e13 + 1e4 + 2e124  +4e34").unwrap();

		let out = SeaSparseMultivectorLowDimensions::<f32, u16>::from_str("-3 -3e4 +3e3 +2e24 +4e34 +3e13 +2e124").unwrap();
		assert_eq!(
			a.inner_product(&b, &m),
			out
		);
	}

	#[test]
	fn hyperbolic() {
		let m = SeaMetric::new(3, 1, 0);

		let a = SeaSparseMultivectorLowDimensions::<f32, u16>::from_str("1 + 2e1 + 1e12 + 2e3").unwrap();
		let b = SeaSparseMultivectorLowDimensions::<f32, u16>::from_str("  + 1e2 + 2e13 + 1e23").unwrap();

		let out = SeaSparseMultivectorLowDimensions::<f32, u16>::from_str("-3e1 -e2 +4e3 +2e13 +e23").unwrap();
		assert_eq!(
			a.inner_product(&b, &m),
			out
		);

		let a = SeaSparseMultivectorLowDimensions::<f32, u16>::from_str(" + e1 + 2e2 - e234 + e1234  ").unwrap();
		let b = SeaSparseMultivectorLowDimensions::<f32, u16>::from_str("  - 3e13 + 1e4 - 2e124  +4e34").unwrap();

		let out = SeaSparseMultivectorLowDimensions::<f32, u16>::from_str("-4e2 -5e3 +4e12 +4e14 +e23 -5e24 -e123").unwrap();
		assert_eq!(
			a.inner_product(&b, &m),
			out
		);

		let a = SeaSparseMultivectorLowDimensions::<f32, u16>::from_str(" + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e34 + 13e134 + 14e234 + 15e1234").unwrap();
		let b = SeaSparseMultivectorLowDimensions::<f32, u16>::from_str("  + 3e13 + 1e4 + 2e124  +4e34").unwrap();

		let out = SeaSparseMultivectorLowDimensions::<f32, u16>::from_str("47 +51e1 +49e2 +53e3 -28e4 +33e12 -10e13 -4e14 -14e23 +47e24 +4e34 -15e123 +2e124").unwrap();
		assert_eq!(
			a.inner_product(&b, &m),
			out
		);
	}

	#[test]
	fn other() {
		let m = SeaMetric::new(3, 2, 0);

		let a = SeaSparseMultivectorLowDimensions::<f32, u32>::from_str(" + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345").unwrap();
		let b = SeaSparseMultivectorLowDimensions::<f32, u32>::from_str("  + 3e13 + 1e4 + 2e124  +4e345").unwrap();

		let out = SeaSparseMultivectorLowDimensions::<f32, u32>::from_str(" -49 -53e1 -63e2 +3e3 -5e4 -87e12 +3e13 -4e14 +2e24 +14e35 -23e45 +2e124 +13e135 +14e235 +45e245 +4e345 +15e1235").unwrap();
		assert_eq!(
			a.inner_product(&b, &m),
			out
		);

		let m = SeaMetric::new(3, 1, 1);
		
		let a = SeaSparseMultivectorLowDimensions::<f32, u32>::from_str(" + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345").unwrap();
		let b = SeaSparseMultivectorLowDimensions::<f32, u32>::from_str("  + 3e13 + 1e4 + 2e124  + 4e345").unwrap();

		let out = SeaSparseMultivectorLowDimensions::<f32, u32>::from_str(" -1 -e1 -7e2 +3e3 -5e4 +14e35 -23e45 -27e12 +3e13 -4e14 +2e24 +13e135 +14e235 +45e245 +4e345 +2e124 +15e1235​").unwrap();
		assert_eq!(
			a.inner_product(&b, &m),
			out
		);
	}

	#[test]
	fn dimensions_6() {
			let m = SeaMetric::new(3, 3, 0);
	
			let a = SeaSparseMultivectorLowDimensions::<f32, u64>::from_str(" + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345 +e6 + 12e126 + 8e46").unwrap();
			let b = SeaSparseMultivectorLowDimensions::<f32, u64>::from_str("  + 3e13 + 1e4 + 2e124  +4e345 +83e1235 + 70e1245 + 4e1345 + 8e2345 + 12e123456").unwrap();
	
			let out = SeaSparseMultivectorLowDimensions::<f32, u64>::from_str("115+115e1 -27e2 +1053e3 +1240e4 +189e5 +188e6 -87e12 +3e13 -4e14 +202e15 +168e16 +2e24 -215e25 -156e26 -351e35 -301e45 +2e124 -228e125 -144e126 -185e135 -156e145 +33e235 +83e245 -120e345 -132e356 -84e456 +194e1235 +70e1245 +4e1345 -120e1356 -72e1456 +8e2345 +108e2356 +60e2456 -36e3456 +12e12345 +96e12356 +48e12456 -24e13456 +12e23456 +12e123456​").unwrap();
			// 32.862µs
			assert_eq!(
				a.inner_product(&b, &m),
				out
			);
	
			let a = SeaSparseMultivectorLowDimensions::<f32, u64>::from_str("1 + 2e1 + 1e12 + 2e3").unwrap();
			let b = SeaSparseMultivectorLowDimensions::<f32, u64>::from_str("  + 1e2 + 2e13 + 1e23").unwrap();
	
			let out = SeaSparseMultivectorLowDimensions::<f32, u64>::from_str("-3e1 -e2 +4e3 +2e13 +e23").unwrap();
			assert_eq!(
				a.inner_product(&b, &m),
				out
			);

			// let now = Instant::now();
			// let num_iters = 10_000;
			// let mut t = a.clone();
			// for _ in 0..num_iters {
			// 	t = t.inner_product(&b, &m);
			// }
			// let elapsed_time = now.elapsed();
			// println!("time time time {:?}", elapsed_time/num_iters);
		}
}