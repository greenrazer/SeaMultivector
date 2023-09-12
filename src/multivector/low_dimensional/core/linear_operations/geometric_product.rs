use super::*;

impl<T:SeaCoefficient, I:SeaBasisBladeSet> SeaLowDimensionalSparseMultivector<T, I>{ 
	pub fn geometric_product(&self, rhs: &Self, metric: &SeaMetric) -> Self {
		let (num_metric_non_zeros, metric_non_zeros, metric_negs) = Self::metric_product_constants(metric);

		let output_basis_blade_bits = {
			self.basis_blade_bits.geometric_product(
				(&rhs.basis_blade_bits, num_metric_non_zeros)
			)
		};

		let mut output_coefficients = unsafe {
			Self::create_coefficients_array(&output_basis_blade_bits)
		};

		for (out_idx, o) in output_basis_blade_bits.iter().enumerate() {
			
			let curr_rhs_bbb = rhs.basis_blade_bits.geometric_product_with_basis_blade(&o, &metric_non_zeros);

			let dot_mask = self.basis_blade_bits.intersection(&curr_rhs_bbb);
			
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
		
		let a = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("1").unwrap();
		let b = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("3e12").unwrap();

		let out = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("3e12").unwrap();
		//774ns
		assert_eq!(
			a.geometric_product(&b, &m),
			out
		);

		let a = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("e1 + e2").unwrap();
		let b = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("e2 + e3").unwrap();

		let out = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("1 + e12 + e13 + e23").unwrap();
		//3.382µs
		assert_eq!(
			a.geometric_product(&b, &m),
			out
		);

		let a = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("1 + e1 + e3").unwrap();
		let b = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("     e1  + 3e12").unwrap();

		let out = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("1 + e1 + 3e2 +3e12 -e13 + 3e123").unwrap();
		assert_eq!(
			a.geometric_product(&b, &m),
			out
		);

		let a = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("1 + 2e1 + 1e12 + 2e3").unwrap();
		let b = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("  + 1e2 + 2e13 + 1e23").unwrap();

		let out = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("-3e1 - e2 + 2e12 + 4e3 + 3e13 -3e23 + 2e123").unwrap();
		//7.542µs
		assert_eq!(
			a.geometric_product(&b, &m),
			out
		);

		let m = SeaMetric::new(4, 0, 0);

		let a = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str(" + 1 + e1 + 2e12 + e13 + e234 ").unwrap();
		let b = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("  + 1e4 + 2e124 + 3e13 +4e34").unwrap();

		let out = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("-3 -4e2 + 3e3 -3e4 +e13 +5e14 -5e23 + 2e24 + 4e34 + 7e124 + 5e134 + 2e234 + 8e1234").unwrap();
		//21.11µs
		//10.495µs
		assert_eq!(
			a.geometric_product(&b, &m),
			out
		);

		let _m = SeaMetric::new(4, 0, 0);

		let _a = SeaLowDimensionalSparseMultivector::<f32, u64>::from_str("1 + 1e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123").unwrap();
		let _b = SeaLowDimensionalSparseMultivector::<f32, u64>::from_str("1 + 1e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123").unwrap();

		// let now = Instant::now();
		// let num_iters = 10_000;
		// let mut t = a.clone();
		// for _ in 0..num_iters {
		// 	t = t.geometric_product(&b, &m);
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

		let a = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("1 + 2e1 + 1e12 + 2e3").unwrap();
		let b = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("  + 1e2 + 2e13 + 1e23").unwrap();

		let out = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("-3e1 - e2 + 2e12 + 4e3 + 3e13 -3e23 + 2e123").unwrap();
		assert_eq!(
			a.geometric_product(&b, &m),
			out
		);

		let a = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str(" + 2e3 + 2e14 + 1e124 ").unwrap();
		let b = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("  + 1e24 + 2e134").unwrap();

		let out = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("-4e14 -2e234").unwrap();
		assert_eq!(
			a.geometric_product(&b, &m),
			out
		);

		let a = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str(" + 1 + e1 + 2e12 + e13 + e234 ").unwrap();
		let b = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("  + 3e13 + 1e4 + 2e124  +4e34").unwrap();

		let out = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("-3 +3e3 + 3e13 -6e23 -3e4  + 5e14 +2e24 + 7e124 + 4e34  +5e134 + 2e234 + 8e1234").unwrap();
		//12.67µs
		assert_eq!(
			a.geometric_product(&b, &m),
			out
		);

		let a = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("  + e1 + 2e2 - e234 + e1234 ").unwrap();
		let b = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("  - 3e13 + 1e4 - 2e124  +4e34").unwrap();

		let out = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("-3e3 +6e123 + 5e14 - 3e24 +3e124 + 4e134 + 8e234 ").unwrap();
		//8.403µs
		assert_eq!(
			a.geometric_product(&b, &m),
			out
		);
	}

	#[test]
	fn hyperbolic() {
		let m = SeaMetric::new(3, 1, 0);

		let a = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("1 + 2e1 + 1e12 + 2e3").unwrap();
		let b = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("  + 1e2 + 2e13 + 1e23").unwrap();

		let out = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("-3e1 - e2 + 2e12 + 4e3 + 3e13 -3e23 + 2e123").unwrap();
		assert_eq!(
			a.geometric_product(&b, &m),
			out
		);

		let a = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str(" + e1 + 2e2 - e234 + e1234  ").unwrap();
		let b = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("  - 3e13 + 1e4 - 2e124  +4e34").unwrap();

		let out = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("-4e2 -5e3 + 4e12 + 2e13 + 5e14 + e23 -3e24 + 5e123 + 3e124 + 4e134 + 8e234").unwrap();
		assert_eq!(
			a.geometric_product(&b, &m),
			out
		);

		let a = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str(" + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e34 + 13e134 + 14e234 + 15e1234").unwrap();
		let b = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("  + 3e13 + 1e4 + 2e124  +4e34").unwrap();

		let out = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("47 +51e1 + 49e2 + 53e3 - 28e4 + 51e12 + 54e13 - 19e14 - 9e23 + 73e24 + 21e34 - e123 + 75e124 + 21e134 -9e234 -3e1234").unwrap();
		//21.547µs
		assert_eq!(
			a.geometric_product(&b, &m),
			out
		);
	}

	#[test]
	fn other() {
		let m = SeaMetric::new(3, 2, 0);

		let a = SeaLowDimensionalSparseMultivector::<f32, u32>::from_str(" + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345").unwrap();
		let b = SeaLowDimensionalSparseMultivector::<f32, u32>::from_str("  + 3e13 + 1e4 + 2e124  +4e345").unwrap();

		let out = SeaLowDimensionalSparseMultivector::<f32, u32>::from_str("-49 -53e1 -63e2 +3e3 -5e4 -69e12 + 3e13 - 3e14 - 9e23 + 4e24 + 17e34 + 14e35 -23e45 -6e123 + 5e124 +17e134 + 21e135 - 16e145 -17e234 + 80e235 + 69e245 + 4e345 - 15e1234 +83e1235 + 70e1245 + 4e1345 + 8e2345 + 12e12345").unwrap();
		assert_eq!(
			a.geometric_product(&b, &m),
			out
		);

		let m = SeaMetric::new(3, 1, 1);
		
		let a = SeaLowDimensionalSparseMultivector::<f32, u32>::from_str(" + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345").unwrap();
		let b = SeaLowDimensionalSparseMultivector::<f32, u32>::from_str("  + 3e13 + 1e4 + 2e124  +4e345").unwrap();

		let out = SeaLowDimensionalSparseMultivector::<f32, u32>::from_str("-1 -e1 -7e2 +3e3 - 5e4 +14e35 -23e45 -9e12 +3e13 - 3e14 -9e23 +4e24 +17e34 +21e135 - 16e145 +80e235 +69e245 + 4e345 -6e123 + 5e124 + 17e134 - 17e234 +83e1235 +70e1245 +4e1345 + 8e2345 -15e1234 +12e12345").unwrap();
		// 45.752µs
		// new 15.492µs
		assert_eq!(
			a.geometric_product(&b, &m),
			out
		);
	}

	#[test]
	fn dimensions_6() {
		let m = SeaMetric::new(3, 3, 0);

		let a = SeaLowDimensionalSparseMultivector::<f32, u64>::from_str(" + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345 +e6 + 12e126 + 8e46").unwrap();
		let b = SeaLowDimensionalSparseMultivector::<f32, u64>::from_str("  + 3e13 + 1e4 + 2e124  +4e345 +83e1235 + 70e1245 + 4e1345 + 8e2345 + 12e123456").unwrap();

		let out = SeaLowDimensionalSparseMultivector::<f32, u64>::from_str("115+115e1-27e2+1053e3+1240e4+189e5+188e6-21e12+983e13+1159e14+202e15+168e16-919e23-1075e24-215e25-156e26+17e34-351e35-301e45+23e46-846e123-991e124-228e125-128e126+17e134-265e135+3e136-228e145-17e234+143e235-36e236+135e245+303e345-1160e356-924e456-15e1234+230e1235+124e1245-14e1246-560e1256+442e1345+24e1346-152e1356-72e1456-393e2345+44e2356+60e2456-40e3456-360e12345+179e12356+118e12456+76e13456-28e23456-700e123456").unwrap();
		// 171.677µs
		// new: 58.328µs
		// with unsafe: 56.79µs
		// 49.433µs
		assert_eq!(
			a.geometric_product(&b, &m),
			out
		);

		let a = SeaLowDimensionalSparseMultivector::<f32, u64>::from_str("1 + 2e1 + 1e12 + 2e3").unwrap();
		let b = SeaLowDimensionalSparseMultivector::<f32, u64>::from_str("  + 1e2 + 2e13 + 1e23").unwrap();

		let out = SeaLowDimensionalSparseMultivector::<f32, u64>::from_str("-3e1 - e2 + 2e12 + 4e3 + 3e13 -3e23 + 2e123").unwrap();
		//13.991µs
		//3.334µs
		assert_eq!(
			a.geometric_product(&b, &m),
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