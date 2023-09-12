use super::*;

//Exp
impl<T:SeaCoefficient, I:SeaBasisBladeSet> SeaMultivectorLowSparseCoefficents<T, I> {
	pub fn exp(&self, metric: &SeaMetric) -> Self {
		
		let mut next_self = self.geometric_product(self, metric);
		
		if next_self.is_scalar() {
			let alpha = next_self.coefficients[0];
			return if alpha < T::zero() {
				let a = (-alpha).sqrt();
				let (sin_a, cos_a) = a.sin_cos();
				self*sin_a/a + cos_a
			}
			else if alpha == T::zero() {
				self + T::one()
			}
			else {
				let a = alpha.sqrt();
				self*a.sinh()/a + a.cosh()
			};
		}

		let s = self.scalar_norm(metric);
		if s < T::from(10_u8).unwrap() {
			let mut curr = self + T::one();
			curr += &(&next_self/T::factorial(2));
			for i in 3_u64..=14 {
				next_self = next_self.geometric_product(self, metric);
				curr += &(&next_self/T::factorial(i));
			}
			curr
		}
		else {
			// exp(A) = pow(exp(A/s), s)
			// where s is a power of 2 for speed
			let s = s.abs().to_u32().unwrap().next_power_of_two() >> 1;
			let norm_self = self/T::from(s).unwrap();
		
			let mut next_self = norm_self.geometric_product(&norm_self, metric);
			let mut curr = &norm_self + T::one();
			curr += &(&next_self/T::factorial(2));
			for i in 3_u64..=14 {
				next_self = next_self.geometric_product(&norm_self, metric);
				curr += &(&next_self/T::factorial(i));
			}
			
			for _ in 0..=s.trailing_zeros() {
				curr = curr.geometric_product(&curr, metric);
			} 
	
			curr
		}
	}
}

#[cfg(test)]
mod tests{
	use std::time::Instant;

	use crate::tests::catch_unwind_silent;

	use super::*;

	#[test]
	fn exp() {
		let m = SeaMetric::new(3, 0, 0);

		// let a = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("3e12").unwrap();
		// let out = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("-0.9899924994 + 0.1411200017e12").unwrap();
		// for (a_bb, o_bb) in a.exp(&m).iter().zip(out.iter()) {
		// 	assert_eq!(a_bb.index, o_bb.index);
		// 	let abs_diff = (a_bb.coefficient - o_bb.coefficient).abs();
		// 	assert!(abs_diff < 1e-6, "left != right(epsilon={})\n{:?}\n{:?}", 1e-6, a.exp(&m), out);
		// }

		let a = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("3e12 + 2e13 + 3e23").unwrap();
		let out = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("-0.021971453 -0.6394477487e12 -0.4262984991e13 -0.6394477487e23").unwrap();
		for (a_bb, o_bb) in a.exp(&m).iter().zip(out.iter()) {
			assert_eq!(a_bb.index, o_bb.index);
			let abs_diff = (a_bb.coefficient - o_bb.coefficient).abs();
			assert!(abs_diff < 1e-6, "left != right(epsilon={})\n{:?}\n{:?}", 1e-6, a.exp(&m), out);
		}

		let a = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("e1 + e2").unwrap();
		let out = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("2.1781836 + 1.3682989e1 + 1.3682989e2").unwrap();
		for (a_bb, o_bb) in a.exp(&m).iter().zip(out.iter()) {
			assert_eq!(a_bb.index, o_bb.index);
			let abs_diff = (a_bb.coefficient - o_bb.coefficient).abs();
			assert!(abs_diff < 1e-6, "left != right(epsilon={})\n{:?}\n{:?}", 1e-6, a.exp(&m), out);
		}

		//unstable
		let a = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str("1 + e1 + e2 + e12 + e3 + e13 + e23 + e123").unwrap();
		let out = SeaMultivectorLowSparseCoefficents::<f32, u8>::from_str(" -1.0374261 -2.0374262e1 +3.3594806e2 -2.0374262e3 +3.3594806e12 -2.0374262e13 +3.3594806e23 +3.3594806e123").unwrap();
		for (a_bb, o_bb) in a.exp(&m).iter().zip(out.iter()) {
			assert_eq!(a_bb.index, o_bb.index);
			let abs_diff = (a_bb.coefficient - o_bb.coefficient).abs();
			assert!(abs_diff < 1e-2, "left != right(epsilon={})\n{:?}\n{:?}", 1e-2, a.exp(&m), out);
		}

		let m = SeaMetric::new(3, 0, 1);

		let a = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("3e12").unwrap();
		let out = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("-0.9899924994 + 0.1411200017e12").unwrap();
		for (a_bb, o_bb) in a.exp(&m).iter().zip(out.iter()) {
			assert_eq!(a_bb.index, o_bb.index);
			let abs_diff = (a_bb.coefficient - o_bb.coefficient).abs();
			assert!(abs_diff < 1e-6, "left != right(epsilon={})\n{:?}\n{:?}", 1e-6, a.exp(&m), out);
		}

		let a = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("3e34").unwrap();
		let out = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("1 + 3e34").unwrap();
		for (a_bb, o_bb) in a.exp(&m).iter().zip(out.iter()) {
			assert_eq!(a_bb.index, o_bb.index);
			let abs_diff = (a_bb.coefficient - o_bb.coefficient).abs();
			assert!(abs_diff < 1e-6, "left != right(epsilon={})\n{:?}\n{:?}", 1e-6, a.exp(&m), out);
		}

		let a = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("e1 + e4").unwrap();
		let out = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("1.5430807 + 1.1752012e1 + 1.1752012e4").unwrap();
		for (a_bb, o_bb) in a.exp(&m).iter().zip(out.iter()) {
			assert_eq!(a_bb.index, o_bb.index);
			let abs_diff = (a_bb.coefficient - o_bb.coefficient).abs();
			assert!(abs_diff < 1e-6, "left != right(epsilon={})\n{:?}\n{:?}", 1e-6, a.exp(&m), out);
		}

		let a = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("e12 + 2e13 + 3e23").unwrap();
		let out = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("-0.8252990842 -0.1509213299e12 -0.3018426597e13 -0.4527639747e23").unwrap();
		for (a_bb, o_bb) in a.exp(&m).iter().zip(out.iter()) {
			assert_eq!(a_bb.index, o_bb.index);
			let abs_diff = (a_bb.coefficient - o_bb.coefficient).abs();
			assert!(abs_diff < 1e-6, "left != right(epsilon={})\n{:?}\n{:?}", 1e-6, a.exp(&m), out);
		}

		let a = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("4e14 + 5e24 + 6e34").unwrap();
		let out = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("1 +4e14 +5e24 +6e34").unwrap();
		for (a_bb, o_bb) in a.exp(&m).iter().zip(out.iter()) {
			assert_eq!(a_bb.index, o_bb.index);
			let abs_diff = (a_bb.coefficient - o_bb.coefficient).abs();
			assert!(abs_diff < 1e-6, "left != right(epsilon={})\n{:?}\n{:?}", 1e-6, a.exp(&m), out);
		}

		//unstable
		let a = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("e12 + 2e13 + 3e23 + 4e14 + 5e24 + 6e34").unwrap();
		let out = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("-0.8252991 -1.7576293e14 +0.015276484e24 -1.2898171e34 -0.15084463e12 -0.30168927e13 -0.452534e23 -1.2073706388e1234").unwrap();
		for (a_bb, o_bb) in a.exp(&m).iter().zip(out.iter()) {
			assert_eq!(a_bb.index, o_bb.index);
			let abs_diff = (a_bb.coefficient - o_bb.coefficient).abs();
			assert!(abs_diff < 1e-3, "left != right(epsilon={})\n{:?}\n{:?}", 1e-3, a.exp(&m), out);
		}

		//unstable
		let a = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("1 + e1 + e2 + e12 + e3 + e13 + e23 + e123 + e4 + e14 + e24 + e124 + e34 + e134 + e234 + e1234").unwrap();
		let out = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("-1.0374261 -2.0374262e1 +3.3594806e2 +3.3594806e12 -2.0374262e3 -2.0374262e13 +3.3594806e23 +3.3594806e123 -7.434334e4 -7.434334e14 +3.3594806e24 +3.3594806e124 +3.3594806e34 +3.3594806e134 +3.3594806e234 +3.3594806e1234").unwrap();
		for (a_bb, o_bb) in a.exp(&m).iter().zip(out.iter()) {
			assert_eq!(a_bb.index, o_bb.index);
			let abs_diff = (a_bb.coefficient - o_bb.coefficient).abs();
			assert!(abs_diff < 1e-3, "left != right(epsilon={})\n{:?}\n{:?}", 1e-3, a.exp(&m), out);
		}

		//unstable
		let a = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + e123 + e4 + e14 + e24 + 11e124 + 12e34 + 13e134 + 14e234 + 15e1234").unwrap();
		let out = SeaMultivectorLowSparseCoefficents::<f32, u16>::from_str("-5.4965935 -5.299999e1 -4.537388e2 -8.054262e12 -13.710886e3 -16.514507e13 -18.604914e23 -22.163525e123 + 165.71802e4 + 134.27966e14 -33.44229e24 -73.07824e124 -53.009445e34 -118.40819e134 -298.4169e234 -302.72876e1234").unwrap();
		for (a_bb, o_bb) in a.exp(&m).iter().zip(out.iter()) {
			assert_eq!(a_bb.index, o_bb.index);
			let abs_diff = (a_bb.coefficient - o_bb.coefficient).abs();
			assert!(abs_diff < 1e-2, "left != right(epsilon={})\n{:?}\n{:?}", 1e-2, a.exp(&m), out);
		}
	}
}