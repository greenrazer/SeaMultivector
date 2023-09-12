use super::*;

// Dual
impl<T:SeaCoefficient, I:SeaBasisBladeSet> SeaLowDimensionalSparseMultivector<T, I> {
	pub fn poincare_dual(&self) -> Self {
		Self { 
			basis_blade_bits: self.basis_blade_bits.dual(), 
			coefficients: self.coefficients.iter().copied().rev().collect()
		}
	}

	pub fn right_complimented(&self) -> Self {
		let mut out = self.poincare_dual();
		out.negate_basis_blades(&I::right_compliment_negations());
		out
	}

	pub fn left_complimented(&self) -> Self {
		let mut out = self.poincare_dual();
		out.negate_basis_blades(&I::left_compliment_negations());
		out
	}
	
	pub fn dual(&self, metric: &SeaMetric) -> Self {
		let mut output = self.poincare_dual();
		output.dualization_helper(metric, false);
		output
	}

	pub fn undual(&self, metric: &SeaMetric) -> Self {
		let mut output = self.poincare_dual();
		output.dualization_helper(metric, true);
		output
	}

}

//Dual Mut
impl<T:SeaCoefficient, I:SeaBasisBladeSet> SeaLowDimensionalSparseMultivector<T, I> {
	pub fn poincare_dualization(&mut self){
		self.basis_blade_bits = self.basis_blade_bits.dual();
		self.coefficients.reverse();
	}

	pub fn right_compliment(&mut self) {
		self.poincare_dualization();
		self.negate_basis_blades(&I::right_compliment_negations());
	}

	pub fn left_compliment(&mut self) {
		self.poincare_dualization();
		self.negate_basis_blades(&I::left_compliment_negations());
	}

	pub fn dualization(&mut self, metric: &SeaMetric) {
		self.poincare_dualization();
		self.dualization_helper(metric, false)
	}

	pub fn undualization(&mut self, metric: &SeaMetric) {
		self.poincare_dualization();
		self.dualization_helper(metric, true)
	}

	fn dualization_helper(&mut self, metric: &SeaMetric, inverse: bool) {
		if metric.has_zero() {
			let negations = if inverse {
				I::left_compliment_negations()
			}
			else {
				I::right_compliment_negations()
			};
			self.negate_basis_blades(&negations);
			return;
		}

		let (_, _, metric_negs) = Self::metric_product_constants(metric);
		let negations = I::dual_negations(
			inverse, 
			metric.pseudoscalar_inverse_negitive() , 
			metric.has_negitive(), 
			&metric_negs
		);

		self.negate_basis_blades(&negations);
	}
}

#[cfg(test)]
mod tests{
	use std::time::Instant;

	use crate::tests::catch_unwind_silent;

	use super::*;

	#[test]
	fn sanity_check() {
		let m = SeaMetric::new(3, 0, 0);
		let a = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("1e13").unwrap();
		let out_dual = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("1e123").unwrap();
		assert_eq!(
			&a ^ &a.dual(&m),
			out_dual
		);

		let m = SeaMetric::new(4, 0, 0);
		let a = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("1e13").unwrap();
		let out_dual = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("-1e1234").unwrap();
		assert_eq!(
			&a ^ &a.dual(&m),
			out_dual
		);

		let m = SeaMetric::new(5, 0, 0);
		let a = SeaLowDimensionalSparseMultivector::<f32, u32>::from_str("1e13").unwrap();
		let out_dual = SeaLowDimensionalSparseMultivector::<f32, u32>::from_str("-1e12345").unwrap();
		assert_eq!(
			&a ^ &a.dual(&m),
			out_dual
		);

		let m = SeaMetric::new(6, 0, 0);
		let a = SeaLowDimensionalSparseMultivector::<f32, u64>::from_str("1e13").unwrap();
		let out_dual = SeaLowDimensionalSparseMultivector::<f32, u64>::from_str("1e123456").unwrap();
		assert_eq!(
			&a ^ &a.dual(&m),
			out_dual
		);

		let m = SeaMetric::new(4, 0, 1);
		let a = SeaLowDimensionalSparseMultivector::<f32, u32>::from_str("1e13").unwrap();
		let out_dual = SeaLowDimensionalSparseMultivector::<f32, u32>::from_str("1e12345").unwrap();
		assert_eq!(
			&a ^ &a.dual(&m),
			out_dual
		);

		let m = SeaMetric::new(4, 0, 1);
		let a = SeaLowDimensionalSparseMultivector::<f32, u32>::from_str("1e124").unwrap();
		let out_dual = SeaLowDimensionalSparseMultivector::<f32, u32>::from_str("1e12345").unwrap();
		assert_eq!(
			&a ^ &a.dual(&m),
			out_dual
		);
	}

	#[test]
	fn positive_metric() {
		let m = SeaMetric::new(3, 0, 0);

		let a = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123").unwrap();
		let out_dual = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("7 + 6e1 - 5e2 - 4e12 + 3e3 + 2e13 - 1e23 - 1e123").unwrap();
		assert_eq!(
			a.dual(&m),
			out_dual
		);
		let out_undual = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("-7 - 6e1 + 5e2 + 4e12 - 3e3 - 2e13 + 1e23 + 1e123").unwrap();
		assert_eq!(
			a.undual(&m),
			out_undual
		);
		assert_eq!(
			a.undual(&m).dual(&m),
			a
		);

		let m = SeaMetric::new(4, 0, 0);

		let a = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e34 + 13e134 + 14e234 + 15e1234").unwrap();
		let out_dual = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("15 +14e1 -13e2 +11e3 -7e4 -12e12 +10e13 -6e14 -9e23 +5e24 -3e34 -8e123 +4e124 -2e134 +e234 +e1234").unwrap();
		assert_eq!(
			a.dual(&m),
			out_dual
		);
		let out_undual = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("15 +14e1 -13e2 +11e3 -7e4 -12e12 +10e13 -6e14 -9e23 +5e24 -3e34 -8e123 +4e124 -2e134 +e234 +e1234").unwrap();
		assert_eq!(
			a.undual(&m),
			out_undual
		);
		assert_eq!(
			a.dual(&m).undual(&m),
			a
		);

		let m = SeaMetric::new(5, 0, 0);

		let a = SeaLowDimensionalSparseMultivector::<f32, u32>::from_str("1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e34 + 13e134 + 14e234 + 15e1234 + e5 + e15 + 2e25 + 3e125 + 4e35 + 5e135 + 6e235 + 7e1235 + 8e45 + 9e145 + 10e245 + 11e1245 + 12e345 + 13e1345 + 14e2345 + 15e12345").unwrap();
		let out_dual = SeaLowDimensionalSparseMultivector::<f32, u32>::from_str("15 +14e1 -13e2 +11e3 -7e4 +15e5 -12e12 +10e13 -6e14 +14e15 -9e23 +5e24 -13e25 -3e34 +11e35 -7e45 -8e123 +4e124 -12e125 -2e134 +10e135 -6e145 +e234 -9e235 +5e245 -3e345 +e1234 -8e1235 +4e1245 -2e1345 +e2345 +e12345").unwrap();
		assert_eq!(
			a.dual(&m),
			out_dual
		);
		let out_undual = SeaLowDimensionalSparseMultivector::<f32, u32>::from_str("15 +14e1 -13e2 +11e3 -7e4 +15e5 -12e12 +10e13 -6e14 +14e15 -9e23 +5e24 -13e25 -3e34 +11e35 -7e45 -8e123 +4e124 -12e125 -2e134 +10e135 -6e145 +e234 -9e235 +5e245 -3e345 +e1234 -8e1235 +4e1245 -2e1345 +e2345 +e12345").unwrap();
		assert_eq!(
			a.undual(&m),
			out_undual
		);
		assert_eq!(
			a.dual(&m).undual(&m),
			a
		);

		let m = SeaMetric::new(6, 0, 0);

		let a = SeaLowDimensionalSparseMultivector::<f32, u64>::from_str("1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e34 + 13e134 + 14e234 + 15e1234 + e5 + e15 + 2e25 + 3e125 + 4e35 + 5e135 + 6e235 + 7e1235 + 8e45 + 9e145 + 10e245 + 11e1245 + 12e345 + 13e1345 + 14e2345 + 15e12345 + e6 + e16 + 2e26 + 3e126 + 4e36 + 5e136 + 6e236 + 7e1236 + 8e46 + 9e146 + 10e246 + 11e1246 + 12e346 + 13e1346 + 14e2346 + 15e12346 + e56 + e156 + 2e256 + 3e1256 + 4e356 + 5e1356 + 6e2356 + 7e12356 + 8e456 + 9e1456 + 10e2456 + 11e12456 + 12e3456 + 13e13456 + 14e23456 + 15e123456").unwrap();
		let out_dual = SeaLowDimensionalSparseMultivector::<f32, u64>::from_str("15 +14e1 -13e2 +11e3 -7e4 +15e5 -15e6 -12e12 +10e13 -6e14 +14e15 -14e16 -9e23 +5e24 -13e25 +13e26 -3e34 +11e35 -11e36 -7e45 +7e46 -15e56 -8e123 +4e124 -12e125 +12e126 -2e134 +10e135 -10e136 -6e145 +6e146 -14e156 +e234 -9e235 +9e236 +5e245 -5e246 +13e256 -3e345 +3e346 -11e356 +7e456 +e1234 -8e1235 +8e1236 +4e1245 -4e1246 +12e1256 -2e1345 +2e1346 -10e1356 +6e1456 +e2345 -e2346 +9e2356 -5e2456 +3e3456 +e12345 -e12346 +8e12356 -4e12456 +2e13456 -e23456 -e123456").unwrap();
		assert_eq!(
			a.dual(&m),
			out_dual
		);
		let out_undual = SeaLowDimensionalSparseMultivector::<f32, u64>::from_str("-15 -14e1 +13e2 -11e3 +7e4 -15e5 +15e6 +12e12 -10e13 +6e14 -14e15 +14e16 +9e23 -5e24 +13e25 -13e26 +3e34 -11e35 +11e36 +7e45 -7e46 +15e56 +8e123 -4e124 +12e125 -12e126 +2e134 -10e135 +10e136 +6e145 -6e146 +14e156 -e234 +9e235 -9e236 -5e245 +5e246 -13e256 +3e345 -3e346 +11e356 -7e456 -e1234 +8e1235 -8e1236 -4e1245 +4e1246 -12e1256 +2e1345 -2e1346 +10e1356 -6e1456 -e2345 +e2346 -9e2356 +5e2456 -3e3456 -e12345 +e12346 -8e12356 +4e12456 -2e13456 +e23456 +e123456").unwrap();
		assert_eq!(
			a.undual(&m),
			out_undual
		);
		assert_eq!(
			a.dual(&m).undual(&m),
			a
		);
	}
	
	#[test]
	fn negitive_metric() {
		let m = SeaMetric::new(3, 1, 0);
		let a = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e34 + 13e134 + 14e234 + 15e1234").unwrap();
		let out_dual = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("15 +14e1 -13e2 +11e3 +7e4 -12e12 +10e13 +6e14 -9e23 -5e24 +3e34 -8e123 -4e124 +2e134 -e234 -e1234").unwrap();
		assert_eq!(
			a.dual(&m),
			out_dual
		);
		let out_undual = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str(" -15 -14e1 +13e2 -11e3 -7e4 +12e12 -10e13 -6e14 +9e23 +5e24 -3e34 +8e123 +4e124 -2e134 +e234 +e1234").unwrap();
		assert_eq!(
			a.undual(&m),
			out_undual
		);
		assert_eq!(
			a.undual(&m).dual(&m),
			a
		);

		let m = SeaMetric::new(2, 2, 0);
		let a = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e34 + 13e134 + 14e234 + 15e1234").unwrap();
		let out_dual = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("15 +14e1 -13e2 -11e3 +7e4 -12e12 -10e13 +6e14 +9e23 -5e24 -3e34 +8e123 -4e124 -2e134 +e234 +e1234").unwrap();
		//846 ns
		assert_eq!(
			a.dual(&m),
			out_dual
		);
		let out_undual = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("15 +14e1 -13e2 -11e3 +7e4 -12e12 -10e13 +6e14 +9e23 -5e24 -3e34 +8e123 -4e124 -2e134 +e234 +e1234").unwrap();
		//835 ns
		assert_eq!(
			a.undual(&m),
			out_undual
		);
		assert_eq!(
			a.dual(&m).undual(&m),
			a
		);

		let m = SeaMetric::new(2, 3, 0);

		let a = SeaLowDimensionalSparseMultivector::<f32, u32>::from_str("1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e34 + 13e134 + 14e234 + 15e1234 + e5 + e15 + 2e25 + 3e125 + 4e35 + 5e135 + 6e235 + 7e1235 + 8e45 + 9e145 + 10e245 + 11e1245 + 12e345 + 13e1345 + 14e2345 + 15e12345").unwrap();
		let out_dual = SeaLowDimensionalSparseMultivector::<f32, u32>::from_str("15 +14e1 -13e2 -11e3 +7e4 -15e5 -12e12 -10e13 +6e14 -14e15 +9e23 -5e24 +13e25 -3e34 +11e35 -7e45 +8e123 -4e124 +12e125 -2e134 +10e135 -6e145 +e234 -9e235 +5e245 +3e345 +e1234 -8e1235 +4e1245 +2e1345 -e2345 -e12345").unwrap();
		assert_eq!(
			a.dual(&m),
			out_dual
		);
		let out_undual = SeaLowDimensionalSparseMultivector::<f32, u32>::from_str(" -15 -14e1 +13e2 +11e3 -7e4 +15e5 +12e12 +10e13 -6e14 +14e15 -9e23 +5e24 -13e25 +3e34 -11e35 +7e45 -8e123 +4e124 -12e125 +2e134 -10e135 +6e145 -e234 +9e235 -5e245 -3e345 -e1234 +8e1235 -4e1245 -2e1345 +e2345 +e12345").unwrap();
		assert_eq!(
			a.undual(&m),
			out_undual
		);
		assert_eq!(
			a.dual(&m).undual(&m),
			a
		);
	}

	#[test]
	fn zero_metric() {
		let m = SeaMetric::new(2, 0, 1);

		let a = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123").unwrap();
		let out = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("7 +3e3 +6e1 -5e2 -2e13 +e23 +4e12 +e123").unwrap();
		assert_eq!(
			a.dual(&m),
			out
		);
		assert_eq!(
			a.dual(&m).undual(&m),
			a
		);

		let m = SeaMetric::new(1, 0, 2);

		let a = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123").unwrap();
		let out = SeaLowDimensionalSparseMultivector::<f32, u8>::from_str("7 -5e2 +3e3 +6e1 +e23 +4e12 -2e13 +e123").unwrap();
		assert_eq!(
			a.dual(&m),
			out
		);
		assert_eq!(
			a.dual(&m).undual(&m),
			a
		);

		let m = SeaMetric::new(3, 0, 1);

		let a = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e34 + 13e134 + 14e234 + 15e1234").unwrap();
		let out = SeaLowDimensionalSparseMultivector::<f32, u16>::from_str("15 + 14e1 - 13e2 + 12e12 + 11e3 - 10e13 + 9e23 + 8e123 - 7e4 + 6e14 - 5e24 - 4e124 + 3e34 + 2e134 - 1e234 + 1e1234").unwrap();
		assert_eq!(
			a.dual(&m),
			out
		);
		assert_eq!(
			a.dual(&m).undual(&m),
			a
		);

		let m = SeaMetric::new(4, 0, 1);

		let a = SeaLowDimensionalSparseMultivector::<f32, u32>::from_str("1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e34 + 13e134 + 14e234 + 15e1234 + e5 + e15 + 2e25 + 3e125 + 4e35 + 5e135 + 6e235 + 7e1235 + 8e45 + 9e145 + 10e245 + 11e1245 + 12e345 + 13e1345 + 14e2345 + 15e12345").unwrap();
		let out = SeaLowDimensionalSparseMultivector::<f32, u32>::from_str("15 +15e5 +14e1 -13e2 +11e3 -7e4 -14e15 +13e25 -11e35 +7e45 +12e12 -10e13 +6e14 +9e23 -5e24 +3e34 +12e125 -10e135 +6e145 +9e235 -5e245 +3e345 +8e123 -4e124 +2e134 -e234 -8e1235 +4e1245 -2e1345 +e2345 +e1234 +e12345").unwrap();
		assert_eq!(
			a.dual(&m),
			out
		);
		assert_eq!(
			a.dual(&m).undual(&m),
			a
		);

		let m = SeaMetric::new(3, 0, 2);

		let a = SeaLowDimensionalSparseMultivector::<f32, u32>::from_str(" + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345").unwrap();
		let out = SeaLowDimensionalSparseMultivector::<f32, u32>::from_str("15 +14e1 -13e2 +7e45 -11e35 +12e12 +6e145 -5e245 +3e345 -10e135 +9e235 +4e1245 -2e1345 +e2345 -8e1235 +e12345").unwrap();
		assert_eq!(
			a.dual(&m),
			out
		);
		assert_eq!(
			a.dual(&m).undual(&m),
			a
		);

		let a = SeaLowDimensionalSparseMultivector::<f32, u32>::from_str("  + 3e13 + 1e4 + 2e124  +4e345").unwrap();
		let out = SeaLowDimensionalSparseMultivector::<f32, u32>::from_str("+4e12 -2e35 -3e245 -e1235").unwrap();
		assert_eq!(
			a.dual(&m),
			out
		);
		assert_eq!(
			a.dual(&m).undual(&m),
			a
		);

		let a = SeaLowDimensionalSparseMultivector::<f32, u32>::from_str(" + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345").unwrap();
		let out = SeaLowDimensionalSparseMultivector::<f32, u32>::from_str("15 +14e1 -13e2 +7e45 -11e35 +12e12 +6e145 -5e245 +3e345 -10e135 +9e235 +4e1245 -2e1345 +e2345 -8e1235 +e12345").unwrap();
		assert_eq!(
			a.dual(&m),
			out
		);
		assert_eq!(
			a.dual(&m).undual(&m),
			a
		);

		let a = SeaLowDimensionalSparseMultivector::<f32, u32>::from_str("  + 3e12 + 1e4 + 2e124  +4e345").unwrap();
		let out = SeaLowDimensionalSparseMultivector::<f32, u32>::from_str("+4e12 -2e35 -e1235 + 3e345").unwrap();
		assert_eq!(
			a.dual(&m),
			out
		);
		assert_eq!(
			a.dual(&m).undual(&m),
			a
		);
	
	}
}
