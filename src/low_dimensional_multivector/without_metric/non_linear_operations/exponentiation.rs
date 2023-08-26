use super::*;

//Exp
impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> SeaMultivectorLowSparseCoefficents<T, I> {
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