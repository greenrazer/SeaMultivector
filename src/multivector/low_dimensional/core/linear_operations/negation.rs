use super::*;

impl<T:SeaCoefficient, I:SeaBasisBladeSet> SeaLowDimensionalSparseMultivector<T, I> {
	pub fn negate(&mut self) {
		for i in 0..self.coefficients.len() {
			self.coefficients[i] = -self.coefficients[i]
		}
	}

	pub fn negated(&self) -> Self {
		let mut out = self.clone();
		out.negate();
		out
	}

	pub fn negate_basis_blades(&mut self, blade_bits: &I){
		let negation_bits = self.basis_blade_bits.intersection(blade_bits);
		for idx in negation_bits.iter() {
			let i = self.basis_blade_bits.count_basis_blades_less_than(&idx); 
			self.coefficients[i as usize] = -self.coefficients[i as usize];
		}
	}

	pub fn negated_basis_blades(&self, blade_bits: &I) -> Self{
		let mut new_coeffs = self.coefficients.clone();
		let negation_bits = self.basis_blade_bits.intersection(blade_bits);
		for idx in negation_bits.iter() {
			let i = self.basis_blade_bits.count_basis_blades_less_than(&idx); 
			new_coeffs[i as usize] = -new_coeffs[i as usize];
		}
	
		Self { 
			basis_blade_bits: self.basis_blade_bits.clone(), 
			coefficients: new_coeffs
		}
	}
}