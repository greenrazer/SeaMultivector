use super::*;

impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> OuterProduct<&Self> for SeaMultivectorLowSparseCoefficents<T, I> {
	type Output = Self;
	fn outer_product(&self, rhs: &Self) -> Self::Output {
		let output_basis_blade_bits = {
			self.basis_blade_bits.outer_product(&rhs.basis_blade_bits)
		};

		let mut output_coefficients = unsafe {
			Self::create_coefficients_array(&output_basis_blade_bits)
		};

		let reverse_rhs_bbbs = rhs.basis_blade_bits.reverse_bits();
		let initial_shift = &SeaBasisBladeBits::<I>::NUM_BASIS_BLADES - &SeaBasisBladeIndexLow::ONE;

		for (out_idx, o) in output_basis_blade_bits.into_iter().enumerate() {
			let cannonical_neg_mask = SeaBasisBladeBits::<I>::cannonical_reordering_neg_mask(&o);

			let mut curr_rhs_bbb = &reverse_rhs_bbbs >> &(&initial_shift - &o);
			curr_rhs_bbb &= &SeaBasisBladeBits::<I>::zero_mask_for_output_index(&o);

			let dot_mask = &self.basis_blade_bits & &curr_rhs_bbb;

			output_coefficients[out_idx] = Self::special_dot_product(
				self,
				rhs,
				&dot_mask,
				&cannonical_neg_mask,
				&o
			);
		}

		Self {
			basis_blade_bits: output_basis_blade_bits,
			coefficients: output_coefficients
		}
	}
}