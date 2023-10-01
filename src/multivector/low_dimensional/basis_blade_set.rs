use std::ops::{Deref, DerefMut, AddAssign, BitOrAssign, Shl, Add, BitOr, BitAndAssign, BitAnd, Not, ShlAssign, ShrAssign, BitXor, Shr, BitXorAssign};
use std::str::FromStr;

use funty::{Unsigned, Integral};

use crate::errors::ParseMultivectorError;

use crate::multivector::grade::SeaGrade;
use crate::multivector::multivector_traits::{GeometricProduct, LeftContraction, RightContraction, InnerProduct, OuterProduct};
use crate::multivector::basis_blade_index::SeaBasisBladeIndexLow;

impl From<&SeaBasisBladeIndexLow> for u8 {
	fn from(index: &SeaBasisBladeIndexLow) -> Self {
		let a = **index;
		Self::ONE << a
	}
}
impl From<SeaBasisBladeIndexLow> for u8 {
	fn from(index: SeaBasisBladeIndexLow) -> Self {
		Self::from(&index)
	}
}
impl From<&SeaBasisBladeIndexLow> for u16 {
	fn from(index: &SeaBasisBladeIndexLow) -> Self {
		let a = **index;
		Self::ONE << a
	}
}
impl From<SeaBasisBladeIndexLow> for u16 {
	fn from(index: SeaBasisBladeIndexLow) -> Self {
		Self::from(&index)
	}
}
impl From<&SeaBasisBladeIndexLow> for u32 {
	fn from(index: &SeaBasisBladeIndexLow) -> Self {
		let a = **index;
		Self::ONE << a
	}
}
impl From<SeaBasisBladeIndexLow> for u32 {
	fn from(index: SeaBasisBladeIndexLow) -> Self {
		Self::from(&index)
	}
}
impl From<&SeaBasisBladeIndexLow> for u64 {
	fn from(index: &SeaBasisBladeIndexLow) -> Self {
		let a = **index;
		Self::ONE << a
	}
}
impl From<SeaBasisBladeIndexLow> for u64 {
	fn from(index: SeaBasisBladeIndexLow) -> Self {
		Self::from(&index)
	}
}


pub trait SeaBasisBladeUnsignedConstants: Unsigned + From<SeaBasisBladeIndexLow> + for<'a> From<&'a SeaBasisBladeIndexLow>{
	const DIMENSIONS: u8;
	const NUM_BASIS_BLADES: usize = Self::BITS as usize;

	// const ONE: Self = Self::ONE;
	
	/// All bits of an even grade,
	/// or all bits with an index with an even number of bits set
	/// or even binary weight.
	/// 
	/// e.g 10010110
	/// because for each index
	/// 0 -> 0 bits set -> even
	/// 1 -> 1 bits set -> odd
	/// 2 -> 1 bits set -> odd
	/// 3 -> 2 bits set -> even
	/// 4 -> 1 bits set -> odd
	/// 5 -> 2 bits set -> odd
	/// 6 -> 2 bits set -> odd
	/// 7 -> 3 bits set -> odd
	const EVEN_GRADE_BITS: Self;
	
	// 1001 0110  0110 1001  0110 1001  1001 0110  0110 1001  1001 0110  1001 0110  0110 1001
	const EVEN_GRADE_UNSIGNED_BITS: u64 = 0x96_69_69_96_69_96_96_69;

	/// All bits of an odd grade,
	/// or all bits with an index with an odd number of bits set
	/// or odd binary weight.
	const ODD_GRADE_BITS: Self;

	// 0110 1001  1001 0110  1001 0110  0110 1001  1001 0110  0110 1001  0110 1001  1001 0110
	const ODD_GRADE_UNSIGNED_BITS: u64 = !Self::EVEN_GRADE_UNSIGNED_BITS;

	/// All indexes that would be negitive after
	/// being cannonically reordered with themselves,
	/// or all indexes that would be negitive after reversing all
	/// it's basis vectors.
	/// 
	/// if 1/2 * binary weight of index * (binary weight of index-1) is even -> 0
	/// if odd -> 1
	const REVERSION_NEGATION_BITS: Self;

	// 1000 0001  0001 0111  0001 0111  0111 1110  0001 0111  0111 1110  0111 1110  1110 1000
	const REVERSION_NEGATION_UNSIGNED_BITS: u64 = 0x81_17_17_7E_17_7E_7E_E8;

	const INVOLUTION_NEGATION_BITS: Self;

	const INVOLUTION_NEGATION_UNSIGNED_BITS: u64 = Self::ODD_GRADE_UNSIGNED_BITS;

	/// xor between reversion_negation_bits and odd_grade_bits
	const CLIFFORD_NEGATION_BITS: Self;

	// 1110 1000  1000 0001  1000 0001  0001 0111  1000 0001  0001 0111  0001 0111  0111 1110
	const CLIFFORD_NEGATION_UNSIGNED_BITS: u64 = 0xE8_81_81_17_81_17_17_7E;

	/// right compliment negation bits
	/// such that x * compliment(x) = I(pseudo scalar)
	/// complement = reverse(x)*these_bits
	const RIGHT_COMPLIMENT_NEGATION_BITS: Self;

	/// left compliment negation bits
	/// such that compliment(x) * x  = I(pseudo scalar)
	/// complement = reverse(x)*these_bits
	const LEFT_COMPLIMENT_NEGATION_BITS: Self;

	/// negation bits corresponding to the sign of
	/// sign(x << I(pseudo scalar))
	/// for positive-only metrics
	const DUAL_NEGATION_BITS: Self;
	
	/// for each grade, the bits in that grade
	/// for example for grade 0 the only bit set is always in the 0's place
	/// for grade 1 every index with a binary weight of 1 is 1
	/// for grade 2 every index with a binary weight of 2 is 1
	/// ...
	const GRADE_BITS: [Self; 7];

	const GRADE_UNSIGNED_BITS: [u64; 7] = [
		0x00_00_00_00_00_00_00_01,
		0x00_00_00_01_00_01_01_16,
		0x00_01_01_16_01_16_16_68,
		0x01_16_16_68_16_68_68_80,
		0x16_68_68_80_68_80_80_00,
		0x68_80_80_00_80_00_00_00,
		0x80_00_00_00_00_00_00_00
	];

	
	/// for each basis_vector(1->0, e1->1, e2->2, e3->3)
	/// 1 at each index that contains that basis vector
	/// e.g. 
	/// basis vector 1
	/// 01010101
	/// basis vector 2
	/// 00110011
	const CONTAINS_BASIS_VEC_BITS: [Self; 7];

	const CONTAINS_BASIS_VEC_UNSIGNED_BITS: [u64; 7] = [
		0xFF_FF_FF_FF_FF_FF_FF_FF,
		0xAA_AA_AA_AA_AA_AA_AA_AA,
		0xCC_CC_CC_CC_CC_CC_CC_CC,
		0xF0_F0_F0_F0_F0_F0_F0_F0,
		0xFF_00_FF_00_FF_00_FF_00,
		0xFF_FF_00_00_FF_FF_00_00,
		0xFF_FF_FF_FF_00_00_00_00
	];
	
	/// opposite of contains basis vector bits
	const EXCLUDES_BASIS_VEC_BITS: [Self; 7];

	const EXCLUDES_BASIS_VEC_UNSIGNED_BITS: [u64; 7] = [
		!Self::CONTAINS_BASIS_VEC_UNSIGNED_BITS[0],
		!Self::CONTAINS_BASIS_VEC_UNSIGNED_BITS[1],
		!Self::CONTAINS_BASIS_VEC_UNSIGNED_BITS[2],
		!Self::CONTAINS_BASIS_VEC_UNSIGNED_BITS[3],
		!Self::CONTAINS_BASIS_VEC_UNSIGNED_BITS[4],
		!Self::CONTAINS_BASIS_VEC_UNSIGNED_BITS[5],
		!Self::CONTAINS_BASIS_VEC_UNSIGNED_BITS[6],
	];

	/// for each basis vector in an output index after some product
	/// this negates the reversion grade bits to get the current
	/// cannonical reordering sign bits
	/// for example for each product output index
	/// the cannonical reordering sign mask is
	/// 1 -> Reversion bits
	/// e1 -> Reversion bits ^ output reversion negation bits[1]
	/// e2 -> Reversion bits ^ output reversion negation bits[2]
	/// e12 -> Reversion bits ^ output reversion negation bits[1] ^ output reversion negation bits[2]
	const OUTPUT_CANNONICAL_REORDERING_SIGN_REVERSION_NEGATION_BITS: [Self; 7];

	const K: [u64; 7] = Self::CONTAINS_BASIS_VEC_UNSIGNED_BITS;
	const OUTPUT_CANNONICAL_REORDERING_SIGN_REVERSION_NEGATION_UNSIGNED_BITS: [u64; 7] = [
		Self::K[1] ^ Self::K[2] ^ Self::K[3] ^ Self::K[4] ^ Self::K[5] ^ Self::K[6],
				Self::K[2] ^ Self::K[3] ^ Self::K[4] ^ Self::K[5] ^ Self::K[6],
							Self::K[3] ^ Self::K[4] ^ Self::K[5] ^ Self::K[6],		
										Self::K[4] ^ Self::K[5] ^ Self::K[6],
													Self::K[5] ^ Self::K[6],
															Self::K[6],
																	0,
	];
}
impl SeaBasisBladeUnsignedConstants for u8 {
	const DIMENSIONS: u8 = 3;

	const EVEN_GRADE_BITS: Self = Self::EVEN_GRADE_UNSIGNED_BITS as Self;
	const ODD_GRADE_BITS: Self = Self::ODD_GRADE_UNSIGNED_BITS as Self;
	const REVERSION_NEGATION_BITS: Self = Self::REVERSION_NEGATION_UNSIGNED_BITS as Self;
	const INVOLUTION_NEGATION_BITS: Self = Self::INVOLUTION_NEGATION_UNSIGNED_BITS as Self;
	const CLIFFORD_NEGATION_BITS: Self = Self::CLIFFORD_NEGATION_UNSIGNED_BITS as Self;
	// 0010 0100
	const RIGHT_COMPLIMENT_NEGATION_BITS: Self = 0x24;
	// 0010 0100
	const LEFT_COMPLIMENT_NEGATION_BITS: Self = 0x24;
	// 0011 0011
	const DUAL_NEGATION_BITS: Self = 0x33;

	const GRADE_BITS: [Self;7] = [
		Self::GRADE_UNSIGNED_BITS[0] as Self,
		Self::GRADE_UNSIGNED_BITS[1] as Self,
		Self::GRADE_UNSIGNED_BITS[2] as Self,
		Self::GRADE_UNSIGNED_BITS[3] as Self,
		Self::GRADE_UNSIGNED_BITS[4] as Self,
		Self::GRADE_UNSIGNED_BITS[5] as Self,
		Self::GRADE_UNSIGNED_BITS[6] as Self,
	];
	const CONTAINS_BASIS_VEC_BITS: [Self; 7] = [
		Self::CONTAINS_BASIS_VEC_UNSIGNED_BITS[0] as Self,
		Self::CONTAINS_BASIS_VEC_UNSIGNED_BITS[1] as Self,
		Self::CONTAINS_BASIS_VEC_UNSIGNED_BITS[2] as Self,
		Self::CONTAINS_BASIS_VEC_UNSIGNED_BITS[3] as Self,
		Self::CONTAINS_BASIS_VEC_UNSIGNED_BITS[4] as Self,
		Self::CONTAINS_BASIS_VEC_UNSIGNED_BITS[5] as Self,
		Self::CONTAINS_BASIS_VEC_UNSIGNED_BITS[6] as Self,
	];
	const EXCLUDES_BASIS_VEC_BITS: [Self; 7] = [
		Self::EXCLUDES_BASIS_VEC_UNSIGNED_BITS[0] as Self,
		Self::EXCLUDES_BASIS_VEC_UNSIGNED_BITS[1] as Self,
		Self::EXCLUDES_BASIS_VEC_UNSIGNED_BITS[2] as Self,
		Self::EXCLUDES_BASIS_VEC_UNSIGNED_BITS[3] as Self,
		Self::EXCLUDES_BASIS_VEC_UNSIGNED_BITS[4] as Self,
		Self::EXCLUDES_BASIS_VEC_UNSIGNED_BITS[5] as Self,
		Self::EXCLUDES_BASIS_VEC_UNSIGNED_BITS[6] as Self,
	];
	const OUTPUT_CANNONICAL_REORDERING_SIGN_REVERSION_NEGATION_BITS: [Self; 7] = [
		Self::OUTPUT_CANNONICAL_REORDERING_SIGN_REVERSION_NEGATION_UNSIGNED_BITS[0] as Self,
		Self::OUTPUT_CANNONICAL_REORDERING_SIGN_REVERSION_NEGATION_UNSIGNED_BITS[1] as Self,
		Self::OUTPUT_CANNONICAL_REORDERING_SIGN_REVERSION_NEGATION_UNSIGNED_BITS[2] as Self,
		Self::OUTPUT_CANNONICAL_REORDERING_SIGN_REVERSION_NEGATION_UNSIGNED_BITS[3] as Self,
		Self::OUTPUT_CANNONICAL_REORDERING_SIGN_REVERSION_NEGATION_UNSIGNED_BITS[4] as Self,
		Self::OUTPUT_CANNONICAL_REORDERING_SIGN_REVERSION_NEGATION_UNSIGNED_BITS[5] as Self,
		Self::OUTPUT_CANNONICAL_REORDERING_SIGN_REVERSION_NEGATION_UNSIGNED_BITS[6] as Self,
	];
}
impl SeaBasisBladeUnsignedConstants for u16 {
	const DIMENSIONS: u8 = 4;
	const EVEN_GRADE_BITS: Self = Self::EVEN_GRADE_UNSIGNED_BITS as Self;
	const ODD_GRADE_BITS: Self = Self::ODD_GRADE_UNSIGNED_BITS as Self;
	const REVERSION_NEGATION_BITS: Self = Self::REVERSION_NEGATION_UNSIGNED_BITS as Self;
	const INVOLUTION_NEGATION_BITS: Self = Self::INVOLUTION_NEGATION_UNSIGNED_BITS as Self;
	const CLIFFORD_NEGATION_BITS: Self = Self::CLIFFORD_NEGATION_UNSIGNED_BITS as Self;
	// 0100 1101 0010 0100
	const RIGHT_COMPLIMENT_NEGATION_BITS: Self = 0x4D_24;
	// 0010 0100 1011 0010
	const LEFT_COMPLIMENT_NEGATION_BITS: Self = 0x24_B2;
	// 0011 0011 1100 1100
	const DUAL_NEGATION_BITS: Self = 0x33_CC;

	const GRADE_BITS: [Self;7] = [
		Self::GRADE_UNSIGNED_BITS[0] as Self,
		Self::GRADE_UNSIGNED_BITS[1] as Self,
		Self::GRADE_UNSIGNED_BITS[2] as Self,
		Self::GRADE_UNSIGNED_BITS[3] as Self,
		Self::GRADE_UNSIGNED_BITS[4] as Self,
		Self::GRADE_UNSIGNED_BITS[5] as Self,
		Self::GRADE_UNSIGNED_BITS[6] as Self,
	];
	const CONTAINS_BASIS_VEC_BITS: [Self; 7] = [
		Self::CONTAINS_BASIS_VEC_UNSIGNED_BITS[0] as Self,
		Self::CONTAINS_BASIS_VEC_UNSIGNED_BITS[1] as Self,
		Self::CONTAINS_BASIS_VEC_UNSIGNED_BITS[2] as Self,
		Self::CONTAINS_BASIS_VEC_UNSIGNED_BITS[3] as Self,
		Self::CONTAINS_BASIS_VEC_UNSIGNED_BITS[4] as Self,
		Self::CONTAINS_BASIS_VEC_UNSIGNED_BITS[5] as Self,
		Self::CONTAINS_BASIS_VEC_UNSIGNED_BITS[6] as Self,
	];
	const EXCLUDES_BASIS_VEC_BITS: [Self; 7] = [
		Self::EXCLUDES_BASIS_VEC_UNSIGNED_BITS[0] as Self,
		Self::EXCLUDES_BASIS_VEC_UNSIGNED_BITS[1] as Self,
		Self::EXCLUDES_BASIS_VEC_UNSIGNED_BITS[2] as Self,
		Self::EXCLUDES_BASIS_VEC_UNSIGNED_BITS[3] as Self,
		Self::EXCLUDES_BASIS_VEC_UNSIGNED_BITS[4] as Self,
		Self::EXCLUDES_BASIS_VEC_UNSIGNED_BITS[5] as Self,
		Self::EXCLUDES_BASIS_VEC_UNSIGNED_BITS[6] as Self,
	];
	const OUTPUT_CANNONICAL_REORDERING_SIGN_REVERSION_NEGATION_BITS: [Self; 7] = [
		Self::OUTPUT_CANNONICAL_REORDERING_SIGN_REVERSION_NEGATION_UNSIGNED_BITS[0] as Self,
		Self::OUTPUT_CANNONICAL_REORDERING_SIGN_REVERSION_NEGATION_UNSIGNED_BITS[1] as Self,
		Self::OUTPUT_CANNONICAL_REORDERING_SIGN_REVERSION_NEGATION_UNSIGNED_BITS[2] as Self,
		Self::OUTPUT_CANNONICAL_REORDERING_SIGN_REVERSION_NEGATION_UNSIGNED_BITS[3] as Self,
		Self::OUTPUT_CANNONICAL_REORDERING_SIGN_REVERSION_NEGATION_UNSIGNED_BITS[4] as Self,
		Self::OUTPUT_CANNONICAL_REORDERING_SIGN_REVERSION_NEGATION_UNSIGNED_BITS[5] as Self,
		Self::OUTPUT_CANNONICAL_REORDERING_SIGN_REVERSION_NEGATION_UNSIGNED_BITS[6] as Self,
	];
}
impl SeaBasisBladeUnsignedConstants for u32 {
	const DIMENSIONS: u8 = 5;
	const EVEN_GRADE_BITS: Self = Self::EVEN_GRADE_UNSIGNED_BITS as Self;
	const ODD_GRADE_BITS: Self = Self::ODD_GRADE_UNSIGNED_BITS as Self;
	const REVERSION_NEGATION_BITS: Self = Self::REVERSION_NEGATION_UNSIGNED_BITS as Self;
	const INVOLUTION_NEGATION_BITS: Self = Self::INVOLUTION_NEGATION_UNSIGNED_BITS as Self;
	const CLIFFORD_NEGATION_BITS: Self = Self::CLIFFORD_NEGATION_UNSIGNED_BITS as Self;
	// 0010 0100 1011 0010 0100 1101 0010 0100
	const RIGHT_COMPLIMENT_NEGATION_BITS: Self = 0x24_B2_4D_24;
	// 0010 0100 1011 0010 0100 1101 0010 0100
	const LEFT_COMPLIMENT_NEGATION_BITS: Self = 0x24_B2_4D_24;
	// 0011 0011 1100 1100 0011 0011 1100 1100
	const DUAL_NEGATION_BITS: Self = 0x33_CC_33_CC;

	const GRADE_BITS: [Self;7] = [
		Self::GRADE_UNSIGNED_BITS[0] as Self,
		Self::GRADE_UNSIGNED_BITS[1] as Self,
		Self::GRADE_UNSIGNED_BITS[2] as Self,
		Self::GRADE_UNSIGNED_BITS[3] as Self,
		Self::GRADE_UNSIGNED_BITS[4] as Self,
		Self::GRADE_UNSIGNED_BITS[5] as Self,
		Self::GRADE_UNSIGNED_BITS[6] as Self,
	];
	const CONTAINS_BASIS_VEC_BITS: [Self; 7] = [
		Self::CONTAINS_BASIS_VEC_UNSIGNED_BITS[0] as Self,
		Self::CONTAINS_BASIS_VEC_UNSIGNED_BITS[1] as Self,
		Self::CONTAINS_BASIS_VEC_UNSIGNED_BITS[2] as Self,
		Self::CONTAINS_BASIS_VEC_UNSIGNED_BITS[3] as Self,
		Self::CONTAINS_BASIS_VEC_UNSIGNED_BITS[4] as Self,
		Self::CONTAINS_BASIS_VEC_UNSIGNED_BITS[5] as Self,
		Self::CONTAINS_BASIS_VEC_UNSIGNED_BITS[6] as Self,
	];
	const EXCLUDES_BASIS_VEC_BITS: [Self; 7] = [
		Self::EXCLUDES_BASIS_VEC_UNSIGNED_BITS[0] as Self,
		Self::EXCLUDES_BASIS_VEC_UNSIGNED_BITS[1] as Self,
		Self::EXCLUDES_BASIS_VEC_UNSIGNED_BITS[2] as Self,
		Self::EXCLUDES_BASIS_VEC_UNSIGNED_BITS[3] as Self,
		Self::EXCLUDES_BASIS_VEC_UNSIGNED_BITS[4] as Self,
		Self::EXCLUDES_BASIS_VEC_UNSIGNED_BITS[5] as Self,
		Self::EXCLUDES_BASIS_VEC_UNSIGNED_BITS[6] as Self,
	];
	const OUTPUT_CANNONICAL_REORDERING_SIGN_REVERSION_NEGATION_BITS: [Self; 7] = [
		Self::OUTPUT_CANNONICAL_REORDERING_SIGN_REVERSION_NEGATION_UNSIGNED_BITS[0] as Self,
		Self::OUTPUT_CANNONICAL_REORDERING_SIGN_REVERSION_NEGATION_UNSIGNED_BITS[1] as Self,
		Self::OUTPUT_CANNONICAL_REORDERING_SIGN_REVERSION_NEGATION_UNSIGNED_BITS[2] as Self,
		Self::OUTPUT_CANNONICAL_REORDERING_SIGN_REVERSION_NEGATION_UNSIGNED_BITS[3] as Self,
		Self::OUTPUT_CANNONICAL_REORDERING_SIGN_REVERSION_NEGATION_UNSIGNED_BITS[4] as Self,
		Self::OUTPUT_CANNONICAL_REORDERING_SIGN_REVERSION_NEGATION_UNSIGNED_BITS[5] as Self,
		Self::OUTPUT_CANNONICAL_REORDERING_SIGN_REVERSION_NEGATION_UNSIGNED_BITS[6] as Self,
	];
}
impl SeaBasisBladeUnsignedConstants for u64 {
	const DIMENSIONS: u8 = 6;
	const EVEN_GRADE_BITS: Self = Self::EVEN_GRADE_UNSIGNED_BITS as Self;
	const ODD_GRADE_BITS: Self = Self::ODD_GRADE_UNSIGNED_BITS as Self;
	const REVERSION_NEGATION_BITS: Self = Self::REVERSION_NEGATION_UNSIGNED_BITS as Self;
	const INVOLUTION_NEGATION_BITS: Self = Self::INVOLUTION_NEGATION_UNSIGNED_BITS as Self;
	const CLIFFORD_NEGATION_BITS: Self = Self::CLIFFORD_NEGATION_UNSIGNED_BITS as Self;
	// 0100 1011 0010 0100 1101 1011 0100 1101 0010 0100 1011 0010 0100 1101 0010 0100
	const RIGHT_COMPLIMENT_NEGATION_BITS: Self = 0x4D_24_DB_4D_24_B2_4D_24;
	// 0010 0100 1011 0010 0100 1101 0010 0100 1011 0010 1101 1011 0010 0100 1101 0010
	const LEFT_COMPLIMENT_NEGATION_BITS: Self = 0x24_D2_4D_24_B2_DB_24_D2;
	// 0011 0011 1100 1100 0011 0011 1100 1100 1100 1100 0011 0011 1100 1100 0011 0011
	const DUAL_NEGATION_BITS: Self = 0x33_CC_33_CC_CC_33_CC_33;

	const GRADE_BITS: [Self;7] = [
		Self::GRADE_UNSIGNED_BITS[0] as Self,
		Self::GRADE_UNSIGNED_BITS[1] as Self,
		Self::GRADE_UNSIGNED_BITS[2] as Self,
		Self::GRADE_UNSIGNED_BITS[3] as Self,
		Self::GRADE_UNSIGNED_BITS[4] as Self,
		Self::GRADE_UNSIGNED_BITS[5] as Self,
		Self::GRADE_UNSIGNED_BITS[6] as Self,
	];
	const CONTAINS_BASIS_VEC_BITS: [Self; 7] = [
		Self::CONTAINS_BASIS_VEC_UNSIGNED_BITS[0] as Self,
		Self::CONTAINS_BASIS_VEC_UNSIGNED_BITS[1] as Self,
		Self::CONTAINS_BASIS_VEC_UNSIGNED_BITS[2] as Self,
		Self::CONTAINS_BASIS_VEC_UNSIGNED_BITS[3] as Self,
		Self::CONTAINS_BASIS_VEC_UNSIGNED_BITS[4] as Self,
		Self::CONTAINS_BASIS_VEC_UNSIGNED_BITS[5] as Self,
		Self::CONTAINS_BASIS_VEC_UNSIGNED_BITS[6] as Self,
	];
	const EXCLUDES_BASIS_VEC_BITS: [Self; 7] = [
		Self::EXCLUDES_BASIS_VEC_UNSIGNED_BITS[0] as Self,
		Self::EXCLUDES_BASIS_VEC_UNSIGNED_BITS[1] as Self,
		Self::EXCLUDES_BASIS_VEC_UNSIGNED_BITS[2] as Self,
		Self::EXCLUDES_BASIS_VEC_UNSIGNED_BITS[3] as Self,
		Self::EXCLUDES_BASIS_VEC_UNSIGNED_BITS[4] as Self,
		Self::EXCLUDES_BASIS_VEC_UNSIGNED_BITS[5] as Self,
		Self::EXCLUDES_BASIS_VEC_UNSIGNED_BITS[6] as Self,
	];
	const OUTPUT_CANNONICAL_REORDERING_SIGN_REVERSION_NEGATION_BITS: [Self; 7] = [
		Self::OUTPUT_CANNONICAL_REORDERING_SIGN_REVERSION_NEGATION_UNSIGNED_BITS[0] as Self,
		Self::OUTPUT_CANNONICAL_REORDERING_SIGN_REVERSION_NEGATION_UNSIGNED_BITS[1] as Self,
		Self::OUTPUT_CANNONICAL_REORDERING_SIGN_REVERSION_NEGATION_UNSIGNED_BITS[2] as Self,
		Self::OUTPUT_CANNONICAL_REORDERING_SIGN_REVERSION_NEGATION_UNSIGNED_BITS[3] as Self,
		Self::OUTPUT_CANNONICAL_REORDERING_SIGN_REVERSION_NEGATION_UNSIGNED_BITS[4] as Self,
		Self::OUTPUT_CANNONICAL_REORDERING_SIGN_REVERSION_NEGATION_UNSIGNED_BITS[5] as Self,
		Self::OUTPUT_CANNONICAL_REORDERING_SIGN_REVERSION_NEGATION_UNSIGNED_BITS[6] as Self,
	];
}

pub trait SeaBasisBladeSetUnsigned: SeaBasisBladeUnsignedConstants + From<SeaBasisBladeIndexLow> + for<'a> From<&'a SeaBasisBladeIndexLow>{ 
	/// sets all least signifigant bits up to n
	/// e.g. 5
	/// ...0011111
	/// 3
	/// ...0000111
	fn basis_blades_less_than_mask(idx: &SeaBasisBladeIndexLow) -> Self {
		let one = Self::ONE;
		one.checked_shl(**idx as u32).unwrap_or(Self::ZERO).overflowing_sub(one).0
	}

	fn basis_blades_greater_than_mask(idx: &SeaBasisBladeIndexLow) -> Self {
		!Self::basis_blades_less_than_mask(idx)
	}

	fn outer_product_zero_mask(output_index: &SeaBasisBladeIndexLow) -> Self {
		let mut out = Self::ONE;
		for bv in output_index {
			out ^= out << (1 << bv);
		}
		out
	}

	fn cannonical_reordering_negations_mask(output_index: &SeaBasisBladeIndexLow) -> Self {
		let mut out = Self::REVERSION_NEGATION_BITS;
		for i in output_index {
			let temp = Self::OUTPUT_CANNONICAL_REORDERING_SIGN_REVERSION_NEGATION_BITS[(i+1) as usize];
			out ^= temp;
		}
		out
	}

	fn metric_negations_mask(output_index: &SeaBasisBladeIndexLow, metric_negs: &SeaBasisBladeIndexLow) -> Self {
		let curr_met_negs = &!output_index & metric_negs;
		let mut out = Self::ZERO;
		for i in &curr_met_negs {
			let temp = Self::CONTAINS_BASIS_VEC_BITS[(i+1) as usize];
			out ^= temp;
		}
		out
	}

	fn dual_metric_negations_mask(index: &SeaBasisBladeIndexLow) -> Self {
		let mut out = Self::ZERO;
		for i in index {
			let temp = !Self::CONTAINS_BASIS_VEC_BITS[(i+1) as usize];
			out ^= temp;
		}
		out
	}

	fn remove_basis_blades_with_any_mask(index: &SeaBasisBladeIndexLow) -> Self {
		let mut out = Self::CONTAINS_BASIS_VEC_BITS[0];
		for i in index {
			out &= Self::EXCLUDES_BASIS_VEC_BITS[(i+1) as usize];
		}
		out
	}

	fn remove_basis_blades_without_all_mask(index: &SeaBasisBladeIndexLow) -> Self {
		let mut out = Self::CONTAINS_BASIS_VEC_BITS[0];
		for i in index {
			out &= Self::CONTAINS_BASIS_VEC_BITS[(i+1) as usize];
		}
		out
	}

	/// swaps adjacent bit groups with length 2^n
	/// e.g. 2 1111011001111001
	/// 2^2 = 4
	/// 1111 0110 0111 1001
	/// 0110 1111 1001 0111
	/// 0110111110010111
	fn geometric_product_with_basis_vec(&mut self, n:u8) {
		assert_ne!(n, 0);
		let bitmask = Self::CONTAINS_BASIS_VEC_BITS[n as usize];
		let shift_amnt = Self::ONE << (n-1);
		let msb = *self & bitmask;
		*self &= !bitmask;
		*self <<= shift_amnt;
		*self |= msb >> shift_amnt;
	}
}
impl<I:SeaBasisBladeUnsignedConstants> SeaBasisBladeSetUnsigned for I {}

pub struct SeaBasisBladeSetIter<I:SeaBasisBladeSet>(I);
impl<I: SeaBasisBladeSetUnsigned> Iterator for SeaBasisBladeSetIter<I> {
	type Item = SeaBasisBladeIndexLow;

	fn next(&mut self) -> Option<Self::Item> {
		if self.0 == I::ZERO {
			return None;
		}
		let out = self.0.trailing_zeros();
		self.0 ^= I::ONE << out;
		Some(SeaBasisBladeIndexLow(out as u8))
	}

	fn size_hint(&self) -> (usize, Option<usize>) {
		(0, Some(self.0.count_ones() as usize))
	}

	fn count(self) -> usize
		where
			Self: Sized, {
		self.0.count_ones() as usize
	} 
}

pub trait SeaBasisBladeSet:
	From<SeaBasisBladeIndexLow> + for<'a> From<&'a SeaBasisBladeIndexLow> 
	+ Clone
	+ LeftContraction<(Self, u8), Output=Self>
	+ for<'a> LeftContraction<(&'a Self, u8), Output=Self>
	+ RightContraction<(Self, u8), Output=Self>
	+ for<'a> RightContraction<(&'a Self, u8), Output=Self>
	+ InnerProduct<(Self, u8), Output=Self>
	+ for<'a> InnerProduct<(&'a Self, u8), Output=Self>
	+ OuterProduct<Self, Output=Self>
	+ for<'a> OuterProduct<&'a Self, Output=Self>
	+ GeometricProduct<(Self, u8), Output=Self>
	+ for<'a> GeometricProduct<(&'a Self, u8), Output=Self>
		where
			Self: Sized
{
	const DIMENSIONS: u8;
	const NUM_BASIS_BLADES: usize;
	const MAX_BASIS_BLADE_INDEX: SeaBasisBladeIndexLow;

	const EMPTY: Self;
	const SCALAR: Self;
	const ALL: Self;

	type IntoIterType: Iterator<Item=SeaBasisBladeIndexLow>;

	fn reversion_negations() -> Self;
	fn involution_negations() -> Self;
	fn clifford_negations() -> Self;
	fn right_compliment_negations() -> Self;
	fn left_compliment_negations() -> Self;
	fn dual_negations(inverse: bool, pseudoscalar_inverse_neg: bool, metric_has_neg: bool, metric_negs: &SeaBasisBladeIndexLow) -> Self;
	
	fn cannonical_reordering_negations(output_index: &SeaBasisBladeIndexLow) -> Self;
	fn metric_negations(output_index: &SeaBasisBladeIndexLow, metric_negs: &SeaBasisBladeIndexLow) -> Self;
	fn total_negations(output_index: &SeaBasisBladeIndexLow, metric_negs: &SeaBasisBladeIndexLow) -> Self;
	fn dual_metric_negations(output_index: &SeaBasisBladeIndexLow) -> Self;

	fn outer_product_zeros(output_index: &SeaBasisBladeIndexLow) -> Self;
	fn basis_blades_less_than(index: &SeaBasisBladeIndexLow) -> Self;
	fn basis_blades_greater_than(index: &SeaBasisBladeIndexLow) -> Self;
	fn grade(grade:&SeaGrade) -> Self;
	
	fn contains_basis_vec(&self, basis_vector_index: &u8) -> Self;
	fn excludes_basis_vec(&self, basis_vector_index: &u8) -> Self;
	fn is_empty(&self) -> bool;
	fn is_scalar(&self) -> bool;
	fn has(&self, index: &SeaBasisBladeIndexLow) -> bool;
	fn count_basis_blades(&self) -> u32;
	fn count_basis_blades_less_than(&self, idx:&SeaBasisBladeIndexLow) -> u32;
	fn dual(&self) -> Self;
	fn dual_index(index:&SeaBasisBladeIndexLow) -> SeaBasisBladeIndexLow;

	fn left_contraction_with_basis_blade(&self, index: &SeaBasisBladeIndexLow, metric_non_zeros: &SeaBasisBladeIndexLow) -> Self;
	fn right_contraction_with_basis_blade(&self, index: &SeaBasisBladeIndexLow, metric_non_zeros: &SeaBasisBladeIndexLow) -> Self;
	fn outer_product_with_basis_blade(&self, index: &SeaBasisBladeIndexLow) -> Self;
	fn geometric_product_with_basis_blade(&self, index:&SeaBasisBladeIndexLow, metric_non_zeros: &SeaBasisBladeIndexLow) -> Self;

	fn iter(&self) -> Self::IntoIterType;
	fn intersection(&self, other:&Self) -> Self;
	fn union(&self, other:&Self) -> Self;
	fn symmetric_difference(&self, other:&Self) -> Self;

	fn set(&mut self, index: &SeaBasisBladeIndexLow, value:bool);
	fn toggle(&mut self, index: &SeaBasisBladeIndexLow);
	fn intersection_mut(&mut self, other:&Self);
	fn union_mut(&mut self, other:&Self);
	fn symmetric_difference_mut(&mut self, other:&Self);
}

impl<I: SeaBasisBladeSetUnsigned> SeaBasisBladeSet for I {
	const DIMENSIONS: u8 = Self::DIMENSIONS;
	const NUM_BASIS_BLADES: usize = Self::NUM_BASIS_BLADES;
	const MAX_BASIS_BLADE_INDEX: SeaBasisBladeIndexLow = SeaBasisBladeIndexLow((Self::NUM_BASIS_BLADES - 1) as u8);

	const EMPTY: Self = Self::ZERO;
	const SCALAR: Self = Self::ONE;
	const ALL: Self = Self::CONTAINS_BASIS_VEC_BITS[0];

	type IntoIterType = SeaBasisBladeSetIter<I>;

	fn reversion_negations() -> Self {
		Self::REVERSION_NEGATION_BITS
	}

	fn involution_negations() -> Self {
		Self::INVOLUTION_NEGATION_BITS
	}

	fn clifford_negations() -> Self {
		Self::CLIFFORD_NEGATION_BITS
	}

	fn right_compliment_negations() -> Self {
		Self::RIGHT_COMPLIMENT_NEGATION_BITS
	}

	fn left_compliment_negations() -> Self {
		Self::LEFT_COMPLIMENT_NEGATION_BITS
	}

	fn dual_negations(inverse: bool, pseudoscalar_inverse_neg: bool, metric_has_neg: bool, metric_negs: &SeaBasisBladeIndexLow) -> Self {
		let mut negations = Self::DUAL_NEGATION_BITS;

		if !inverse && pseudoscalar_inverse_neg {
			negations = !negations
		}

		if metric_has_neg {
			let metric_neg_mask = Self::dual_metric_negations_mask(metric_negs);
			negations ^= &metric_neg_mask;
		}

		negations
	}

	fn cannonical_reordering_negations(output_index: &SeaBasisBladeIndexLow) -> Self {
		Self::cannonical_reordering_negations_mask(output_index)
	}

	fn metric_negations(output_index: &SeaBasisBladeIndexLow, metric_negs: &SeaBasisBladeIndexLow) -> Self {
		Self::metric_negations_mask(output_index, metric_negs)
	}
	fn total_negations(output_index: &SeaBasisBladeIndexLow,  metric_negs: &SeaBasisBladeIndexLow) -> Self {
		Self::cannonical_reordering_negations_mask(output_index) ^ Self::metric_negations_mask(output_index, metric_negs)
	}
	fn dual_metric_negations(output_index: &SeaBasisBladeIndexLow) -> Self{
		Self::dual_metric_negations_mask(output_index)
	}


	fn outer_product_zeros(output_index: &SeaBasisBladeIndexLow) -> Self {
		Self::outer_product_zero_mask(output_index)
	}

	fn basis_blades_less_than(index: &SeaBasisBladeIndexLow) -> Self{
		Self::basis_blades_less_than_mask(index)
	}

	fn basis_blades_greater_than(index: &SeaBasisBladeIndexLow) -> Self{
		Self::basis_blades_greater_than_mask(index)
	}

	fn grade(grade:&SeaGrade) -> Self {
		Self::GRADE_BITS[**grade as usize]
	}

	fn contains_basis_vec(&self, basis_vector_index: &u8) -> Self {
		*self & Self::CONTAINS_BASIS_VEC_BITS[*basis_vector_index as usize]
	}

	fn excludes_basis_vec(&self, basis_vector_index: &u8) -> Self {
		*self & Self::EXCLUDES_BASIS_VEC_BITS[*basis_vector_index as usize]
	}

	fn is_empty(&self) -> bool {
		*self == Self::ZERO
	}

	fn is_scalar(&self) -> bool {
		*self == Self::ONE
	}

	fn has(&self, index: &SeaBasisBladeIndexLow) -> bool {
		!(*self & I::from(index)).is_empty()
	}

	fn count_basis_blades(&self) -> u32 {
		self.count_ones()
	}

	fn count_basis_blades_less_than(&self, idx: &SeaBasisBladeIndexLow) -> u32{
		let mask = Self::basis_blades_less_than_mask(idx);
		(*self & mask).count_ones()
	}

	fn dual(&self) -> Self {
		self.reverse_bits()
	}

	fn dual_index(index:&SeaBasisBladeIndexLow) -> SeaBasisBladeIndexLow {
		&Self::MAX_BASIS_BLADE_INDEX - index
	}

	fn geometric_product_with_basis_blade(&self, index: &SeaBasisBladeIndexLow, metric_non_zeros: &SeaBasisBladeIndexLow) -> Self{
		let mut curr = *self;
		for i in index {
			let bv_idx = i+1;
			curr.geometric_product_with_basis_vec(bv_idx);
		}
		curr & I::outer_product_zeros(&(index | metric_non_zeros))
	}

	fn left_contraction_with_basis_blade(&self, index: &SeaBasisBladeIndexLow, metric_non_zeros: &SeaBasisBladeIndexLow) -> Self{
		let mut curr = *self >> **index;
		curr &= Self::remove_basis_blades_with_any_mask(index);
		curr & I::outer_product_zeros(&(index | metric_non_zeros))
	}

	fn right_contraction_with_basis_blade(&self, index: &SeaBasisBladeIndexLow, metric_non_zeros: &SeaBasisBladeIndexLow) -> Self{
		let mut curr = *self << **index;
		curr &= Self::remove_basis_blades_without_all_mask(index);
		curr & I::outer_product_zeros(&(index | metric_non_zeros))
	}

	fn outer_product_with_basis_blade(&self, index: &SeaBasisBladeIndexLow) -> Self {
		let curr = *self >> *Self::dual_index(index);
		curr & I::outer_product_zeros(index)
	}

	fn set(&mut self, index: &SeaBasisBladeIndexLow, value:bool) {
		if value {
			*self |= Self::from(index)
		}
		else {
			*self &= Self::from(!index)
		}
	}

	fn toggle(&mut self, index: &SeaBasisBladeIndexLow) {
		*self ^= Self::from(index);
	}

	fn iter(&self) -> Self::IntoIterType {
		SeaBasisBladeSetIter(*self)
	}

	fn intersection(&self, other:&Self) -> Self {
		*self & other
	}

	fn union(&self, other:&Self) -> Self {
		*self | other
	}

	fn symmetric_difference(&self, other:&Self) -> Self {
		*self ^ other
	}

	fn intersection_mut(&mut self, other:&Self) {
		*self &= other
	}

	fn union_mut(&mut self, other:&Self) {
		*self |= other
	}

	fn symmetric_difference_mut(&mut self, other:&Self) {
		*self ^= other
	}
}

impl<I: SeaBasisBladeSetUnsigned> LeftContraction<(&SeaBasisBladeIndexLow, u8)> for I {
	type Output = Self;

	fn left_contraction(&self, (index, num_metric_non_zeros): (&SeaBasisBladeIndexLow, u8)) -> Self::Output {
		if index.contains_basis_vectors_greater_than(num_metric_non_zeros) {
			return Self::ZERO
		}
		let mut curr = *self >> **index;
		for i in index {
			let bv_idx = i+1;
			curr &= !Self::CONTAINS_BASIS_VEC_BITS[bv_idx as usize];
		}
		curr
	}
}
impl<I: SeaBasisBladeSetUnsigned> LeftContraction<(&Self, u8)> for I {
	type Output = Self;

	fn left_contraction(&self, (rhs, num_metric_non_zeros): (&Self, u8)) -> Self {
		let mut curr = Self::ZERO;
		for i in self.iter() {
			curr |= rhs.left_contraction((&i, num_metric_non_zeros));
		}
		curr
	}
}
impl<I: SeaBasisBladeSetUnsigned> LeftContraction<(Self, u8)> for I {
	type Output = Self;

	fn left_contraction(&self, (rhs, num_metric_non_zeros): (Self, u8)) -> Self {
		self.left_contraction((&rhs, num_metric_non_zeros))
	}
}

impl<I: SeaBasisBladeSetUnsigned> RightContraction<(&SeaBasisBladeIndexLow, u8)> for I {
	type Output = Self;

	fn right_contraction(&self, (index, num_metric_non_zeros): (&SeaBasisBladeIndexLow, u8)) -> Self {
		let mut curr = self.dual() >> (I::NUM_BASIS_BLADES - 1 - **index as usize);
		for i in 0..I::DIMENSIONS {
			let bv_idx = i+1;
			if index.has(bv_idx) {
				if bv_idx > num_metric_non_zeros {
					curr &= I::CONTAINS_BASIS_VEC_BITS[bv_idx as usize];
				}
			}
			else {
				curr &= !I::CONTAINS_BASIS_VEC_BITS[bv_idx as usize];
			}
		}		
		curr
	}
}
impl<I: SeaBasisBladeSetUnsigned> RightContraction<(&Self, u8)> for I {
	type Output = Self;

	fn right_contraction(&self, (rhs, num_metric_non_zeros): (&Self, u8)) -> Self {
		rhs.left_contraction((self, num_metric_non_zeros))
	}
}
impl<I: SeaBasisBladeSetUnsigned> RightContraction<(Self, u8)> for I {
	type Output = Self;

	fn right_contraction(&self, (rhs, num_metric_non_zeros): (Self, u8)) -> Self {
		self.right_contraction((&rhs, num_metric_non_zeros))
	}
}

impl<I: SeaBasisBladeSetUnsigned> InnerProduct<(&SeaBasisBladeIndexLow, u8)> for I {
	type Output = Self;

	fn inner_product(&self, (index, num_metric_non_zeros): (&SeaBasisBladeIndexLow, u8)) -> Self {
		self.left_contraction((index, num_metric_non_zeros)) | self.right_contraction((index, num_metric_non_zeros))
	}
}
impl<I: SeaBasisBladeSetUnsigned> InnerProduct<(&Self, u8)> for I {
	type Output = Self;

	fn inner_product(&self, (rhs, num_metric_non_zeros): (&Self, u8)) -> Self {
		self.left_contraction((rhs, num_metric_non_zeros)) | self.right_contraction((rhs, num_metric_non_zeros))
	}
}
impl<I: SeaBasisBladeSetUnsigned> InnerProduct<(Self, u8)> for I {
	type Output = Self;

	fn inner_product(&self, (rhs, num_metric_non_zeros): (Self, u8)) -> Self {
		self.inner_product((&rhs, num_metric_non_zeros))
	}
}

impl<I: SeaBasisBladeSetUnsigned> OuterProduct<&SeaBasisBladeIndexLow> for I {
	type Output = Self;

	fn outer_product(&self, index: &SeaBasisBladeIndexLow) -> Self {
		let mut curr = *self << index.0;
		for i in index {
			let bv_idx = i+1;
			curr &= I::CONTAINS_BASIS_VEC_BITS[bv_idx as usize];
		}
		curr
	}
}
impl<I: SeaBasisBladeSetUnsigned> OuterProduct<&Self> for I {
	type Output = Self;

	fn outer_product(&self, rhs: &Self) -> Self {
		let mut curr = Self::ZERO;
		for i in self.iter() {
			curr |= rhs.outer_product(&i);
		}
		curr
	}
}
impl<I: SeaBasisBladeSetUnsigned> OuterProduct<Self> for I {
	type Output = Self;

	fn outer_product(&self, rhs: Self) -> Self {
		self.outer_product(&rhs)
	}
}

impl<I: SeaBasisBladeSetUnsigned> GeometricProduct<(&SeaBasisBladeIndexLow, u8)> for I {
	type Output = Self;

	fn geometric_product(&self, (index, num_metric_non_zeros): (&SeaBasisBladeIndexLow, u8)) -> Self {
		let mut curr = *self;
		for i in index {
			let bv_idx = i+1;
			curr.geometric_product_with_basis_vec(bv_idx);
			if bv_idx > num_metric_non_zeros {
				curr &= I::CONTAINS_BASIS_VEC_BITS[bv_idx as usize];
			}
		}
		curr
	}
}
impl<I: SeaBasisBladeSetUnsigned> GeometricProduct<(&Self, u8)> for I {
	type Output = Self;

	fn geometric_product(&self, (rhs, num_metric_non_zeros): (&Self, u8)) -> Self {
		let mut curr = Self::ZERO;
		for i in self.iter() {
			curr |= rhs.geometric_product((&i, num_metric_non_zeros));
		}
		curr
	}
}
impl<I: SeaBasisBladeSetUnsigned> GeometricProduct<(Self, u8)> for I {
	type Output = Self;
	
	fn geometric_product(&self, (rhs, num_metric_non_zeros): (Self, u8)) -> Self {
		self.geometric_product((&rhs, num_metric_non_zeros))
	}
}

#[cfg(test)]
mod tests {
    use std::time::Instant;

    use super::*;

	// #[test]
	// fn thing() {

	// 	// for i in 0..32{
	// 	// 	let a = SeaBasisBladeSet::<u32>::cannonical_reordering_neg_mask(&SeaBasisBladeIndexLow(i));
	// 	// 	println!("{:#x}", a.0);
	// 	// }
	
	// 	// let a = SeaBasisBladeSet::<u32>::geometric_product(&SeaBasisBladeSet(513 as u32), (&SeaBasisBladeSet(26752 as u32), 3 as u8));
	// 	// println!("{:?}", a);

	// 	let a = SeaBasisBladeSet::<u32>::outer_product(&SeaBasisBladeSet(5633 as u32), &SeaBasisBladeSet(16 as u32));
	// 	println!("{:?}", a);
	// }

	#[test]
	fn left_contraction() {
		// 1 + e1 + e2 + e13 + e123
		let a = 167_u8;
		// e12
		let b = 8_u8;
		// e1 + e2 + e12
		let out = 14_u8;
		assert_eq!(
			a.left_contraction((&b, 3)),
			out
		);

		// 1 + e13 + e123
		let a = 161_u8;
		// e1 + e12 + e13 + e23 
		let b = 106_u8;
		// 1 + e1 + e12 + e13 + e23
		let out = 107_u8;
		assert_eq!(
			a.left_contraction((&b, 3)),
			out
		);

		// e1 + e2
		let a = 6_u8;
		// e2 + e3
		let b = 20_u8;
		// 1
		let out = 1_u8;
		assert_eq!(
			a.left_contraction((&b, 3)),
			out
		);

		// e1 + e12 + e3 + e134 + e234
		let a = 24602_u16;
		// e34
		let b = 1_u16 << 12;
		// e4
		let out = 256_u16;
		assert_eq!(
			a.left_contraction((&b, 3)),
			out
		);

		// 1 + e13 + e4 + e123
		let a = 417_u16;
		// e1 + e12 + e134 + e234
		let b = 24586_u16;
		// e1 + e12 + e4 + e134 + e234 
		let out = 24842_u16;
		assert_eq!(
			a.left_contraction((&b, 3)),
			out
		);

		// e3 + e14 + e124
		let a = 2576_u16;
		// e24 + e134
		let b = 9216_u16;
		// e14
		let out = 512_u16;
		assert_eq!(
			a.left_contraction((&b, 3)),
			out
		);
	}

	#[test]
	fn right_contraction() {
		// 1 + e1 + e2 + e13 + e123
		let a = 167_u8;
		// e12
		let b = 8_u8;
		// e3
		let out = 16_u8;
		assert_eq!(
			a.right_contraction((&b, 3)),
			out
		);

		// 1 + e13 + e123
		let a = 161_u8;
		// e1 + e12 + e13 + e23 
		let b = 106_u8;
		// 1 + e1 + e2 + e3 + e23
		let out = 87_u8;
		assert_eq!(
			a.right_contraction((&b, 3)),
			out
		);

		// e1 + e2
		let a = 6_u8;
		// e2 + e3
		let b = 20_u8;
		// 1
		let out = 1_u8;
		assert_eq!(
			a.right_contraction((&b, 3)),
			out
		);

		// e1 + e12 + e3 + e134 + e234
		let a = 24602_u16;
		// e34
		let b = 1_u16 << 12;
		// 0
		let out = 0_u16;
		assert_eq!(
			a.right_contraction((&b, 3)),
			out
		);

		// 1 + e13 + e4 + e123
		let a = 417_u16;
		// e1 + e12 + e134 + e234
		let b = 24586_u16;
		// e3 + e23
		let out = 80_u16;
		assert_eq!(
			a.right_contraction((&b, 3)),
			out
		);

		// e3 + e14 + e124
		let a = 2576_u16;
		// e24 + e134
		let b = 9216_u16;
		// 0
		let out = 0_u16;
		assert_eq!(
			a.right_contraction((&b, 3)),
			out
		);
	}

	#[test]
	fn inner_product() {
		// 1 + e1 + e2 + e13 + e123
		let a = 167_u8;
		// e12
		let b = 8_u8;
		// e1 + e2 + e12 + e3
		let out = 30_u8;
		
		assert_eq!(
			a.inner_product((&b, 3)),
			out
		);

		// 1 + e13 + e123
		let a = 161_u8;
		// e1 + e12 + e13 + e23 
		let b = 106_u8;
		// 1 + e2 + e12 + e3 + e13 + e23
		let out = 127_u8;
		assert_eq!(
			a.inner_product((&b, 3)),
			out
		);

		// e1 + e2
		let a = 6_u8;
		// e2 + e3
		let b = 20_u8;
		// 1
		let out = 1_u8;
		assert_eq!(
			a.inner_product((&b, 3)),
			out
		);

		// e1 + e12 + e3 + e134 + e234
		let a = 24602_u16;
		// e34
		let b = 1_u16 << 12;
		// e4
		let out = 256_u16;
		assert_eq!(
			a.inner_product((&b, 3)),
			out
		);

		// 1 + e13 + e4 + e123
		let a = 417_u16;
		// e1 + e12 + e134 + e234
		let b = 24586_u16;
		// e1 + e12 + e3  + e23 + e4 +  e134 + e234
		let out = 24922_u16;
		assert_eq!(
			a.inner_product((&b, 3)),
			out
		);

		// e3 + e14 + e124
		let a = 2576_u16;
		// e24 + e134
		let b = 9216_u16;
		// e14
		let out = 512_u16;
		assert_eq!(
			a.inner_product((&b, 3)),
			out
		);

		// 3d
		// usual 1.228µs
		// no right 1.162µs (105% faster)
		// 4d
		// usual 2.589µs
		// no right 1.948µs (132% faster)
		// 5d
		// usual 5.299µs
		// no right 3.488µs (151% faster)
		// 6d
		// usual 11.643µs
		// no right 6.624µs (175% faster)
		// let a = (0_u64).overflowing_sub(1).0;
		// let b = (0_u64).overflowing_sub(1).0;

		// let now = Instant::now();
		// let num_iters = 10_000;
		// let mut t = a.clone();
		// for _ in 0..num_iters {
		// 	t = t.inner_product((&b, 3));
		// }
		// let elapsed_time = now.elapsed();
		// println!("time time time {:?}", elapsed_time/num_iters);
	}

	#[test]
	fn outer_product() {
		// 1 + e1 + e2 + e13 + e123
		let a = 167_u8;
		// e12
		let b = 8_u8;
		// e12
		let out = 8_u8;
		assert_eq!(
			a.outer_product(&b),
			out
		);

		// 1 + e13 + e123
		let a = 161_u8;
		// e1 + e12 + e13 + e23 
		let b = 106_u8;
		// e1 + e12 + e13 + e23
		let out = 106_u8;
		assert_eq!(
			a.outer_product(&b),
			out
		);

		// e1 + e2
		let a = 6_u8;
		// e2 + e3
		let b = 20_u8;
		// e12 + e13 + e23
		let out = 104_u8;
		assert_eq!(
			a.outer_product(&b),
			out
		);

		// e1 + e12 + e3 + e134 + e234
		let a = 24602_u16;
		// e34
		let b = 1_u16 << 12;
		// e134 + e1234
		let out = 40960_u16;
		assert_eq!(
			a.outer_product(&b),
			out
		);

		// 1 + e13 + e4 + e123
		let a = 417_u16;
		// e1 + e12 + e134 + e234
		let b = 24586_u16;
		// e1 + e14 + e12 + e124 + e134 + e234
		let out = 27146_u16;
		assert_eq!(
			a.outer_product(&b),
			out
		);

		// e3 + e14 + e124
		let a = 2576_u16;
		// e24 + e134
		let b = 9216_u16;
		// e234
		let out = 16384_u16;
		assert_eq!(
			a.outer_product(&b),
			out
		);
	}

	#[test]
	fn geometric_product() {
		// 1 + e1 + e2 + e13 + e123
		let a = 167_u8;
		// e12
		let b = 8_u8;
		// e1 + e2 - e12 + e3 + e23 
		let out = 94_u8;
		assert_eq!(
			a.geometric_product((&b, 3)),
			out
		);

		// 1 + e13 + e123
		let a = 161_u8;
		// e1 + e12 + e13 + e23 
		let b = 106_u8;
		// 1 + e1 + e2 + e12 + e3 + e13 + e23
		let out = 127_u8;
		assert_eq!(
			a.geometric_product((&b, 3)),
			out
		);

		// e1 + e2
		let a = 6_u8;
		// e2 + e3
		let b = 20_u8;
		// 1 + e12 + e13 + e23
		let out = 105_u8;
		assert_eq!(
			a.geometric_product((&b, 3)),
			out
		);

		// e1 + e12 + e3 + e134 + e234
		let a = 24602_u16;
		// e34
		let b = 1_u16 << 12;
		// e4 + e134 + e1234
		let out = 41216_u16;
		assert_eq!(
			a.geometric_product((&b, 3)),
			out
		);

		// 1 + e13 + e4 + e123
		let a = 417_u16;
		// e1 + e12 + e134 + e234
		let b = 24586_u16;
		// e4 + e1 + e3 + e14 + e24 + e12 + e23 + e124 + e134 + e234
		let out = 28506_u16;
		assert_eq!(
			a.geometric_product((&b, 3)),
			out
		);

		// e3 + e14 + e124
		let a = 2576_u16;
		// e24 + e134
		let b = 9216_u16;
		// e14 + e234
		let out = 16896_u16;
		assert_eq!(
			a.geometric_product((&b, 3)),
			out
		);
	}
}
