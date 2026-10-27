//! The membership and sign masks behind every product, written once as `u64` for 6D.
//!
//! Bit i of a mask describes the blade with index i. An index means the same blade in every
//! dimension, so the mask for an n-dimensional algebra is the low 2^n bits of the 6D mask, which
//! `as u8`, `as u16` or `as u32` produces. The exceptions are the complement and dual signs, which
//! depend on the pseudoscalar and so on the dimension; they are computed per dimension from the
//! reordering sign below.

/// Blades of even grade (an even number of set bits). Its complement, the odd grades, is the
/// Thue–Morse sequence.
pub(crate) const EVEN: u64 = 0x9669_6996_6996_9669;

/// Blades that reversion negates: grade g with g(g - 1)/2 odd, i.e. grades 2 and 3 modulo 4.
pub(crate) const REVERSE_SIGNS: u64 = 0x8117_177E_177E_7EE8;

/// Blades that Clifford conjugation negates: grade g with g(g + 1)/2 odd, i.e. grades 1 and 2
/// modulo 4. It is reversion combined with grade involution.
pub(crate) const CONJUGATE_SIGNS: u64 = 0xE881_8117_8117_177E;

/// `GRADE[g]` holds the blades of grade g.
pub(crate) const GRADE: [u64; 7] = [
    0x0000_0000_0000_0001,
    0x0000_0001_0001_0116,
    0x0001_0116_0116_1668,
    0x0116_1668_1668_6880,
    0x1668_6880_6880_8000,
    0x6880_8000_8000_0000,
    0x8000_0000_0000_0000,
];

/// `CONTAINS[k]` holds the blades that contain basis vector e(k) for k in 1..=6: bit k - 1 of
/// the index is set, so the pattern is runs of 2^(k-1) zeros and ones. `CONTAINS[0]` is every
/// blade.
pub(crate) const CONTAINS: [u64; 7] = [
    0xFFFF_FFFF_FFFF_FFFF,
    0xAAAA_AAAA_AAAA_AAAA,
    0xCCCC_CCCC_CCCC_CCCC,
    0xF0F0_F0F0_F0F0_F0F0,
    0xFF00_FF00_FF00_FF00,
    0xFFFF_0000_FFFF_0000,
    0xFFFF_FFFF_0000_0000,
];

/// `REORDER[k]` holds the blades that contain a basis vector numbered above k.
///
/// The reordering sign of the pair (x, x XOR o) is the reversion sign of x, flipped once for every
/// pair of a vector in o and a higher vector in x. XORing `REORDER[k]` into the reversion mask for
/// each vector e(k) of o therefore gives every pair's sign for output blade o at once.
pub(crate) const REORDER: [u64; 7] = {
    let mut out = [0; 7];
    let mut k = 0;
    while k < 7 {
        let mut higher = k + 1;
        while higher < 7 {
            out[k] ^= CONTAINS[higher];
            higher += 1;
        }
        k += 1;
    }
    out
};

/// Whether moving the vectors of blade `a` past those of blade `b` into ascending order takes an
/// odd number of swaps, i.e. whether `a * b` has a negative sign before the metric is applied.
pub(crate) const fn reorder_negative(a: u64, b: u64) -> bool {
    let mut a = a >> 1;
    let mut swaps = 0;
    while a != 0 {
        swaps += (a & b).count_ones();
        a >>= 1;
    }
    swaps % 2 == 1
}

/// The three sign tables that depend on the dimension.
#[derive(Clone, Copy)]
enum Signs {
    RightComplement,
    LeftComplement,
    Dual,
}

/// Collects the output blades c of a `dim`-dimensional algebra whose sign is negative.
const fn signs_by_output(dim: u8, signs: Signs) -> u64 {
    let max = (1_u64 << dim) - 1;
    let mut out = 0;
    let mut c = 0;
    while c <= max {
        let negative = match signs {
            Signs::RightComplement => reorder_negative(c, max - c),
            Signs::LeftComplement => reorder_negative(max - c, c),
            Signs::Dual => reorder_negative(max - c, max),
        };
        if negative {
            out |= 1 << c;
        }
        c += 1;
    }
    out
}

/// Signs of the right complement: after bit reversal moves x to its complement c, c is negated
/// where `c ^ x = -I`, so the result r satisfies `r ^ x = I`.
pub(crate) const fn right_complement_signs(dim: u8) -> u64 {
    signs_by_output(dim, Signs::RightComplement)
}

/// Signs of the left complement: the result l satisfies `x ^ l = I`.
pub(crate) const fn left_complement_signs(dim: u8) -> u64 {
    signs_by_output(dim, Signs::LeftComplement)
}

/// Signs of `x * I` in a metric with only positive vectors: x lands on its complement c, with the
/// reordering sign of moving x past the pseudoscalar.
pub(crate) const fn dual_signs(dim: u8) -> u64 {
    signs_by_output(dim, Signs::Dual)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn grade_mask(keep: impl Fn(u32) -> bool) -> u64 {
        (0..64)
            .filter(|&i: &u64| keep(i.count_ones()))
            .fold(0, |m, i| m | 1 << i)
    }

    #[test]
    fn grade_masks_match_their_definitions() {
        assert_eq!(EVEN, grade_mask(|g| g % 2 == 0));
        assert_eq!(
            REVERSE_SIGNS,
            grade_mask(|g| g * g.saturating_sub(1) / 2 % 2 == 1)
        );
        assert_eq!(CONJUGATE_SIGNS, grade_mask(|g| g * (g + 1) / 2 % 2 == 1));
        for (g, mask) in GRADE.iter().enumerate() {
            assert_eq!(*mask, grade_mask(|h| h == g as u32), "grade {g}");
        }
    }

    #[test]
    fn contains_masks_match_their_definitions() {
        assert_eq!(CONTAINS[0], u64::MAX);
        for (k, mask) in CONTAINS.iter().enumerate().skip(1) {
            let expected = (0..64_u64)
                .filter(|i| i >> (k - 1) & 1 == 1)
                .fold(0, |m, i| m | 1 << i);
            assert_eq!(*mask, expected, "e{k}");
        }
    }

    /// The literals the crate used before the tables were derived, with the 6D left complement
    /// corrected (it had four wrong bits: blades 5, 6, 53 and 54).
    #[test]
    fn derived_complement_and_dual_signs_match_the_literals() {
        let literals = [
            (3, 0x24, 0x24, 0x33),
            (4, 0x4D24, 0x24B2, 0x33CC),
            (5, 0x24B2_4D24, 0x24B2_4D24, 0x33CC_33CC),
            (
                6,
                0x4D24_DB4D_24B2_4D24,
                0x24B2_4D24_B2DB_24B2,
                0x33CC_33CC_CC33_CC33,
            ),
        ];
        for (dim, right, left, dual) in literals {
            assert_eq!(right_complement_signs(dim), right, "right {dim}D");
            assert_eq!(left_complement_signs(dim), left, "left {dim}D");
            assert_eq!(dual_signs(dim), dual, "dual {dim}D");
        }
    }
}
