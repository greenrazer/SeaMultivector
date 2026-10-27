//! Basis blades: their indices, sets of them, and the sign tables behind every product.

use std::fmt;
use std::ops::{BitAnd, BitOr, BitXor, Not};
use std::str::FromStr;

use crate::error::ParseError;

mod set;
mod tables;

pub use set::{BladeSet, SetBlades};

/// The index of a basis blade.
///
/// Bit k is set when the blade contains basis vector e(k+1): e1 is 1, e2 is 2, e12 is 3, e3 is 4.
/// Multiplying two basis blades lands on the blade whose index is the XOR of theirs. An index
/// means the same blade in every dimension, so 3D blades are the indices below 8.
///
/// # Examples
///
/// ```
/// use seagrass_math::BladeIndex;
///
/// assert_eq!(BladeIndex::E12.bits(), 0b011);
/// assert_eq!(BladeIndex::E12.grade(), 2);
/// assert_eq!(BladeIndex::E1.bits() ^ BladeIndex::E12.bits(), BladeIndex::E2.bits());
/// assert_eq!(BladeIndex::E134.to_string(), "e134");
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Default)]
pub struct BladeIndex(u8);

macro_rules! blade {
    ($name:ident, $bits:literal) => {
        #[doc = concat!("The basis blade `", stringify!($name), "`.")]
        pub const $name: Self = Self($bits);
    };
}

impl BladeIndex {
    /// The scalar blade, 1.
    pub const SCALAR: Self = Self(0);
    blade!(E1, 1);
    blade!(E2, 2);
    blade!(E12, 3);
    blade!(E3, 4);
    blade!(E13, 5);
    blade!(E23, 6);
    blade!(E123, 7);
    blade!(E4, 8);
    blade!(E14, 9);
    blade!(E24, 10);
    blade!(E124, 11);
    blade!(E34, 12);
    blade!(E134, 13);
    blade!(E234, 14);
    blade!(E1234, 15);
    blade!(E5, 16);
    blade!(E15, 17);
    blade!(E25, 18);
    blade!(E125, 19);
    blade!(E35, 20);
    blade!(E135, 21);
    blade!(E235, 22);
    blade!(E1235, 23);
    blade!(E45, 24);
    blade!(E145, 25);
    blade!(E245, 26);
    blade!(E1245, 27);
    blade!(E345, 28);
    blade!(E1345, 29);
    blade!(E2345, 30);
    blade!(E12345, 31);
    blade!(E6, 32);
    blade!(E16, 33);
    blade!(E26, 34);
    blade!(E126, 35);
    blade!(E36, 36);
    blade!(E136, 37);
    blade!(E236, 38);
    blade!(E1236, 39);
    blade!(E46, 40);
    blade!(E146, 41);
    blade!(E246, 42);
    blade!(E1246, 43);
    blade!(E346, 44);
    blade!(E1346, 45);
    blade!(E2346, 46);
    blade!(E12346, 47);
    blade!(E56, 48);
    blade!(E156, 49);
    blade!(E256, 50);
    blade!(E1256, 51);
    blade!(E356, 52);
    blade!(E1356, 53);
    blade!(E2356, 54);
    blade!(E12356, 55);
    blade!(E456, 56);
    blade!(E1456, 57);
    blade!(E2456, 58);
    blade!(E12456, 59);
    blade!(E3456, 60);
    blade!(E13456, 61);
    blade!(E23456, 62);
    blade!(E123456, 63);

    /// The blade with the given bitmask.
    ///
    /// # Panics
    ///
    /// Panics if `bits` is 64 or more, which no algebra of up to 6 dimensions contains.
    pub const fn new(bits: u8) -> Self {
        assert!(bits < 64, "blade bitmask must be below 64");
        Self(bits)
    }

    /// The bitmask: bit k is set when the blade contains e(k+1).
    pub const fn bits(self) -> u8 {
        self.0
    }

    /// The number of basis vectors in the blade.
    pub const fn grade(self) -> u32 {
        self.0.count_ones()
    }

    /// Whether the blade contains basis vector `e{vector}`, numbered from 1.
    pub const fn contains(self, vector: u8) -> bool {
        vector >= 1 && vector <= 8 && self.0 >> (vector - 1) & 1 == 1
    }

    /// The blade's basis vectors in ascending order, numbered from 1.
    pub fn vectors(self) -> impl Iterator<Item = u8> {
        self.positions().map(|k| k + 1)
    }

    /// The blade's bit positions in ascending order: vector e(k+1) yields k.
    pub(crate) fn positions(self) -> Positions {
        Positions(self.0)
    }

    /// The blade e1..en holding the first `n` basis vectors.
    pub(crate) const fn first_vectors(n: u8) -> Self {
        if n >= 8 {
            Self(u8::MAX)
        } else {
            Self((1 << n) - 1)
        }
    }

    /// Whether the blade holds a basis vector numbered above `n`.
    pub(crate) const fn has_vector_above(self, n: u8) -> bool {
        self.0 & !Self::first_vectors(n).0 != 0
    }

    /// Parses the digits after `e`, such as `"124"` for e124.
    pub(crate) fn from_digits(digits: &str) -> Result<Self, ParseError> {
        if digits.is_empty() {
            return Err(ParseError::InvalidBasisBladeIndex);
        }
        let mut bits = 0_u8;
        let mut last = 0;
        for c in digits.chars() {
            let vector = match c.to_digit(10) {
                Some(0) => return Err(ParseError::InvalidBasisVectorZero),
                Some(7..) => return Err(ParseError::InvalidBasisVectorOutOfRange),
                Some(d) => d as u8,
                None => return Err(ParseError::InvalidBasisVector),
            };
            if vector == last {
                return Err(ParseError::DuplicateBasisVector);
            }
            if vector < last {
                return Err(ParseError::InvalidBasisVectorOrdering);
            }
            last = vector;
            bits |= 1 << (vector - 1);
        }
        Ok(Self(bits))
    }
}

impl fmt::Display for BladeIndex {
    /// Writes `1` for the scalar and `e` followed by the vector numbers otherwise.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.0 == 0 {
            return f.write_str("1");
        }
        f.write_str("e")?;
        for vector in self.vectors() {
            write!(f, "{vector}")?;
        }
        Ok(())
    }
}

impl FromStr for BladeIndex {
    type Err = ParseError;

    /// Reads the [`Display`](fmt::Display) form: `1` or `e` followed by ascending vector numbers.
    fn from_str(s: &str) -> Result<Self, ParseError> {
        match s {
            "1" => Ok(Self::SCALAR),
            _ => match s.strip_prefix('e') {
                Some(digits) => Self::from_digits(digits),
                None => Err(ParseError::UnexpectedToken { position: 0 }),
            },
        }
    }
}

impl BitXor for BladeIndex {
    type Output = Self;

    fn bitxor(self, rhs: Self) -> Self {
        Self(self.0 ^ rhs.0)
    }
}

impl BitOr for BladeIndex {
    type Output = Self;

    fn bitor(self, rhs: Self) -> Self {
        Self(self.0 | rhs.0)
    }
}

impl BitAnd for BladeIndex {
    type Output = Self;

    fn bitand(self, rhs: Self) -> Self {
        Self(self.0 & rhs.0)
    }
}

impl Not for BladeIndex {
    type Output = Self;

    /// Flips all 8 bits; only useful masked by another index.
    fn not(self) -> Self {
        Self(!self.0)
    }
}

/// Iterator over the set bit positions of a [`BladeIndex`], lowest first.
#[derive(Debug, Clone)]
pub(crate) struct Positions(u8);

impl Iterator for Positions {
    type Item = u8;

    fn next(&mut self) -> Option<u8> {
        if self.0 == 0 {
            return None;
        }
        let position = self.0.trailing_zeros() as u8;
        self.0 &= self.0 - 1;
        Some(position)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_digits_and_rejects_malformed_blades() {
        assert_eq!(BladeIndex::from_digits("124"), Ok(BladeIndex::E124));
        let cases = [
            ("", ParseError::InvalidBasisBladeIndex),
            ("102", ParseError::InvalidBasisVectorZero),
            ("17", ParseError::InvalidBasisVectorOutOfRange),
            ("1x", ParseError::InvalidBasisVector),
            ("11", ParseError::DuplicateBasisVector),
            ("21", ParseError::InvalidBasisVectorOrdering),
        ];
        for (digits, err) in cases {
            assert_eq!(BladeIndex::from_digits(digits), Err(err), "{digits}");
        }
    }

    #[test]
    fn round_trips_through_display() {
        for bits in 0..64 {
            let blade = BladeIndex::new(bits);
            assert_eq!(blade.to_string().parse::<BladeIndex>(), Ok(blade));
        }
    }

    #[test]
    fn first_vectors_and_vectors_above() {
        assert_eq!(BladeIndex::first_vectors(3), BladeIndex::E123);
        assert_eq!(BladeIndex::first_vectors(8).bits(), 0xFF);
        assert!(BladeIndex::E14.has_vector_above(3));
        assert!(!BladeIndex::E13.has_vector_above(3));
    }
}
