//! Metric signatures: how many basis vectors square to +1, -1 and 0.

use std::fmt;
use std::marker::PhantomData;
use std::str::FromStr;

use crate::blade::{BladeIndex, BladeSet};
use crate::error::{MetricError, ParseError};

mod sealed {
    pub trait Sealed {}
}

/// A metric signature (p, q, r): p basis vectors square to +1, q to -1 and r to 0.
///
/// Basis vectors are ordered positive, then negative, then null, so in (3, 0, 1) the null vector
/// is e4 (the vector many PGA texts call e0). The dimension p + q + r picks the bitset type
/// [`Metric::Bits`]. The trait is implemented by [`Signature`], fixed at compile time, and
/// [`Dynamic`], chosen at runtime, and cannot be implemented outside the crate.
///
/// # Examples
///
/// ```
/// use seagrass_math::prelude::*;
///
/// let sta = Dyn4::new(1, 3, 0).unwrap();
/// assert_eq!((sta.positive(), sta.negative(), sta.zero()), (1, 3, 0));
/// assert_eq!(Pga3::default().to_string(), "3,0,1");
/// ```
pub trait Metric:
    Copy + PartialEq + fmt::Debug + fmt::Display + Send + Sync + 'static + sealed::Sealed
{
    /// The bitset type that holds a set of blades of this dimension.
    type Bits: BladeSet;

    /// The number of basis vectors that square to +1.
    fn positive(&self) -> u8;

    /// The number of basis vectors that square to -1.
    fn negative(&self) -> u8;

    /// The number of basis vectors that square to 0.
    fn zero(&self) -> u8;

    /// The number of basis vectors, p + q + r.
    fn dimensions(&self) -> u8 {
        Self::Bits::DIM
    }

    /// The metric of a parsed value, from the `p,q,r` text before `>` if there is one.
    ///
    /// # Errors
    ///
    /// [`ParseError::MissingMetric`] when a runtime metric has no prefix,
    /// [`ParseError::InvalidMetric`] or [`ParseError::InvalidMetricDimensions`] when the prefix is
    /// not a signature of this dimension, and [`ParseError::MismatchedMetric`] when it differs
    /// from a compile-time signature.
    fn from_prefix(prefix: Option<&str>) -> Result<Self, ParseError>;
}

/// The counts and masks the product loops need from a metric.
#[derive(Debug, Clone, Copy)]
pub(crate) struct ProductConstants {
    /// p + q: vectors numbered above this are null.
    pub(crate) non_null_count: u8,
    /// The blade holding every non-null vector.
    pub(crate) non_null: BladeIndex,
    /// The blade holding every negative vector.
    pub(crate) negative: BladeIndex,
}

impl ProductConstants {
    #[inline]
    pub(crate) fn of<M: Metric>(metric: &M) -> Self {
        let non_null_count = metric.positive() + metric.negative();
        let non_null = BladeIndex::first_vectors(non_null_count);
        let positive = BladeIndex::first_vectors(metric.positive());
        Self {
            non_null_count,
            non_null,
            negative: non_null ^ positive,
        }
    }
}

/// Whether the pseudoscalar's inverse is its negation: I * I = -1 when the reversion swaps of I
/// plus its negative vectors are odd.
pub(crate) fn pseudoscalar_inverse_negative<M: Metric>(metric: &M) -> bool {
    let n = u32::from(metric.dimensions());
    (n * (n - 1) / 2 + u32::from(metric.negative())) % 2 == 1
}

/// Parses `p`, `p,q` or `p,q,r`; missing counts are zero.
pub(crate) fn parse_signature(s: &str) -> Result<(u8, u8, u8), MetricError> {
    if s.trim().is_empty() {
        return Err(MetricError::MissingNumPositive);
    }
    let mut parts = s.split(',').map(str::trim);
    let mut count = |err| {
        parts
            .next()
            .map(|t| t.parse::<u8>().map_err(|_| err))
            .transpose()
    };
    let p = count(MetricError::InvalidNumPositive)?.unwrap_or(0);
    let q = count(MetricError::InvalidNumNegative)?.unwrap_or(0);
    let r = count(MetricError::InvalidNumZero)?.unwrap_or(0);
    if parts.next().is_some() {
        return Err(MetricError::InvalidMetricLength);
    }
    Ok((p, q, r))
}

/// A signature fixed at compile time; a zero-sized type.
///
/// Implemented for every p + q + r from 3 to 6. Operations between values of different
/// signatures do not compile, and the metric costs nothing at runtime.
///
/// # Examples
///
/// ```
/// use seagrass_math::prelude::*;
///
/// let a = Multivector::<f64, Signature<4, 1, 0>>::from_str("e4 + e5").unwrap();
/// assert_eq!((&a * &a).scalar(), 0.0); // e4 + e5 is a null vector in conformal 3D
/// assert_eq!(std::mem::size_of::<Signature<4, 1, 0>>(), 0);
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub struct Signature<const P: u8, const Q: u8, const R: u8>;

impl<const P: u8, const Q: u8, const R: u8> fmt::Display for Signature<P, Q, R> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{P},{Q},{R}")
    }
}

impl<const P: u8, const Q: u8, const R: u8> sealed::Sealed for Signature<P, Q, R> {}

macro_rules! signatures {
    ($bits:ty: $(($p:literal, $q:literal, $r:literal))*) => {$(
        impl Metric for Signature<$p, $q, $r> {
            type Bits = $bits;

            #[inline]
            fn positive(&self) -> u8 {
                $p
            }

            #[inline]
            fn negative(&self) -> u8 {
                $q
            }

            #[inline]
            fn zero(&self) -> u8 {
                $r
            }

            fn from_prefix(prefix: Option<&str>) -> Result<Self, ParseError> {
                match prefix {
                    None => Ok(Self),
                    Some(s) if parse_signature(s)? == ($p, $q, $r) => Ok(Self),
                    Some(_) => Err(ParseError::MismatchedMetric),
                }
            }
        }
    )*};
}

signatures!(u8: (3, 0, 0) (2, 1, 0) (2, 0, 1) (1, 2, 0) (1, 1, 1) (1, 0, 2) (0, 3, 0) (0, 2, 1) (0, 1, 2) (0, 0, 3));
signatures!(u16: (4, 0, 0) (3, 1, 0) (3, 0, 1) (2, 2, 0) (2, 1, 1) (2, 0, 2) (1, 3, 0) (1, 2, 1) (1, 1, 2) (1, 0, 3) (0, 4, 0) (0, 3, 1) (0, 2, 2) (0, 1, 3) (0, 0, 4));
signatures!(u32: (5, 0, 0) (4, 1, 0) (4, 0, 1) (3, 2, 0) (3, 1, 1) (3, 0, 2) (2, 3, 0) (2, 2, 1) (2, 1, 2) (2, 0, 3) (1, 4, 0) (1, 3, 1) (1, 2, 2) (1, 1, 3) (1, 0, 4) (0, 5, 0) (0, 4, 1) (0, 3, 2) (0, 2, 3) (0, 1, 4) (0, 0, 5));
signatures!(u64: (6, 0, 0) (5, 1, 0) (5, 0, 1) (4, 2, 0) (4, 1, 1) (4, 0, 2) (3, 3, 0) (3, 2, 1) (3, 1, 2) (3, 0, 3) (2, 4, 0) (2, 3, 1) (2, 2, 2) (2, 1, 3) (2, 0, 4) (1, 5, 0) (1, 4, 1) (1, 3, 2) (1, 2, 3) (1, 1, 4) (1, 0, 5) (0, 6, 0) (0, 5, 1) (0, 4, 2) (0, 3, 3) (0, 2, 4) (0, 1, 5) (0, 0, 6));

/// A signature chosen at runtime, for algebras of the dimension of the bitset type `B`.
///
/// Three bytes. Binary operations check that both operands share the signature and panic if not.
///
/// # Examples
///
/// ```
/// use seagrass_math::prelude::*;
///
/// let metric: Dyn4 = "3,0,1".parse().unwrap();
/// assert_eq!(metric, Dyn4::new(3, 0, 1).unwrap());
/// let e1 = Multivector::<f64, Dyn4>::basis(metric, BladeIndex::E1);
/// assert_eq!(e1.to_string(), "3,0,1>e1");
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Dynamic<B> {
    positive: u8,
    negative: u8,
    zero: u8,
    bits: PhantomData<B>,
}

impl<B: BladeSet> Dynamic<B> {
    /// The signature (p, q, r).
    ///
    /// # Errors
    ///
    /// [`MetricError::WrongDimension`] if p + q + r differs from the dimension of `B`.
    pub fn new(p: u8, q: u8, r: u8) -> Result<Self, MetricError> {
        let found = u16::from(p) + u16::from(q) + u16::from(r);
        if found != u16::from(B::DIM) {
            return Err(MetricError::WrongDimension {
                expected: B::DIM,
                found: found.min(255) as u8,
            });
        }
        Ok(Self {
            positive: p,
            negative: q,
            zero: r,
            bits: PhantomData,
        })
    }
}

impl<B> fmt::Display for Dynamic<B> {
    /// Writes `p,q,r`, the syntax [`FromStr`] reads.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{},{},{}", self.positive, self.negative, self.zero)
    }
}

impl<B: BladeSet> FromStr for Dynamic<B> {
    type Err = MetricError;

    /// Reads `p`, `p,q` or `p,q,r`; missing counts are zero.
    fn from_str(s: &str) -> Result<Self, MetricError> {
        let (p, q, r) = parse_signature(s)?;
        Self::new(p, q, r)
    }
}

impl<B> sealed::Sealed for Dynamic<B> {}

impl<B: BladeSet> Metric for Dynamic<B> {
    type Bits = B;

    #[inline]
    fn positive(&self) -> u8 {
        self.positive
    }

    #[inline]
    fn negative(&self) -> u8 {
        self.negative
    }

    #[inline]
    fn zero(&self) -> u8 {
        self.zero
    }

    fn from_prefix(prefix: Option<&str>) -> Result<Self, ParseError> {
        let (p, q, r) = parse_signature(prefix.ok_or(ParseError::MissingMetric)?)?;
        Self::new(p, q, r).map_err(|_| ParseError::InvalidMetricDimensions)
    }
}

/// Euclidean 3D, (3, 0, 0).
pub type Euclid3 = Signature<3, 0, 0>;
/// Euclidean 4D, (4, 0, 0).
pub type Euclid4 = Signature<4, 0, 0>;
/// Euclidean 5D, (5, 0, 0).
pub type Euclid5 = Signature<5, 0, 0>;
/// Euclidean 6D, (6, 0, 0).
pub type Euclid6 = Signature<6, 0, 0>;
/// Projective geometric algebra of the plane, (2, 0, 1).
pub type Pga2 = Signature<2, 0, 1>;
/// Projective geometric algebra of space, (3, 0, 1); e4 is the null vector.
pub type Pga3 = Signature<3, 0, 1>;
/// Projective geometric algebra of 4D space, (4, 0, 1).
pub type Pga4 = Signature<4, 0, 1>;
/// Conformal geometric algebra of the plane, (3, 1, 0).
pub type Cga2 = Signature<3, 1, 0>;
/// Conformal geometric algebra of space, (4, 1, 0).
pub type Cga3 = Signature<4, 1, 0>;
/// Spacetime algebra, (1, 3, 0).
pub type Sta = Signature<1, 3, 0>;
/// A runtime signature of dimension 3.
pub type Dyn3 = Dynamic<u8>;
/// A runtime signature of dimension 4.
pub type Dyn4 = Dynamic<u16>;
/// A runtime signature of dimension 5.
pub type Dyn5 = Dynamic<u32>;
/// A runtime signature of dimension 6.
pub type Dyn6 = Dynamic<u64>;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_signatures_and_reports_each_error() {
        assert_eq!(parse_signature("3"), Ok((3, 0, 0)));
        assert_eq!(parse_signature(" 1, 3 "), Ok((1, 3, 0)));
        assert_eq!(parse_signature("3,0,1"), Ok((3, 0, 1)));
        let cases = [
            ("", MetricError::MissingNumPositive),
            ("a", MetricError::InvalidNumPositive),
            ("3,b", MetricError::InvalidNumNegative),
            ("3,0,c", MetricError::InvalidNumZero),
            ("3,0,1,0", MetricError::InvalidMetricLength),
        ];
        for (s, err) in cases {
            assert_eq!(parse_signature(s), Err(err), "{s:?}");
        }
    }

    #[test]
    fn product_constants_mark_null_and_negative_vectors() {
        let c = ProductConstants::of(&Signature::<2, 1, 1>);
        assert_eq!(c.non_null_count, 3);
        assert_eq!(c.non_null, BladeIndex::E123);
        assert_eq!(c.negative, BladeIndex::E3);
    }

    #[test]
    fn compile_time_signature_accepts_only_its_own_prefix() {
        assert_eq!(Pga3::from_prefix(None), Ok(Pga3::default()));
        assert_eq!(Pga3::from_prefix(Some("3,0,1")), Ok(Pga3::default()));
        assert_eq!(
            Pga3::from_prefix(Some("4")),
            Err(ParseError::MismatchedMetric)
        );
    }
}
