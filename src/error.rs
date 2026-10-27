//! Error types.

use std::error::Error;
use std::fmt;

/// Why a multivector could not be parsed from a string.
///
/// # Examples
///
/// ```
/// use seagrass_math::prelude::*;
///
/// let err = Multivector::<f64, Dyn3>::from_str("3>e31").unwrap_err();
/// assert_eq!(err, ParseError::InvalidBasisVectorOrdering);
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ParseError {
    /// A coefficient is not a number.
    InvalidBasisBladeCoefficient,
    /// A blade is written as `e` with no basis vectors after it.
    InvalidBasisBladeIndex,
    /// A blade contains a character that is not a digit.
    InvalidBasisVector,
    /// A blade contains the basis vector `0`; vectors are numbered from 1.
    InvalidBasisVectorZero,
    /// A blade lists its basis vectors in descending order, as in `e21`.
    InvalidBasisVectorOrdering,
    /// A blade contains a basis vector above 6.
    InvalidBasisVectorOutOfRange,
    /// A blade repeats a basis vector, as in `e11`.
    DuplicateBasisVector,
    /// A blade uses a basis vector beyond the algebra's dimension.
    InvalidBasisBladeForAlgebra,
    /// The same blade appears twice.
    DuplicateBasisBlade,
    /// A runtime-metric value has no `p,q,r>` prefix.
    MissingMetric,
    /// The prefix's dimension `p + q + r` differs from the type's.
    InvalidMetricDimensions,
    /// The prefix is not a valid signature.
    InvalidMetric(MetricError),
    /// The prefix names a different signature than the compile-time one of the type.
    MismatchedMetric,
    /// A character that does not fit the grammar, at this byte offset.
    UnexpectedToken {
        /// Byte offset of the character in the input.
        position: usize,
    },
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidBasisBladeCoefficient => f.write_str("coefficient is not a number"),
            Self::InvalidBasisBladeIndex => f.write_str("blade `e` has no basis vectors"),
            Self::InvalidBasisVector => f.write_str("basis vector is not a digit"),
            Self::InvalidBasisVectorZero => f.write_str("basis vectors are numbered from 1"),
            Self::InvalidBasisVectorOrdering => f.write_str("basis vectors must be ascending"),
            Self::InvalidBasisVectorOutOfRange => f.write_str("basis vector above 6"),
            Self::DuplicateBasisVector => f.write_str("basis vector repeated within a blade"),
            Self::InvalidBasisBladeForAlgebra => f.write_str("blade is outside the algebra"),
            Self::DuplicateBasisBlade => f.write_str("blade appears more than once"),
            Self::MissingMetric => f.write_str("missing `p,q,r>` metric prefix"),
            Self::InvalidMetricDimensions => {
                f.write_str("metric dimension differs from the type's")
            }
            Self::InvalidMetric(_) => f.write_str("invalid metric prefix"),
            Self::MismatchedMetric => {
                f.write_str("metric prefix differs from the type's signature")
            }
            Self::UnexpectedToken { position } => {
                write!(f, "unexpected character at byte {position}")
            }
        }
    }
}

impl Error for ParseError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self {
            Self::InvalidMetric(inner) => Some(inner),
            _ => None,
        }
    }
}

impl From<MetricError> for ParseError {
    fn from(err: MetricError) -> Self {
        Self::InvalidMetric(err)
    }
}

/// Why a signature `p,q,r` is invalid.
///
/// # Examples
///
/// ```
/// use seagrass_math::prelude::*;
///
/// assert_eq!(Dyn4::new(3, 0, 0).unwrap_err(), MetricError::WrongDimension { expected: 4, found: 3 });
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MetricError {
    /// The signature is empty.
    MissingNumPositive,
    /// The count of positive basis vectors is not a number.
    InvalidNumPositive,
    /// The count of negative basis vectors is not a number.
    InvalidNumNegative,
    /// The count of null basis vectors is not a number.
    InvalidNumZero,
    /// More than three counts are given.
    InvalidMetricLength,
    /// `p + q + r` differs from the dimension of the bitset type.
    WrongDimension {
        /// The dimension of the bitset type.
        expected: u8,
        /// `p + q + r`.
        found: u8,
    },
}

impl fmt::Display for MetricError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingNumPositive => f.write_str("empty signature"),
            Self::InvalidNumPositive => f.write_str("count of positive vectors is not a number"),
            Self::InvalidNumNegative => f.write_str("count of negative vectors is not a number"),
            Self::InvalidNumZero => f.write_str("count of null vectors is not a number"),
            Self::InvalidMetricLength => f.write_str("signature has more than three counts"),
            Self::WrongDimension { expected, found } => {
                write!(f, "signature has dimension {found}, expected {expected}")
            }
        }
    }
}

impl Error for MetricError {}

/// Why a multivector could not be converted to a [`Typed`](crate::Typed) value of a given kind.
///
/// # Examples
///
/// ```
/// use seagrass_math::prelude::*;
///
/// let vector = Multivector::<f64, Pga3>::from_str("e1").unwrap();
/// assert_eq!(Motor::try_from(vector).unwrap_err(), KindError::WrongParity);
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum KindError {
    /// The value stores a blade whose grade parity the kind excludes.
    WrongParity,
}

impl fmt::Display for KindError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str("multivector stores blades of the wrong grade parity for this kind")
    }
}

impl Error for KindError {}

/// Why blades and coefficients could not be assembled into a multivector.
///
/// # Examples
///
/// ```
/// use seagrass_math::prelude::*;
///
/// let err = Multivector::<f64, Pga3>::from_blades(Pga3::default(), [(BladeIndex::E12345, 1.0)]);
/// assert_eq!(err.unwrap_err(), BladeError::OutOfRange { index: BladeIndex::E12345 });
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum BladeError {
    /// The blade uses a basis vector beyond the algebra's dimension.
    OutOfRange {
        /// The offending blade.
        index: crate::BladeIndex,
    },
    /// A value stores a blade that a [`Batch`](crate::Batch) does not.
    NotStored {
        /// The blade the batch does not store.
        index: crate::BladeIndex,
    },
    /// The number of coefficients differs from the number of blades in the bitset.
    CoefficientCount {
        /// Blades in the bitset.
        expected: usize,
        /// Coefficients given.
        found: usize,
    },
}

impl fmt::Display for BladeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::OutOfRange { index } => write!(f, "blade {index} is outside the algebra"),
            Self::NotStored { index } => write!(f, "blade {index} is not stored by the batch"),
            Self::CoefficientCount { expected, found } => {
                write!(f, "{found} coefficients for {expected} blades")
            }
        }
    }
}

impl Error for BladeError {}
