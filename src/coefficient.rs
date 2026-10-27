//! The scalar type of multivector coefficients.

use std::fmt::{Debug, Display};
use std::ops::{AddAssign, DivAssign, MulAssign, SubAssign};
use std::str::FromStr;

use num_traits::real::Real;

/// A real number type usable as a multivector coefficient, in practice `f32` or `f64`.
///
/// Implemented for every type with the listed capabilities, so it never needs an explicit impl.
///
/// # Examples
///
/// ```
/// use seagrass_math::Coefficient;
///
/// fn halve<T: Coefficient>(x: T) -> T {
///     x / (T::one() + T::one())
/// }
/// assert_eq!(halve(3.0_f32), 1.5);
/// ```
pub trait Coefficient:
    Real
    + Debug
    + Display
    + FromStr
    + AddAssign
    + SubAssign
    + MulAssign
    + DivAssign
    + Send
    + Sync
    + 'static
{
}

impl<T> Coefficient for T where
    T: Real
        + Debug
        + Display
        + FromStr
        + AddAssign
        + SubAssign
        + MulAssign
        + DivAssign
        + Send
        + Sync
        + 'static
{
}

/// Converts a constant such as `0.5` or a factorial to `T`.
pub(crate) fn constant<T: Coefficient>(x: f64) -> T {
    T::from(x).expect("f32 and f64 represent every constant used by the crate")
}
