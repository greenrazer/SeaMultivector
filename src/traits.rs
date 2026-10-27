//! The products as named operations, for generic code. Each operator forwards to one of these.
//!
//! | Operator | Trait |
//! |---|---|
//! | `a * b` | [`GeometricProduct`] |
//! | `a ^ b` | [`OuterProduct`] |
//! | `a & b` | [`RegressiveProduct`] |
//! | `a \| b` | [`InnerProduct`] |
//! | `a << b` | [`LeftContraction`] |
//! | `a >> b` | [`RightContraction`] |
//! | `m / x` | [`SandwichProduct`] (on [`Typed`](crate::Typed) values) |

/// The geometric product, `a * b`.
///
/// # Examples
///
/// ```
/// use seagrass_math::prelude::*;
///
/// let e1 = Multivector::<f64, Euclid3>::basis(Euclid3::default(), BladeIndex::E1);
/// let e2 = Multivector::<f64, Euclid3>::basis(Euclid3::default(), BladeIndex::E2);
/// assert_eq!(e1.geometric_product(&e2).to_string(), "3,0,0>e12");
/// ```
pub trait GeometricProduct<Rhs = Self> {
    /// The type of the product.
    type Output;

    /// The geometric product of `self` and `rhs`.
    fn geometric_product(&self, rhs: &Rhs) -> Self::Output;
}

/// The outer (wedge) product, `a ^ b`; the meet in projective geometric algebra.
///
/// # Examples
///
/// ```
/// use seagrass_math::prelude::*;
///
/// let mv = |s: &str| Multivector::<f64, Euclid3>::from_str(s).unwrap();
/// assert_eq!((mv("e1") ^ mv("e2")).to_string(), "3,0,0>e12");
/// assert_eq!((mv("e1") ^ mv("e12")).to_string(), "3,0,0>");
/// ```
pub trait OuterProduct<Rhs = Self> {
    /// The type of the product.
    type Output;

    /// The outer product of `self` and `rhs`.
    fn outer_product(&self, rhs: &Rhs) -> Self::Output;
}

/// The inner product, `a | b`: the sum of the left and right contractions, counting the scalar
/// overlap once.
///
/// # Examples
///
/// ```
/// use seagrass_math::prelude::*;
///
/// let mv = |s: &str| Multivector::<f64, Euclid3>::from_str(s).unwrap();
/// assert_eq!((mv("e1") | mv("e12")).get(BladeIndex::E2), 1.0);
/// assert_eq!((mv("e12") | mv("e1")).get(BladeIndex::E2), -1.0);
/// ```
pub trait InnerProduct<Rhs = Self> {
    /// The type of the product.
    type Output;

    /// The inner product of `self` and `rhs`.
    fn inner_product(&self, rhs: &Rhs) -> Self::Output;
}

/// The left contraction, `a << b`.
///
/// # Examples
///
/// ```
/// use seagrass_math::prelude::*;
///
/// let mv = |s: &str| Multivector::<f64, Euclid3>::from_str(s).unwrap();
/// assert_eq!((mv("e1") << mv("e12")).to_string(), "3,0,0>e2");
/// assert!((mv("e12") << mv("e1")).is_zero());
/// ```
pub trait LeftContraction<Rhs = Self> {
    /// The type of the product.
    type Output;

    /// The left contraction of `rhs` by `self`.
    fn left_contraction(&self, rhs: &Rhs) -> Self::Output;
}

/// The right contraction, `a >> b`.
///
/// # Examples
///
/// ```
/// use seagrass_math::prelude::*;
///
/// let mv = |s: &str| Multivector::<f64, Euclid3>::from_str(s).unwrap();
/// assert_eq!((mv("e12") >> mv("e2")).to_string(), "3,0,0>e1");
/// ```
pub trait RightContraction<Rhs = Self> {
    /// The type of the product.
    type Output;

    /// The right contraction of `self` by `rhs`.
    fn right_contraction(&self, rhs: &Rhs) -> Self::Output;
}

/// The regressive (vee) product, `a & b`: the outer product of the duals, undualized. The join
/// in projective geometric algebra.
///
/// # Examples
///
/// ```
/// use seagrass_math::prelude::*;
///
/// let mv = |s: &str| Multivector::<f64, Euclid3>::from_str(s).unwrap();
/// // Two planes through the origin meet in a line: e12 and e23 share the e2 direction.
/// assert_eq!((mv("e12") & mv("e23")).to_string(), "3,0,0>-e2");
/// ```
pub trait RegressiveProduct<Rhs = Self> {
    /// The type of the product.
    type Output;

    /// The regressive product of `self` and `rhs`.
    fn regressive_product(&self, rhs: &Rhs) -> Self::Output;
}

/// The scalar product: the scalar part of the geometric product.
///
/// # Examples
///
/// ```
/// use seagrass_math::prelude::*;
///
/// let mv = |s: &str| Multivector::<f64, Euclid3>::from_str(s).unwrap();
/// assert_eq!(mv("1 + 2e1 + e12").scalar_product(&mv("3 + 4e1 + e12")), 10.0);
/// ```
pub trait ScalarProduct<Rhs = Self> {
    /// The scalar type.
    type Output;

    /// The scalar part of `self * rhs`.
    fn scalar_product(&self, rhs: &Rhs) -> Self::Output;
}

/// The sandwich product `m x reverse(m)`, which applies a versor m to x.
///
/// # Examples
///
/// ```
/// use seagrass_math::prelude::*;
///
/// let mv = |s: &str| Multivector::<f64, Euclid3>::from_str(s).unwrap();
/// // e1 e2 e1 = -e2; the reflection in the plane orthogonal to e1 is the negated sandwich.
/// assert_eq!(mv("e1").sandwich_product(&mv("e2")).to_string(), "3,0,0>-e2");
/// ```
pub trait SandwichProduct<Rhs = Self> {
    /// The type of the transformed value.
    type Output;

    /// `self * rhs * reverse(self)`, computed without building the reverse.
    fn sandwich_product(&self, rhs: &Rhs) -> Self::Output;
}
