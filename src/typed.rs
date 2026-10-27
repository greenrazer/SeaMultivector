//! Multivectors with a known grade parity: motors (even), flectors (odd) and versors (either).

use std::fmt;
use std::marker::PhantomData;

use crate::blade::{BladeIndex, BladeSet};
use crate::coefficient::Coefficient;
use crate::error::KindError;
use crate::metric::Metric;
use crate::multivector::Multivector;
use crate::traits::{SandwichProduct, ScalarProduct};

mod ops;

mod sealed {
    pub trait Sealed {}
}

/// A grade parity that a [`Typed`] value guarantees: [`Even`], [`Odd`] or [`Any`].
///
/// The kind of a product follows from the kinds of its operands (see [`KindMul`]), so a motor
/// times a motor is a motor and a motor applied to a flector is a flector, checked at compile
/// time. The trait cannot be implemented outside the crate.
pub trait Kind:
    Copy + Default + fmt::Debug + PartialEq + Eq + Send + Sync + 'static + sealed::Sealed
{
    /// Whether a value storing `blades` has this kind.
    fn admits<B: BladeSet>(blades: B) -> bool;
}

/// Only even grades: rotors, motors, and lines in projective geometric algebra.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub struct Even;

/// Only odd grades: reflections, flectors, and points and planes in projective geometric algebra.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub struct Odd;

/// No parity guarantee: a versor whose parity is not tracked, such as a general product.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub struct Any;

impl sealed::Sealed for Even {}
impl sealed::Sealed for Odd {}
impl sealed::Sealed for Any {}

impl Kind for Even {
    fn admits<B: BladeSet>(blades: B) -> bool {
        blades & B::ODD == B::EMPTY
    }
}

impl Kind for Odd {
    fn admits<B: BladeSet>(blades: B) -> bool {
        blades & B::EVEN == B::EMPTY
    }
}

impl Kind for Any {
    fn admits<B: BladeSet>(_: B) -> bool {
        true
    }
}

/// The kind of a product: parities add, so even times odd is odd and odd times odd is even.
///
/// This holds for the geometric, outer and inner products and the contractions, because their
/// grades are `r + s` or `|r - s|`, which have the parity of `r + s`. [`Any`] absorbs everything.
///
/// # Examples
///
/// ```
/// use seagrass_math::prelude::*;
///
/// let point: Flector<f64, Pga3> = Multivector::new_point(1.0, 0.0, 0.0);
/// let plane: Flector<f64, Pga3> = Multivector::new_plane(0.0, 0.0, 1.0, 0.0);
/// // A reflection times a reflection is a motion: odd times odd is even.
/// let _: Motor<f64, Pga3> = &point * &plane;
/// ```
pub trait KindMul<Rhs: Kind>: Kind {
    /// The kind of the product.
    type Output: Kind;
}

macro_rules! kind_table {
    ($($lhs:ty, $rhs:ty => $out:ty;)*) => {$(
        impl KindMul<$rhs> for $lhs {
            type Output = $out;
        }
    )*};
}

kind_table! {
    Even, Even => Even;
    Even, Odd => Odd;
    Odd, Even => Odd;
    Odd, Odd => Even;
    Even, Any => Any;
    Odd, Any => Any;
    Any, Even => Any;
    Any, Odd => Any;
    Any, Any => Any;
}

/// The kind of a product of kinds `K1` and `K2`.
pub(crate) type Product<K1, K2> = <K1 as KindMul<K2>>::Output;

/// The parity the dimension of metric `M` adds to duals and regressive products: the dual of a
/// grade-k blade has grade n - k.
pub(crate) type DimensionParity<M> = <<M as Metric>::Bits as BladeSet>::Parity;

/// A multivector whose stored blades all have the grade parity `K`.
///
/// Products of typed values are typed by [`KindMul`], and `/` applies one as a versor with the
/// sandwich product, keeping the kind of the value it is applied to. Mixing a typed value with a
/// plain [`Multivector`] gives a plain multivector. [`Motor`], [`Flector`] and [`Versor`] name the
/// three kinds.
///
/// # Examples
///
/// ```
/// use seagrass_math::prelude::*;
///
/// let rotate = Motor::<f64, Pga3>::new_rotation_motor(0.0, 0.0, std::f64::consts::FRAC_PI_2);
/// let point: Flector<f64, Pga3> = Multivector::<f64, Pga3>::new_point(1.0, 0.0, 0.0);
/// let moved: Flector<f64, Pga3> = &rotate / &point;
/// assert!((moved.point_to_euclidean() - Vector3::new(0.0, 1.0, 0.0)).norm() < 1e-12);
/// ```
#[derive(Debug, PartialEq)]
pub struct Typed<T, M: Metric, K: Kind> {
    inner: Multivector<T, M>,
    kind: PhantomData<K>,
}

/// An even multivector: a rotor, a motor, or a line in projective geometric algebra.
pub type Motor<T, M> = Typed<T, M, Even>;

/// An odd multivector: a flector, or a point or plane in projective geometric algebra.
pub type Flector<T, M> = Typed<T, M, Odd>;

/// A multivector whose grade parity is not tracked.
pub type Versor<T, M> = Typed<T, M, Any>;

impl<T: Coefficient, M: Metric, K: Kind> Typed<T, M, K> {
    /// Wraps a value known to have kind `K`.
    #[inline]
    pub(crate) fn wrap(inner: Multivector<T, M>) -> Self {
        debug_assert!(
            K::admits(inner.blades()),
            "value does not have the kind it is wrapped as"
        );
        Self {
            inner,
            kind: PhantomData,
        }
    }

    /// The value as a typed value of kind `K`.
    ///
    /// # Errors
    ///
    /// [`KindError::WrongParity`] if the value stores a blade of the other parity, even with a
    /// zero coefficient.
    pub fn new(inner: Multivector<T, M>) -> Result<Self, KindError> {
        if K::admits(inner.blades()) {
            Ok(Self {
                inner,
                kind: PhantomData,
            })
        } else {
            Err(KindError::WrongParity)
        }
    }

    /// The untyped value.
    #[inline]
    pub fn as_multivector(&self) -> &Multivector<T, M> {
        &self.inner
    }

    /// The untyped value, by value.
    #[inline]
    pub fn into_inner(self) -> Multivector<T, M> {
        self.inner
    }

    /// The metric.
    #[inline]
    pub fn metric(&self) -> M {
        self.inner.metric()
    }

    /// The set of stored blades.
    #[inline]
    pub fn blades(&self) -> M::Bits {
        self.inner.blades()
    }

    /// The coefficients of the stored blades, in blade-index order.
    #[inline]
    pub fn coefficients(&self) -> &[T] {
        self.inner.coefficients()
    }

    /// The coefficients of the stored blades, mutably. The stored blades, and so the kind, cannot
    /// change.
    #[inline]
    pub fn coefficients_mut(&mut self) -> &mut [T] {
        self.inner.coefficients_mut()
    }

    /// The coefficient of blade `index`, zero when it is not stored.
    #[inline]
    pub fn get(&self, index: BladeIndex) -> T {
        self.inner.get(index)
    }

    /// The scalar part.
    #[inline]
    pub fn scalar(&self) -> T {
        self.inner.scalar()
    }

    /// The coefficient of the pseudoscalar.
    #[inline]
    pub fn pseudoscalar(&self) -> T {
        self.inner.pseudoscalar()
    }

    /// Whether no blade is stored.
    #[inline]
    pub fn is_zero(&self) -> bool {
        self.inner.is_zero()
    }

    /// The grade-`k` part, as a plain multivector.
    ///
    /// # Panics
    ///
    /// Panics if `k` is above 6.
    pub fn grade(&self, k: u32) -> Multivector<T, M> {
        self.inner.grade(k)
    }

    /// The stored blades and their coefficients, in blade-index order.
    pub fn iter(&self) -> impl ExactSizeIterator<Item = (BladeIndex, &T)> {
        self.inner.iter()
    }

    /// The stored blades and their coefficients, mutably.
    pub fn iter_mut(&mut self) -> impl ExactSizeIterator<Item = (BladeIndex, &mut T)> {
        self.inner.iter_mut()
    }

    /// The stored blades of grade `k` and their coefficients.
    pub fn grade_iter(&self, k: u32) -> impl Iterator<Item = (BladeIndex, &T)> {
        self.inner.grade_iter(k)
    }

    /// The stored blades of grade `k` and their coefficients, mutably.
    pub fn grade_iter_mut(&mut self, k: u32) -> impl Iterator<Item = (BladeIndex, &mut T)> {
        self.inner.grade_iter_mut(k)
    }

    /// The norm; see [`Multivector::norm`].
    pub fn norm(&self) -> T {
        self.inner.norm()
    }

    /// The inverse `reverse(self) / (self * reverse(self))`, exact for any versor (whose product
    /// with its reverse is a scalar). For a unit versor it equals the reverse.
    ///
    /// # Examples
    ///
    /// ```
    /// use seagrass_math::prelude::*;
    ///
    /// // A rotor scaled by 2 is still invertible; its reverse alone would not undo it.
    /// let rotor = Motor::<f64, Euclid3>::new(Multivector::from_str("1.2 + 1.6e12").unwrap()).unwrap();
    /// let mut identity = &rotor * &rotor.inverse();
    /// identity.simplify();
    /// assert_eq!(identity.to_string(), "3,0,0>1");
    /// ```
    #[must_use]
    pub fn inverse(&self) -> Self {
        let reversed = self.inner.reversed();
        let scale = self.inner.scalar_product(&reversed);
        Self::wrap(reversed / scale)
    }

    /// Applies `self` as a versor to `x` with the sandwich product, like `self / x`.
    ///
    /// # Examples
    ///
    /// ```
    /// use seagrass_math::prelude::*;
    ///
    /// let shift = Motor::<f64, Pga3>::new_translation_motor(0.0, 0.0, 2.0);
    /// let origin = Multivector::<f64, Pga3>::new_point(0.0, 0.0, 0.0);
    /// assert_eq!(shift.apply(&origin), &shift / &origin);
    /// assert_eq!(shift.apply(&origin).point_to_euclidean(), Vector3::new(0.0, 0.0, 2.0));
    /// ```
    pub fn apply<X>(&self, x: &X) -> <Self as SandwichProduct<X>>::Output
    where
        Self: SandwichProduct<X>,
    {
        self.sandwich_product(x)
    }
}

/// Methods that keep the kind: they only change signs or magnitudes of stored blades.
macro_rules! kind_preserving {
    ($($(#[$doc:meta])* $returning:ident / $in_place:ident;)*) => {
        impl<T: Coefficient, M: Metric, K: Kind> Typed<T, M, K> {
            $(
                $(#[$doc])*
                #[must_use]
                pub fn $returning(&self) -> Self {
                    Self::wrap(self.inner.$returning())
                }

                $(#[$doc])*
                ///
                /// In place.
                pub fn $in_place(&mut self) {
                    self.inner.$in_place();
                }
            )*

            /// Drops the stored blades whose coefficient is zero.
            pub fn remove_zeros(&mut self) {
                self.inner.remove_zeros();
            }

            /// Drops the stored blades whose coefficient is smaller than `epsilon` in magnitude.
            pub fn remove_approx_zeros(&mut self, epsilon: T) {
                self.inner.remove_approx_zeros(epsilon);
            }

            /// Rounds every coefficient to `decimals` decimal places.
            pub fn round_to_decimals(&mut self, decimals: i32) {
                self.inner.round_to_decimals(decimals);
            }

            /// Drops coefficients below 1e-15 and rounds the rest to 14 decimals.
            pub fn simplify(&mut self) {
                self.inner.simplify();
            }

            /// Drops coefficients below `epsilon` and rounds the rest to 14 decimals.
            pub fn simplify_epsilon(&mut self, epsilon: T) {
                self.inner.simplify_epsilon(epsilon);
            }
        }
    };
}

kind_preserving! {
    /// The reverse; see [`Multivector::reversed`].
    reversed / reverse;
    /// The grade involution.
    involuted / involute;
    /// The Clifford conjugate.
    conjugated / conjugate;
    /// The negation.
    negated / negate;
    /// The value divided by its norm.
    normalized / normalize;
}

/// Duals and complements map grade k to n - k, which flips the parity in odd dimensions.
macro_rules! dimension_dual {
    ($($(#[$doc:meta])* $name:ident;)*) => {
        impl<T: Coefficient, M: Metric, K: KindMul<DimensionParity<M>>> Typed<T, M, K> {
            $(
                $(#[$doc])*
                #[must_use]
                pub fn $name(&self) -> Typed<T, M, Product<K, DimensionParity<M>>> {
                    Typed::wrap(self.inner.$name())
                }
            )*
        }
    };
}

dimension_dual! {
    /// The dual; see [`Multivector::dual`].
    dual;
    /// The undual; see [`Multivector::undual`].
    undual;
    /// The metric-free dual; see [`Multivector::poincare_dual`].
    poincare_dual;
    /// The right complement; see [`Multivector::right_complemented`].
    right_complemented;
    /// The left complement; see [`Multivector::left_complemented`].
    left_complemented;
}

impl<T: Coefficient, M: Metric> Typed<T, M, Even> {
    /// The exponential, which is even for an even argument: the exponential of a bivector is a
    /// rotor or motor. See [`Multivector::exp`].
    #[must_use]
    pub fn exp(&self) -> Self {
        Self::wrap(self.inner.exp())
    }
}

impl<T: Coefficient, M: Metric, K: Kind> Clone for Typed<T, M, K> {
    fn clone(&self) -> Self {
        Self::wrap(self.inner.clone())
    }
}

impl<T: Coefficient, M: Metric, K: Kind> TryFrom<Multivector<T, M>> for Typed<T, M, K> {
    type Error = KindError;

    /// See [`Typed::new`].
    fn try_from(inner: Multivector<T, M>) -> Result<Self, KindError> {
        Self::new(inner)
    }
}

impl<T: Coefficient, M: Metric, K: Kind> From<Typed<T, M, K>> for Multivector<T, M> {
    fn from(value: Typed<T, M, K>) -> Self {
        value.inner
    }
}

impl<T, M: Metric, K: Kind> AsRef<Multivector<T, M>> for Typed<T, M, K> {
    fn as_ref(&self) -> &Multivector<T, M> {
        &self.inner
    }
}

impl<T: Coefficient, M: Metric + Default, K: Kind> Default for Typed<T, M, K> {
    /// Zero, which stores no blades and so has every kind.
    fn default() -> Self {
        Self::wrap(Multivector::default())
    }
}

impl<T: Coefficient, M: Metric + Default> From<T> for Typed<T, M, Even> {
    /// The scalar `value`.
    fn from(value: T) -> Self {
        Self::wrap(Multivector::from(value))
    }
}

impl<T: Coefficient, M: Metric, K: Kind> fmt::Display for Typed<T, M, K> {
    /// The same syntax as [`Multivector`]'s `Display`.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt::Display::fmt(&self.inner, f)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::metric::Pga3;

    #[test]
    fn kinds_admit_only_their_parity() {
        assert!(Even::admits(0b1001_0110_0110_1001_u16));
        assert!(!Even::admits(0b10_u16));
        assert!(Odd::admits(0b10_u16));
        assert!(!Odd::admits(1_u16));
        assert!(Any::admits(u16::MAX));
    }

    #[test]
    fn rejects_values_of_the_wrong_parity() {
        let vector = Multivector::<f64, Pga3>::basis(Pga3::default(), BladeIndex::E1);
        assert_eq!(Motor::try_from(vector.clone()), Err(KindError::WrongParity));
        assert!(Flector::try_from(vector).is_ok());
    }
}
