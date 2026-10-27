//! The multivector type: a bitset of the stored blades plus their packed coefficients.

use crate::blade::{BladeIndex, BladeSet};
use crate::coefficient::{Coefficient, constant};
use crate::error::BladeError;
use crate::metric::{Metric, ProductConstants};
use crate::plan::{Coefficients, Op, with_product_plan};

mod arithmetic;
mod display;
mod iter;
mod ops;
mod parse;
mod products;
#[cfg(test)]
mod reference;
mod unary;

pub use iter::{Iter, IterMut};

/// A multivector of the algebra with metric `M`, storing only the blades it holds.
///
/// A value is a bitset of its stored blades ([`Metric::Bits`]) plus one coefficient per stored
/// blade, packed in blade-index order; a blade's slot is the number of stored blades below it.
/// Products work out which blades their result can hold from the operands' bitsets before they
/// touch a coefficient, so the work grows with the stored blades rather than with 2^n.
///
/// Equality is structural: the stored blades are part of the value, so `1 + 0e12 != 1`. Use
/// [`remove_zeros`](Self::remove_zeros) to drop stored zeros, or compare coefficients with a
/// tolerance.
///
/// # Examples
///
/// ```
/// use seagrass_math::prelude::*;
///
/// // A runtime metric: the signature comes from the string.
/// let a = Multivector::<f64, Dyn6>::from_str("6>e1 + 2e2").unwrap();
/// let b = Multivector::<f64, Dyn6>::from_str("6>3e5 + e6").unwrap();
/// assert_eq!((&a * &b).to_string(), "6,0,0>3e15 + 6e25 + e16 + 2e26");
///
/// // A compile-time metric needs no prefix.
/// let e1 = Multivector::<f64, Pga3>::basis(Pga3::default(), BladeIndex::E1);
/// assert_eq!((&e1 * &e1).scalar(), 1.0);
/// ```
#[derive(Debug, PartialEq)]
pub struct Multivector<T, M: Metric> {
    metric: M,
    blades: M::Bits,
    coefficients: Coefficients<T>,
}

impl<T: Coefficient, M: Metric> Multivector<T, M> {
    /// Assembles a value whose coefficient count is known to match `blades`.
    #[inline]
    pub(crate) fn from_raw(metric: M, blades: M::Bits, coefficients: Coefficients<T>) -> Self {
        debug_assert_eq!(blades.count_ones() as usize, coefficients.len());
        Self {
            metric,
            blades,
            coefficients,
        }
    }

    /// The zero multivector, which stores no blades.
    ///
    /// # Examples
    ///
    /// ```
    /// use seagrass_math::prelude::*;
    ///
    /// let zero = Multivector::<f64, Dyn4>::zero(Dyn4::new(1, 3, 0).unwrap());
    /// assert!(zero.is_zero());
    /// assert_eq!(zero.to_string(), "1,3,0>");
    /// ```
    pub fn zero(metric: M) -> Self {
        Self::from_raw(metric, M::Bits::EMPTY, Coefficients::new())
    }

    /// The scalar `value`.
    pub fn from_scalar(metric: M, value: T) -> Self {
        Self::from_raw(metric, M::Bits::SCALAR, Coefficients::from_slice(&[value]))
    }

    /// The basis blade `index` with coefficient 1.
    ///
    /// # Panics
    ///
    /// Panics if `index` uses a basis vector beyond the metric's dimension.
    ///
    /// # Examples
    ///
    /// ```
    /// use seagrass_math::prelude::*;
    ///
    /// // In spacetime algebra (1, 3) the time vector e1 squares to +1 and e2 to -1.
    /// let sta = Sta::default();
    /// let e1 = Multivector::<f64, Sta>::basis(sta, BladeIndex::E1);
    /// let e2 = Multivector::<f64, Sta>::basis(sta, BladeIndex::E2);
    /// assert_eq!(((&e1 * &e1).scalar(), (&e2 * &e2).scalar()), (1.0, -1.0));
    /// ```
    pub fn basis(metric: M, index: BladeIndex) -> Self {
        assert!(
            u32::from(index.bits()) < M::Bits::BLADES,
            "blade {index} is outside the algebra"
        );
        Self::from_raw(
            metric,
            M::Bits::from_index(index),
            Coefficients::from_slice(&[T::one()]),
        )
    }

    /// The sum of the given blades; a blade listed twice gets the sum of its coefficients.
    ///
    /// Building values from numbers this way skips the parser, which costs as much as dozens of
    /// products.
    ///
    /// # Errors
    ///
    /// [`BladeError::OutOfRange`] if a blade uses a basis vector beyond the metric's dimension.
    ///
    /// # Examples
    ///
    /// ```
    /// use seagrass_math::prelude::*;
    ///
    /// let line = Multivector::<f64, Pga3>::from_blades(
    ///     Pga3::default(),
    ///     [(BladeIndex::E12, 2.0), (BladeIndex::E34, -1.0)],
    /// )
    /// .unwrap();
    /// assert_eq!(line.to_string(), "3,0,1>2e12 - e34");
    /// ```
    pub fn from_blades(
        metric: M,
        blades: impl IntoIterator<Item = (BladeIndex, T)>,
    ) -> Result<Self, BladeError> {
        let mut dense = [T::zero(); 64];
        let mut stored = M::Bits::EMPTY;
        for (index, value) in blades {
            if u32::from(index.bits()) >= M::Bits::BLADES {
                return Err(BladeError::OutOfRange { index });
            }
            stored |= M::Bits::from_index(index);
            dense[index.bits() as usize] += value;
        }
        let coefficients = stored
            .iter()
            .map(|index| dense[index.bits() as usize])
            .collect();
        Ok(Self::from_raw(metric, stored, coefficients))
    }

    /// The value with the given stored blades and their coefficients in blade-index order.
    ///
    /// # Errors
    ///
    /// [`BladeError::CoefficientCount`] if the number of coefficients differs from the number of
    /// blades in `blades`.
    ///
    /// # Examples
    ///
    /// ```
    /// use seagrass_math::prelude::*;
    ///
    /// // Blades 1 (scalar) and 3 (e12) of a 3D algebra.
    /// let rotor = Multivector::<f64, Euclid3>::from_parts(Euclid3::default(), 0b1001, vec![0.6, 0.8]);
    /// assert_eq!(rotor.unwrap().to_string(), "3,0,0>0.6 + 0.8e12");
    /// let wrong = Multivector::<f64, Euclid3>::from_parts(Euclid3::default(), 0b1001, vec![1.0]);
    /// assert_eq!(wrong.unwrap_err(), BladeError::CoefficientCount { expected: 2, found: 1 });
    /// ```
    pub fn from_parts(
        metric: M,
        blades: M::Bits,
        coefficients: Vec<T>,
    ) -> Result<Self, BladeError> {
        let expected = blades.count_ones() as usize;
        if coefficients.len() != expected {
            return Err(BladeError::CoefficientCount {
                expected,
                found: coefficients.len(),
            });
        }
        Ok(Self::from_raw(
            metric,
            blades,
            Coefficients::from_vec(coefficients),
        ))
    }

    /// The metric.
    #[inline]
    pub fn metric(&self) -> M {
        self.metric
    }

    /// The set of stored blades.
    #[inline]
    pub fn blades(&self) -> M::Bits {
        self.blades
    }

    /// The coefficients of the stored blades, in blade-index order.
    #[inline]
    pub fn coefficients(&self) -> &[T] {
        &self.coefficients
    }

    /// The coefficients of the stored blades, mutably. The stored blades cannot change.
    #[inline]
    pub fn coefficients_mut(&mut self) -> &mut [T] {
        &mut self.coefficients
    }

    /// The metric, the stored blades and their coefficients.
    pub fn into_parts(self) -> (M, M::Bits, Vec<T>) {
        (self.metric, self.blades, self.coefficients.into_vec())
    }

    /// The coefficient of blade `index`, zero when it is not stored.
    ///
    /// # Examples
    ///
    /// ```
    /// use seagrass_math::prelude::*;
    ///
    /// let a = Multivector::<f64, Pga3>::from_str("2 - 3e14").unwrap();
    /// assert_eq!((a.get(BladeIndex::E14), a.get(BladeIndex::E23), a.scalar()), (-3.0, 0.0, 2.0));
    /// ```
    #[inline]
    pub fn get(&self, index: BladeIndex) -> T {
        if !self.blades.has(index) {
            return T::zero();
        }
        self.coefficients[self.blades.count_below(index) as usize]
    }

    /// The scalar part.
    #[inline]
    pub fn scalar(&self) -> T {
        self.get(BladeIndex::SCALAR)
    }

    /// The coefficient of the pseudoscalar e1..en.
    #[inline]
    pub fn pseudoscalar(&self) -> T {
        self.get(BladeIndex::first_vectors(M::Bits::DIM))
    }

    /// Whether no blade is stored. A value that stores zeros is not zero in this sense.
    #[inline]
    pub fn is_zero(&self) -> bool {
        self.blades.is_empty()
    }

    /// The grade-`k` part: the stored blades with `k` basis vectors.
    ///
    /// # Panics
    ///
    /// Panics if `k` is above 6.
    ///
    /// # Examples
    ///
    /// ```
    /// use seagrass_math::prelude::*;
    ///
    /// let a = Multivector::<f64, Euclid3>::from_str("1 + 2e1 + 3e12 - e23 + e123").unwrap();
    /// assert_eq!(a.grade(2).to_string(), "3,0,0>3e12 - e23");
    /// ```
    pub fn grade(&self, k: u32) -> Self {
        self.masked(M::Bits::GRADE[k as usize])
    }

    /// The stored blades that are in `mask`, with their coefficients.
    pub(crate) fn masked(&self, mask: M::Bits) -> Self {
        let coefficients = self
            .iter()
            .filter(|(index, _)| mask.has(*index))
            .map(|(_, c)| *c)
            .collect();
        Self::from_raw(self.metric, self.blades & mask, coefficients)
    }

    /// Drops the stored blades whose coefficient is zero.
    ///
    /// # Examples
    ///
    /// ```
    /// use seagrass_math::prelude::*;
    ///
    /// let mut a = Multivector::<f64, Euclid3>::from_str("1 + 0e12").unwrap();
    /// assert_ne!(a, Multivector::from(1.0)); // equality compares stored blades
    /// a.remove_zeros();
    /// assert_eq!(a, Multivector::from(1.0));
    /// ```
    pub fn remove_zeros(&mut self) {
        let mut blades = self.blades;
        for (index, c) in self.blades.iter().zip(&self.coefficients) {
            if c.is_zero() {
                blades ^= M::Bits::from_index(index);
            }
        }
        self.blades = blades;
        self.coefficients.retain(|c| !c.is_zero());
    }

    /// Drops the stored blades whose coefficient is smaller than `epsilon` in magnitude.
    pub fn remove_approx_zeros(&mut self, epsilon: T) {
        let mut blades = self.blades;
        for (index, c) in self.blades.iter().zip(&self.coefficients) {
            if c.abs() < epsilon {
                blades ^= M::Bits::from_index(index);
            }
        }
        self.blades = blades;
        self.coefficients.retain(|c| c.abs() >= epsilon);
    }

    /// Rounds every coefficient to `decimals` decimal places.
    pub fn round_to_decimals(&mut self, decimals: i32) {
        let scale = constant::<T>(10.0).powi(decimals);
        for c in &mut self.coefficients {
            *c = (*c * scale).round() / scale;
        }
    }

    /// Drops coefficients below 1e-15 and rounds the rest to 14 decimals, which cleans up the
    /// rounding noise of a chain of products.
    pub fn simplify(&mut self) {
        self.simplify_epsilon(constant(1e-15));
    }

    /// Drops coefficients below `epsilon` and rounds the rest to 14 decimals.
    pub fn simplify_epsilon(&mut self, epsilon: T) {
        self.remove_approx_zeros(epsilon);
        self.round_to_decimals(14);
    }

    /// Panics unless both operands share a metric; only runtime metrics can differ.
    #[inline]
    pub(crate) fn assert_same_metric(&self, rhs: &Self) {
        assert_eq!(self.metric, rhs.metric, "operands have different metrics");
    }

    /// A bilinear product through its cached [plan](crate::plan::ProductPlan).
    ///
    /// # Panics
    ///
    /// Panics if the operands have different runtime metrics.
    #[inline]
    pub(crate) fn planned_product(&self, op: Op, rhs: &Self) -> Self {
        self.assert_same_metric(rhs);
        with_product_plan(op, &self.metric, self.blades, rhs.blades, |plan| {
            Self::from_raw(
                self.metric,
                plan.output(),
                plan.run(&self.coefficients, &rhs.coefficients),
            )
        })
    }

    #[inline]
    pub(crate) fn constants(&self) -> ProductConstants {
        ProductConstants::of(&self.metric)
    }

    /// Whether the value is a scalar: only the scalar blade, or other blades that are all zero.
    pub(crate) fn is_scalar(&self) -> bool {
        if self.blades == M::Bits::SCALAR {
            true
        } else if self.blades.has(BladeIndex::SCALAR) {
            self.coefficients[1..].iter().all(|c| c.is_zero())
        } else {
            false
        }
    }

    /// Negates the stored blades that are in `mask`.
    pub(crate) fn negate_blades(&mut self, mask: M::Bits) {
        let coefficients = self.coefficients.as_mut_slice();
        for index in (self.blades & mask).iter() {
            let slot = self.blades.count_below(index) as usize;
            coefficients[slot] = -coefficients[slot];
        }
    }

    /// A copy with the stored blades in `mask` negated.
    pub(crate) fn negated_blades(&self, mask: M::Bits) -> Self {
        let mut out = self.clone();
        out.negate_blades(mask);
        out
    }
}

impl<T: Coefficient, M: Metric> Clone for Multivector<T, M> {
    /// Copies the coefficients with one `memcpy`, also when they do not fit inline.
    fn clone(&self) -> Self {
        Self::from_raw(
            self.metric,
            self.blades,
            Coefficients::from_slice(&self.coefficients),
        )
    }
}

impl<T: Coefficient, M: Metric + Default> Default for Multivector<T, M> {
    /// The zero multivector.
    fn default() -> Self {
        Self::zero(M::default())
    }
}

impl<T: Coefficient, M: Metric + Default> From<T> for Multivector<T, M> {
    /// The scalar `value` in a compile-time metric.
    fn from(value: T) -> Self {
        Self::from_scalar(M::default(), value)
    }
}

impl<T, M: Metric> AsRef<Multivector<T, M>> for Multivector<T, M> {
    fn as_ref(&self) -> &Self {
        self
    }
}
