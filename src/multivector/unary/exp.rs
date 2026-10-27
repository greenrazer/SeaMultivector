use super::*;

use crate::traits::GeometricProduct;

impl<T: Coefficient, M: Metric> Multivector<T, M> {
    /// The exponential.
    ///
    /// When the argument squares to a scalar `a²` this is a closed form: `cos a + sin a / a · x`
    /// for `a² < 0`, `cosh a + sinh a / a · x` for `a² > 0`, and `1 + x` for `a² = 0`, which covers
    /// rotors, boosts, translators and simple PGA lines. The result stores the blades of `1 + x`.
    ///
    /// Otherwise it is a 14-term Taylor series after scaling and squaring: the argument is halved
    /// until the sum of its absolute coefficients is at most 1/2, and the series result is squared
    /// as often. That sum bounds every power of the argument (each blade pair lands on one blade
    /// with sign ±1 or 0), so the first omitted term is below 0.5¹⁵/15! ≈ 2·10⁻¹⁷ relative to the
    /// result and the series is exact to `f64` precision. Each squaring can double the rounding
    /// error, so arguments with a large coefficient sum lose about one bit per doubling of that
    /// sum; this is validated against a 50-digit reference for coefficient sums up to 91.
    ///
    /// # Examples
    ///
    /// ```
    /// use seagrass_math::prelude::*;
    ///
    /// // A screw motion: rotation and translation about the same axis.
    /// let screw = Multivector::<f64, Pga3>::from_str("-0.5e12 + 2e34").unwrap();
    /// let motor = screw.exp();
    /// let point = Multivector::<f64, Pga3>::new_point(1.0, 0.0, 0.0);
    /// let moved = motor.sandwich(point.as_multivector());
    /// let p = moved.point_to_euclidean();
    /// assert!((p - Vector3::new(1f64.cos(), 1f64.sin(), 4.0)).norm() < 1e-12);
    /// ```
    #[must_use]
    pub fn exp(&self) -> Self {
        let square = self.geometric_product(self);
        // A square that stores no blades is zero: every pair met on a null vector, as for an
        // ideal line in projective geometric algebra.
        if square.is_zero() || square.is_scalar() {
            let alpha = square.scalar();
            return if alpha < T::zero() {
                let a = (-alpha).sqrt();
                let (sin_a, cos_a) = a.sin_cos();
                self * sin_a / a + cos_a
            } else if alpha == T::zero() {
                self + T::one()
            } else {
                let a = alpha.sqrt();
                self * a.sinh() / a + a.cosh()
            };
        }

        let half = constant::<T>(0.5);
        let mut l1 = self
            .coefficients
            .iter()
            .fold(T::zero(), |sum, c| sum + c.abs());
        let mut squarings = 0;
        while l1 > half && l1 < T::max_value() {
            l1 *= half;
            squarings += 1;
        }
        let scaled = self * half.powi(squarings);

        let mut power = scaled.geometric_product(&scaled);
        let mut out = &scaled + T::one();
        out += &(&power / constant::<T>(2.0));
        let mut factorial = 2.0;
        for n in 3..=14 {
            factorial *= f64::from(n);
            power = power.geometric_product(&scaled);
            out += &(&power / constant::<T>(factorial));
        }

        for _ in 0..squarings {
            out = out.geometric_product(&out);
        }
        out
    }
}
