use super::*;

impl<T: Coefficient, M: Metric> Multivector<T, M> {
    /// The norm `sqrt(|<self * reverse(self)>₀|)`.
    ///
    /// Blades that contain a null vector square to zero, so in degenerate metrics they do not
    /// contribute: the norm of a PGA motor is that of its rotation part.
    ///
    /// # Examples
    ///
    /// ```
    /// use seagrass_math::prelude::*;
    ///
    /// let a = Multivector::<f64, Euclid3>::from_str("3e1 + 4e12").unwrap();
    /// assert_eq!(a.norm(), 5.0);
    /// ```
    pub fn norm(&self) -> T {
        let constants = self.constants();
        let signs =
            M::Bits::product_signs(BladeIndex::SCALAR, constants.negative) ^ M::Bits::REVERSE_SIGNS;
        let squares_to_nonzero = M::Bits::sub_blades(constants.non_null);

        let mut out = T::zero();
        for (x, &c) in self.iter() {
            if squares_to_nonzero.has(x) {
                out += if signs.has(x) { -(c * c) } else { c * c };
            }
        }
        out.abs().sqrt()
    }

    /// The value divided by its [`norm`](Self::norm).
    #[must_use]
    pub fn normalized(&self) -> Self {
        self / self.norm()
    }

    /// Divides the value by its [`norm`](Self::norm) in place.
    pub fn normalize(&mut self) {
        let norm = self.norm();
        *self /= norm;
    }
}
