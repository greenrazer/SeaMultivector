use super::*;

/// Reversion, grade involution, Clifford conjugation and negation: each flips the signs of the
/// stored blades in one constant mask.
impl<T: Coefficient, M: Metric> Multivector<T, M> {
    /// The reverse: each blade's basis vectors in the opposite order, which negates grades 2 and 3
    /// modulo 4.
    ///
    /// # Examples
    ///
    /// ```
    /// use seagrass_math::prelude::*;
    ///
    /// let a = Multivector::<f64, Euclid3>::from_str("1 + e1 + e12 + e123").unwrap();
    /// assert_eq!(a.reversed().to_string(), "3,0,0>1 + e1 - e12 - e123");
    /// ```
    #[must_use]
    pub fn reversed(&self) -> Self {
        self.negated_blades(M::Bits::REVERSE_SIGNS)
    }

    /// Reverses in place; see [`reversed`](Self::reversed).
    pub fn reverse(&mut self) {
        self.negate_blades(M::Bits::REVERSE_SIGNS);
    }

    /// The grade involution, which negates the odd grades.
    #[must_use]
    pub fn involuted(&self) -> Self {
        self.negated_blades(M::Bits::INVOLUTE_SIGNS)
    }

    /// Applies the grade involution in place.
    pub fn involute(&mut self) {
        self.negate_blades(M::Bits::INVOLUTE_SIGNS);
    }

    /// The Clifford conjugate, the reverse of the grade involution, which negates grades 1 and 2
    /// modulo 4.
    #[must_use]
    pub fn conjugated(&self) -> Self {
        self.negated_blades(M::Bits::CONJUGATE_SIGNS)
    }

    /// Applies the Clifford conjugate in place.
    pub fn conjugate(&mut self) {
        self.negate_blades(M::Bits::CONJUGATE_SIGNS);
    }

    /// The negation, `-self`.
    #[must_use]
    pub fn negated(&self) -> Self {
        let mut out = self.clone();
        out.negate();
        out
    }

    /// Negates in place.
    pub fn negate(&mut self) {
        for c in &mut self.coefficients {
            *c = -*c;
        }
    }
}
