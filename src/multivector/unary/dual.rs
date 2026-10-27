use super::*;

use crate::metric::pseudoscalar_inverse_negative;

/// Duals and complements are a bit reversal plus a sign mask.
///
/// Blade i's complement has index `2^n - 1 - i`, every bit flipped, and since the bitset holds
/// exactly 2^n bits, `reverse_bits` sends position i there. Coefficients are packed in index
/// order, so reversing them puts each in its complement's slot. No product with the
/// pseudoscalar is involved, which is why the dual also works when the pseudoscalar has no
/// inverse, as in projective geometric algebra.
impl<T: Coefficient, M: Metric> Multivector<T, M> {
    /// The metric-free dual: every blade moved to its complement, with no sign changes.
    #[must_use]
    pub fn poincare_dual(&self) -> Self {
        let mut out = self.clone();
        out.poincare_dualize();
        out
    }

    /// Replaces the value with its [`poincare_dual`](Self::poincare_dual).
    pub fn poincare_dualize(&mut self) {
        self.blades = self.blades.reverse_bits();
        self.coefficients.reverse();
    }

    /// The right complement r of each blade x, defined by `r ^ x = I`.
    ///
    /// # Examples
    ///
    /// ```
    /// use seagrass_math::prelude::*;
    ///
    /// let e1 = Multivector::<f64, Pga3>::basis(Pga3::default(), BladeIndex::E1);
    /// let r = e1.right_complemented();
    /// assert_eq!(r.to_string(), "3,0,1>-e234");
    /// assert_eq!((&r ^ &e1).to_string(), "3,0,1>e1234");
    /// ```
    #[must_use]
    pub fn right_complemented(&self) -> Self {
        let mut out = self.poincare_dual();
        out.negate_blades(M::Bits::RIGHT_COMPLEMENT_SIGNS);
        out
    }

    /// Replaces the value with its [`right_complemented`](Self::right_complemented) form.
    pub fn right_complement(&mut self) {
        self.poincare_dualize();
        self.negate_blades(M::Bits::RIGHT_COMPLEMENT_SIGNS);
    }

    /// The left complement l of each blade x, defined by `x ^ l = I`.
    #[must_use]
    pub fn left_complemented(&self) -> Self {
        let mut out = self.poincare_dual();
        out.negate_blades(M::Bits::LEFT_COMPLEMENT_SIGNS);
        out
    }

    /// Replaces the value with its [`left_complemented`](Self::left_complemented) form.
    pub fn left_complement(&mut self) {
        self.poincare_dualize();
        self.negate_blades(M::Bits::LEFT_COMPLEMENT_SIGNS);
    }

    /// The dual. Without null vectors it is `self * I⁻¹`; with null vectors, where I has no
    /// inverse, it is the [right complement](Self::right_complemented).
    #[must_use]
    pub fn dual(&self) -> Self {
        let mut out = self.poincare_dual();
        out.apply_dual_signs(false);
        out
    }

    /// Replaces the value with its [`dual`](Self::dual).
    pub fn dualize(&mut self) {
        self.poincare_dualize();
        self.apply_dual_signs(false);
    }

    /// The inverse of [`dual`](Self::dual): `self * I` without null vectors, the left complement
    /// with them.
    #[must_use]
    pub fn undual(&self) -> Self {
        let mut out = self.poincare_dual();
        out.apply_dual_signs(true);
        out
    }

    /// Replaces the value with its [`undual`](Self::undual).
    pub fn undualize(&mut self) {
        self.poincare_dualize();
        self.apply_dual_signs(true);
    }

    fn apply_dual_signs(&mut self, inverse: bool) {
        if self.metric.zero() > 0 {
            let signs = if inverse {
                M::Bits::LEFT_COMPLEMENT_SIGNS
            } else {
                M::Bits::RIGHT_COMPLEMENT_SIGNS
            };
            self.negate_blades(signs);
            return;
        }
        let signs = M::Bits::dual_signs(
            inverse,
            pseudoscalar_inverse_negative(&self.metric),
            self.constants().negative,
        );
        self.negate_blades(signs);
    }
}
