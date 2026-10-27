use super::*;

use crate::plan::with_sandwich_plan;
use crate::traits::SandwichProduct;

impl<T: Coefficient, M: Metric> SandwichProduct for Multivector<T, M> {
    type Output = Self;

    /// `self * x * reverse(self)` through a cached sandwich plan,
    /// without building the reverse. Blades whose coefficients cancel are not stored: terms
    /// that cancel symbolically are dropped from the plan, and coefficients that come out
    /// exactly zero are removed.
    ///
    /// # Panics
    ///
    /// Panics if the operands have different runtime metrics.
    fn sandwich_product(&self, x: &Self) -> Self {
        self.assert_same_metric(x);
        let mut out = with_sandwich_plan(&self.metric, self.blades, x.blades, |plan| {
            Self::from_raw(
                self.metric,
                plan.output(),
                plan.run(&self.coefficients, &x.coefficients),
            )
        });
        if out.coefficients.iter().any(|c| c.is_zero()) {
            out.remove_zeros();
        }
        out
    }
}

impl<T: Coefficient, M: Metric> Multivector<T, M> {
    /// The sandwich product `self * x * reverse(self)`, which applies a versor to `x`.
    ///
    /// # Panics
    ///
    /// Panics if the operands have different runtime metrics.
    ///
    /// # Examples
    ///
    /// ```
    /// use seagrass_math::prelude::*;
    ///
    /// // A rotor by 90 degrees in the e12 plane turns e1 into e2.
    /// let (s, c) = std::f64::consts::FRAC_PI_4.sin_cos();
    /// let rotor = Multivector::<f64, Euclid3>::from_str(&format!("{c} - {s}e12")).unwrap();
    /// let e1 = Multivector::<f64, Euclid3>::basis(Euclid3::default(), BladeIndex::E1);
    /// let mut turned = rotor.sandwich(&e1);
    /// turned.simplify();
    /// assert_eq!(turned.to_string(), "3,0,0>e2");
    /// ```
    pub fn sandwich(&self, x: &Self) -> Self {
        self.sandwich_product(x)
    }
}
