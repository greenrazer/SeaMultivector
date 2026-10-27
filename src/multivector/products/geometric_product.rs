use super::*;

use crate::traits::GeometricProduct;

impl<T: Coefficient, M: Metric> GeometricProduct for Multivector<T, M> {
    type Output = Self;

    /// # Panics
    ///
    /// Panics if the operands have different runtime metrics.
    #[inline]
    fn geometric_product(&self, rhs: &Self) -> Self {
        self.planned_product(Op::Geometric, rhs)
    }
}
