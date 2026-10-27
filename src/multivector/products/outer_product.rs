use super::*;

use crate::traits::OuterProduct;

impl<T: Coefficient, M: Metric> OuterProduct for Multivector<T, M> {
    type Output = Self;

    /// # Panics
    ///
    /// Panics if the operands have different runtime metrics.
    #[inline]
    fn outer_product(&self, rhs: &Self) -> Self {
        self.planned_product(Op::Outer, rhs)
    }
}
