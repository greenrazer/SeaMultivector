use super::*;

use crate::traits::InnerProduct;

impl<T: Coefficient, M: Metric> InnerProduct for Multivector<T, M> {
    type Output = Self;

    /// # Panics
    ///
    /// Panics if the operands have different runtime metrics.
    #[inline]
    fn inner_product(&self, rhs: &Self) -> Self {
        self.planned_product(Op::Inner, rhs)
    }
}
