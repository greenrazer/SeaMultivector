use super::*;

use crate::traits::LeftContraction;

impl<T: Coefficient, M: Metric> LeftContraction for Multivector<T, M> {
    type Output = Self;

    /// # Panics
    ///
    /// Panics if the operands have different runtime metrics.
    #[inline]
    fn left_contraction(&self, rhs: &Self) -> Self {
        self.planned_product(Op::LeftContraction, rhs)
    }
}
