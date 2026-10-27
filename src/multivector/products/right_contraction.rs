use super::*;

use crate::traits::RightContraction;

impl<T: Coefficient, M: Metric> RightContraction for Multivector<T, M> {
    type Output = Self;

    /// # Panics
    ///
    /// Panics if the operands have different runtime metrics.
    #[inline]
    fn right_contraction(&self, rhs: &Self) -> Self {
        self.planned_product(Op::RightContraction, rhs)
    }
}
