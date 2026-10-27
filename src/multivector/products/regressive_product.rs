use super::*;

use crate::traits::{OuterProduct, RegressiveProduct};

impl<T: Coefficient, M: Metric> RegressiveProduct for Multivector<T, M> {
    type Output = Self;

    /// The outer product of the duals, undualized.
    ///
    /// # Panics
    ///
    /// Panics if the operands have different runtime metrics.
    fn regressive_product(&self, rhs: &Self) -> Self {
        self.assert_same_metric(rhs);
        let mut out = self.dual().outer_product(&rhs.dual());
        out.undualize();
        out
    }
}
