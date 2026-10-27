use super::*;

use crate::traits::ScalarProduct;

impl<T: Coefficient, M: Metric> ScalarProduct for Multivector<T, M> {
    type Output = T;

    /// # Panics
    ///
    /// Panics if the operands have different runtime metrics.
    fn scalar_product(&self, rhs: &Self) -> T {
        self.assert_same_metric(rhs);
        let constants = self.constants();
        // Only a blade times itself lands on the scalar, and only when it holds no null vector.
        let pairs = self.blades & rhs.blades & M::Bits::sub_blades(constants.non_null);
        let signs = M::Bits::product_signs(BladeIndex::SCALAR, constants.negative);

        let mut out = T::zero();
        for x in pairs.iter() {
            let c = self.coefficients[self.blades.count_below(x) as usize]
                * rhs.coefficients[rhs.blades.count_below(x) as usize];
            out += if signs.has(x) { -c } else { c };
        }
        out
    }
}
