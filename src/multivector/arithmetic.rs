use super::*;

use std::ops::{AddAssign, DivAssign, MulAssign, Neg, SubAssign};

impl<T: Coefficient, M: Metric> Multivector<T, M> {
    /// Combines `rhs` into `self` blade by blade with `op`, storing the union of both blade sets.
    fn merge(&mut self, rhs: &Self, op: impl Fn(&mut T, T)) {
        self.assert_same_metric(rhs);
        if self.blades == rhs.blades {
            for (c, r) in self.coefficients.iter_mut().zip(&rhs.coefficients) {
                op(c, *r);
            }
            return;
        }
        let blades = self.blades | rhs.blades;
        let mut coefficients = Coefficients::from_elem(T::zero(), blades.count_ones() as usize);
        let (lhs, rhs_coefficients) = (self.coefficients.as_slice(), rhs.coefficients.as_slice());
        let (mut lhs_slot, mut rhs_slot) = (0, 0);
        for (index, c) in blades.iter().zip(coefficients.iter_mut()) {
            if self.blades.has(index) {
                *c = lhs[lhs_slot];
                lhs_slot += 1;
            }
            if rhs.blades.has(index) {
                op(c, rhs_coefficients[rhs_slot]);
                rhs_slot += 1;
            }
        }
        self.blades = blades;
        self.coefficients = coefficients;
    }

    /// Adds `value` to the coefficient of `index`, storing the blade if it was not stored.
    fn add_to_blade(&mut self, index: BladeIndex, value: T) {
        let slot = self.blades.count_below(index) as usize;
        if !self.blades.has(index) {
            self.coefficients.insert(slot, T::zero());
            self.blades |= M::Bits::from_index(index);
        }
        self.coefficients[slot] += value;
    }
}

impl<T: Coefficient, M: Metric> AddAssign<&Multivector<T, M>> for Multivector<T, M> {
    /// # Panics
    ///
    /// Panics if the operands have different runtime metrics.
    fn add_assign(&mut self, rhs: &Self) {
        self.merge(rhs, |c, r| *c += r);
    }
}

impl<T: Coefficient, M: Metric> SubAssign<&Multivector<T, M>> for Multivector<T, M> {
    /// # Panics
    ///
    /// Panics if the operands have different runtime metrics.
    fn sub_assign(&mut self, rhs: &Self) {
        self.merge(rhs, |c, r| *c -= r);
    }
}

impl<T: Coefficient, M: Metric> AddAssign<T> for Multivector<T, M> {
    /// Adds a scalar.
    fn add_assign(&mut self, rhs: T) {
        self.add_to_blade(BladeIndex::SCALAR, rhs);
    }
}

impl<T: Coefficient, M: Metric> SubAssign<T> for Multivector<T, M> {
    /// Subtracts a scalar.
    fn sub_assign(&mut self, rhs: T) {
        self.add_to_blade(BladeIndex::SCALAR, -rhs);
    }
}

impl<T: Coefficient, M: Metric> MulAssign<T> for Multivector<T, M> {
    /// Scales every coefficient.
    fn mul_assign(&mut self, rhs: T) {
        for c in &mut self.coefficients {
            *c *= rhs;
        }
    }
}

impl<T: Coefficient, M: Metric> DivAssign<T> for Multivector<T, M> {
    /// Divides every coefficient.
    fn div_assign(&mut self, rhs: T) {
        for c in &mut self.coefficients {
            *c /= rhs;
        }
    }
}

impl<T: Coefficient, M: Metric> Neg for &Multivector<T, M> {
    type Output = Multivector<T, M>;

    fn neg(self) -> Multivector<T, M> {
        self.negated()
    }
}

impl<T: Coefficient, M: Metric> Neg for Multivector<T, M> {
    type Output = Self;

    fn neg(mut self) -> Self {
        self.negate();
        self
    }
}
