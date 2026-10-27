use super::*;

use std::iter::Zip;
use std::slice;

use crate::blade::SetBlades;

/// Iterator over the stored blades of a multivector and their coefficients.
pub type Iter<'a, T, B> = Zip<SetBlades<B>, slice::Iter<'a, T>>;

/// Iterator over the stored blades of a multivector and their coefficients, mutably.
pub type IterMut<'a, T, B> = Zip<SetBlades<B>, slice::IterMut<'a, T>>;

impl<T: Coefficient, M: Metric> Multivector<T, M> {
    /// The stored blades and their coefficients, in blade-index order.
    ///
    /// # Examples
    ///
    /// ```
    /// use seagrass_math::prelude::*;
    ///
    /// let a = Multivector::<f64, Euclid3>::from_str("2 + 3e12").unwrap();
    /// let blades: Vec<_> = a.iter().map(|(index, c)| (index, *c)).collect();
    /// assert_eq!(blades, [(BladeIndex::SCALAR, 2.0), (BladeIndex::E12, 3.0)]);
    /// ```
    pub fn iter(&self) -> Iter<'_, T, M::Bits> {
        self.blades.iter().zip(self.coefficients.iter())
    }

    /// The stored blades and their coefficients, mutably.
    pub fn iter_mut(&mut self) -> IterMut<'_, T, M::Bits> {
        self.blades.iter().zip(self.coefficients.iter_mut())
    }

    /// The stored blades of grade `k` and their coefficients.
    ///
    /// # Panics
    ///
    /// Panics if `k` is above 6.
    pub fn grade_iter(&self, k: u32) -> impl Iterator<Item = (BladeIndex, &T)> {
        let mask = M::Bits::GRADE[k as usize];
        self.iter().filter(move |(index, _)| mask.has(*index))
    }

    /// The stored blades of grade `k` and their coefficients, mutably.
    ///
    /// # Panics
    ///
    /// Panics if `k` is above 6.
    pub fn grade_iter_mut(&mut self, k: u32) -> impl Iterator<Item = (BladeIndex, &mut T)> {
        let mask = M::Bits::GRADE[k as usize];
        self.iter_mut().filter(move |(index, _)| mask.has(*index))
    }
}

impl<'a, T: Coefficient, M: Metric> IntoIterator for &'a Multivector<T, M> {
    type Item = (BladeIndex, &'a T);
    type IntoIter = Iter<'a, T, M::Bits>;

    fn into_iter(self) -> Self::IntoIter {
        self.iter()
    }
}

impl<'a, T: Coefficient, M: Metric> IntoIterator for &'a mut Multivector<T, M> {
    type Item = (BladeIndex, &'a mut T);
    type IntoIter = IterMut<'a, T, M::Bits>;

    fn into_iter(self) -> Self::IntoIter {
        self.iter_mut()
    }
}
