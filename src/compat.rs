//! Names from 0.1, kept for one release as deprecated aliases.

use crate::blade::BladeIndex;
use crate::coefficient::Coefficient;
use crate::error::{MetricError, ParseError};
use crate::metric::{Dyn3, Dyn4, Dyn5, Dyn6, Metric, Pga3};
use crate::multivector::Multivector;
use crate::typed::{Flector, Kind, Motor, Typed, Versor};

/// A 3D multivector with a runtime metric.
#[deprecated(since = "0.2.0", note = "use `Multivector<T, Dyn3>`")]
pub type SeaMultivector3D<T> = Multivector<T, Dyn3>;
/// A 4D multivector with a runtime metric.
#[deprecated(since = "0.2.0", note = "use `Multivector<T, Dyn4>`")]
pub type SeaMultivector4D<T> = Multivector<T, Dyn4>;
/// A 5D multivector with a runtime metric.
#[deprecated(since = "0.2.0", note = "use `Multivector<T, Dyn5>`")]
pub type SeaMultivector5D<T> = Multivector<T, Dyn5>;
/// A 6D multivector with a runtime metric.
#[deprecated(since = "0.2.0", note = "use `Multivector<T, Dyn6>`")]
pub type SeaMultivector6D<T> = Multivector<T, Dyn6>;
/// A PGA3D multivector.
#[deprecated(since = "0.2.0", note = "use `Multivector<T, Pga3>`")]
pub type SeaMultivectorPGA3D<T> = Multivector<T, Pga3>;
/// A PGA3D motor.
#[deprecated(since = "0.2.0", note = "use `Motor<T, Pga3>`")]
pub type SeaMotorPGA3D<T> = Motor<T, Pga3>;
/// A PGA3D flector.
#[deprecated(since = "0.2.0", note = "use `Flector<T, Pga3>`")]
pub type SeaFlectorPGA3D<T> = Flector<T, Pga3>;
/// A PGA3D versor.
#[deprecated(since = "0.2.0", note = "use `Versor<T, Pga3>`")]
pub type SeaVersorPGA3D<T> = Versor<T, Pga3>;
/// A basis blade index.
#[deprecated(since = "0.2.0", note = "use `BladeIndex`")]
pub type SeaBasisBladeIndexLow = BladeIndex;
/// A parse error.
#[deprecated(since = "0.2.0", note = "use `ParseError`")]
pub type ParseMultivectorError = ParseError;
/// A metric parse error.
#[deprecated(since = "0.2.0", note = "use `MetricError`")]
pub type ParseMetricError = MetricError;

/// Method names from 0.1 that only changed spelling or form.
macro_rules! renamed {
    (&self: $($old:ident => $new:ident -> $Out:ty, $note:literal;)*) => {$(
        #[deprecated(since = "0.2.0", note = $note)]
        #[doc = concat!("Renamed to `", stringify!($new), "`.")]
        pub fn $old(&self) -> $Out {
            self.$new()
        }
    )*};
    (&mut self: $($old:ident => $new:ident, $note:literal;)*) => {$(
        #[deprecated(since = "0.2.0", note = $note)]
        #[doc = concat!("Renamed to `", stringify!($new), "`.")]
        pub fn $old(&mut self) {
            self.$new()
        }
    )*};
}

impl<T: Coefficient, M: Metric> Multivector<T, M> {
    renamed! {
        &self:
        reversion_conjugated => reversed -> Self, "use `reversed`";
        involution_conjugated => involuted -> Self, "use `involuted`";
        clifford_conjugated => conjugated -> Self, "use `conjugated`";
        left_complimented => left_complemented -> Self, "use `left_complemented`";
        right_complimented => right_complemented -> Self, "use `right_complemented`";
        scalar_normalized => normalized -> Self, "use `normalized`";
        scalar_norm => norm -> T, "use `norm`";
    }
    renamed! {
        &mut self:
        reversion_conjugation => reverse, "use `reverse`";
        involution_conjugation => involute, "use `involute`";
        clifford_conjugation => conjugate, "use `conjugate`";
        dualization => dualize, "use `dualize`";
        undualization => undualize, "use `undualize`";
        poincare_dualization => poincare_dualize, "use `poincare_dualize`";
        left_compliment => left_complement, "use `left_complement`";
        right_compliment => right_complement, "use `right_complement`";
        scalar_normalize => normalize, "use `normalize`";
    }

    /// Renamed to `get`.
    #[deprecated(since = "0.2.0", note = "use `get`")]
    pub fn basis_blade_coefficient(&self, index: &BladeIndex) -> T {
        self.get(*index)
    }
}

impl<T: Coefficient, M: Metric, K: Kind> Typed<T, M, K> {
    renamed! {
        &self:
        reversion_conjugated => reversed -> Self, "use `reversed`";
        involution_conjugated => involuted -> Self, "use `involuted`";
        clifford_conjugated => conjugated -> Self, "use `conjugated`";
        scalar_normalized => normalized -> Self, "use `normalized`";
        scalar_norm => norm -> T, "use `norm`";
    }
    renamed! {
        &mut self:
        reversion_conjugation => reverse, "use `reverse`";
        involution_conjugation => involute, "use `involute`";
        clifford_conjugation => conjugate, "use `conjugate`";
        scalar_normalize => normalize, "use `normalize`";
    }

    /// Renamed to `get`.
    #[deprecated(since = "0.2.0", note = "use `get`")]
    pub fn basis_blade_coefficient(&self, index: &BladeIndex) -> T {
        self.get(*index)
    }
}
