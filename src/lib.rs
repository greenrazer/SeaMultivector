#![doc = include_str!("../readme.md")]

mod backend;
mod batch;
mod blade;
mod coefficient;
mod compat;
mod error;
mod macros;
mod metric;
mod multivector;
mod pga3d;
mod plan;
mod traits;
mod typed;

#[cfg(feature = "simd")]
pub use backend::SimdBackend;
pub use backend::{Backend, CpuBackend};
#[cfg(feature = "gpu")]
pub use backend::{GpuBackend, GpuBatch, GpuError, GpuOptions, ShaderMode};
pub use batch::{Batch, Product};
pub use blade::{BladeIndex, BladeSet, SetBlades};
pub use coefficient::Coefficient;
// Re-exported so 0.1 code keeps compiling, with a warning at the use site.
#[allow(deprecated)]
pub use compat::{
    ParseMetricError, ParseMultivectorError, SeaBasisBladeIndexLow, SeaFlectorPGA3D, SeaMotorPGA3D,
    SeaMultivector3D, SeaMultivector4D, SeaMultivector5D, SeaMultivector6D, SeaMultivectorPGA3D,
    SeaVersorPGA3D,
};
pub use error::{BladeError, KindError, MetricError, ParseError};
pub use metric::{
    Cga2, Cga3, Dyn3, Dyn4, Dyn5, Dyn6, Dynamic, Euclid3, Euclid4, Euclid5, Euclid6, Metric, Pga2,
    Pga3, Pga4, Signature, Sta,
};
pub use multivector::{Iter, IterMut, Multivector};
pub use pga3d::MotorCreationOrder;
pub use traits::{
    GeometricProduct, InnerProduct, LeftContraction, OuterProduct, RegressiveProduct,
    RightContraction, SandwichProduct, ScalarProduct,
};
pub use typed::{Any, Even, Flector, Kind, KindMul, Motor, Odd, Typed, Versor};

/// Runs the code in `ARCHITECTURE.md` as doctests, as the crate docs do for the readme.
#[cfg(doctest)]
#[doc = include_str!("../ARCHITECTURE.md")]
pub struct ArchitectureDoctests;

/// The types and traits a typical program needs, plus [`FromStr`](std::str::FromStr) and the
/// nalgebra types used by the PGA3D conversions.
pub mod prelude {
    pub use std::str::FromStr;

    #[cfg(feature = "nalgebra")]
    pub use nalgebra::{Matrix4, Vector3};

    #[cfg(feature = "gpu")]
    pub use crate::{GpuBackend, GpuOptions, ShaderMode};

    #[cfg(feature = "simd")]
    pub use crate::SimdBackend;

    pub use crate::{
        Any, Backend, Batch, BladeError, BladeIndex, Cga2, Cga3, Coefficient, CpuBackend, Dyn3,
        Dyn4, Dyn5, Dyn6, Dynamic, Euclid3, Euclid4, Euclid5, Euclid6, Even, Flector,
        GeometricProduct, InnerProduct, Kind, KindError, LeftContraction, Metric, MetricError,
        Motor, MotorCreationOrder, Multivector, Odd, OuterProduct, ParseError, Pga2, Pga3, Pga4,
        Product, RegressiveProduct, RightContraction, SandwichProduct, ScalarProduct, Signature,
        Sta, Typed, Versor,
    };
    // Re-exported so 0.1 code keeps compiling, with a warning at the use site.
    #[allow(deprecated)]
    pub use crate::{
        ParseMetricError, ParseMultivectorError, SeaBasisBladeIndexLow, SeaFlectorPGA3D,
        SeaMotorPGA3D, SeaMultivector3D, SeaMultivector4D, SeaMultivector5D, SeaMultivector6D,
        SeaMultivectorPGA3D, SeaVersorPGA3D,
    };
}
