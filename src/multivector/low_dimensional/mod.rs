// use std::ops::{Deref, DerefMut};
// use std::str::FromStr;

// use crate::coefficient::SeaCoefficient;
// pub use basis_blade_set::SeaBasisBladeSet;
// pub use basis_blade_index::SeaBasisBladeIndexLow;

mod core;
pub mod custom_metric;
pub mod pga3d;
pub mod basis_blade_set;

use custom_metric::SeaLowDimensionalMultivectorCustomMetric;

pub type SeaMultivector3D<T> = SeaLowDimensionalMultivectorCustomMetric::<T, u8>;
pub type SeaMultivector4D<T> = SeaLowDimensionalMultivectorCustomMetric::<T, u16>;
pub type SeaMultivector5D<T> = SeaLowDimensionalMultivectorCustomMetric::<T, u32>;
pub type SeaMultivector6D<T> = SeaLowDimensionalMultivectorCustomMetric::<T, u64>;