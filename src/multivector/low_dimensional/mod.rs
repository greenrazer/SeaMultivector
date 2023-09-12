// use std::ops::{Deref, DerefMut};
// use std::str::FromStr;

// use crate::coefficient::SeaCoefficient;
// pub use basis_blade_set::SeaBasisBladeSet;
// pub use basis_blade_index::SeaBasisBladeIndexLow;

mod core;
mod custom_metric;
pub mod pga3d;
pub mod basis_blade_set;

use custom_metric::SeaSparseMultivectorLowDimensionsWithMetric;

pub type Multivector3D<T> = SeaSparseMultivectorLowDimensionsWithMetric::<T, u8>;
pub type Multivector4D<T> = SeaSparseMultivectorLowDimensionsWithMetric::<T, u16>;
pub type Multivector5D<T> = SeaSparseMultivectorLowDimensionsWithMetric::<T, u32>;
pub type Multivector6D<T> = SeaSparseMultivectorLowDimensionsWithMetric::<T, u64>;