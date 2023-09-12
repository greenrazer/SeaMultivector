// use std::ops::{Deref, DerefMut};
// use std::str::FromStr;

// use crate::coefficient::SeaCoefficient;
// pub use basis_blade_set::SeaBasisBladeSet;
// pub use basis_blade_index::SeaBasisBladeIndexLow;

mod core;
mod custom_metric;
pub mod pga3d;
pub mod basis_blade_set;

use custom_metric::SeaMultivectorLowSparseCoefficentsWithMetric;

pub type Multivector3D<T> = SeaMultivectorLowSparseCoefficentsWithMetric::<T, u8>;
pub type Multivector4D<T> = SeaMultivectorLowSparseCoefficentsWithMetric::<T, u16>;
pub type Multivector5D<T> = SeaMultivectorLowSparseCoefficentsWithMetric::<T, u32>;
pub type Multivector6D<T> = SeaMultivectorLowSparseCoefficentsWithMetric::<T, u64>;