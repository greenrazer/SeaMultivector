
use std::ops::{Deref, DerefMut};
use std::str::FromStr;

use crate::coefficient::SeaCoefficient;

mod basis_blade_index;
mod basis_blade;
mod basis_blade_bits;
mod without_metric;
mod with_metric;

use with_metric::SeaMultivectorLowSparseCoefficents;


pub type Multivector3D<T> = SeaMultivectorLowSparseCoefficents::<T, u8>;
pub type Multivector4D<T> = SeaMultivectorLowSparseCoefficents::<T, u16>;
pub type Multivector5D<T> = SeaMultivectorLowSparseCoefficents::<T, u32>;
pub type Multivector6D<T> = SeaMultivectorLowSparseCoefficents::<T, u64>;