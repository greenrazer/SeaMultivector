
use std::ops::{Deref, DerefMut};
use std::str::FromStr;

use crate::coefficient::SeaCoefficient;

mod basis_blade_index;
mod basis_blade;
mod basis_blade_bits;
mod no_metric;
mod with_metric;

use basis_blade_bits::SeaBasisBladeBitsWrapper;
use with_metric::SeaMultivectorLowSparseCoefficents;

pub type Multivector3D<T> = SeaMultivectorLowSparseCoefficents::<T, SeaBasisBladeBitsWrapper<u8>>;
pub type Multivector4D<T> = SeaMultivectorLowSparseCoefficents::<T, SeaBasisBladeBitsWrapper<u16>>;
pub type Multivector5D<T> = SeaMultivectorLowSparseCoefficents::<T, SeaBasisBladeBitsWrapper<u32>>;
pub type Multivector6D<T> = SeaMultivectorLowSparseCoefficents::<T, SeaBasisBladeBitsWrapper<u64>>;