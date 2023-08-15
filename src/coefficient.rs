use std::ops::{Add, Mul, Neg, AddAssign, Deref, DerefMut, BitXorAssign, BitAndAssign, BitOrAssign, Index, BitXor, BitAnd, BitOr, IndexMut, Sub, Shl, Not, Shr, ShrAssign, ShlAssign};
use std::fmt::{Display, Debug};
use std::str::FromStr;

use num_traits::{Zero, One};

pub trait SeaCoefficient: 
	Clone + Copy + Display + Debug + FromStr
	+ One + Zero 
	+ Neg<Output = Self> + Mul<Output = Self>  + Add<Output = Self> 
	+ AddAssign  
{}
impl<T> SeaCoefficient for T 
	where 
		T: 
			Clone + Copy + Display + Debug + FromStr
			+ One + Zero 
			+ Neg<Output = Self> + Mul<Output = Self>  + Add<Output = Self> 
			+ AddAssign 
{}