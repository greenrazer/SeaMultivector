use std::ops::{Add, Mul, Neg, AddAssign, Deref, DerefMut, BitXorAssign, BitAndAssign, BitOrAssign, Index, BitXor, BitAnd, BitOr, IndexMut, Sub, Shl, Not, Shr, ShrAssign, ShlAssign, MulAssign, Div};
use std::fmt::{Display, Debug};
use std::str::FromStr;

use num_traits::real::Real;
use num_traits::{Zero, One};

pub trait SeaCoefficient: 
	Clone + Copy + Display + Debug + FromStr
	+ One + Zero + Real
	+ AddAssign
	+ MulAssign
{}
impl<T> SeaCoefficient for T 
	where 
		T: 
			Clone + Copy + Display + Debug + FromStr
			+ One + Zero + Real 
			+ AddAssign
			+ MulAssign
{}