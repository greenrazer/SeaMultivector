use std::ops::{Add, Mul, Neg, AddAssign, Deref, DerefMut, BitXorAssign, BitAndAssign, BitOrAssign, Index, BitXor, BitAnd, BitOr, IndexMut, Sub, Shl, Not, Shr, ShrAssign, ShlAssign, MulAssign, Div, DivAssign, SubAssign};
use std::fmt::{Display, Debug};
use std::str::FromStr;

use num_traits::real::Real;
use num_traits::{Zero, One};

pub trait SeaCoefficient: 
	Clone + Copy + Display + Debug + FromStr
	+ One + Zero + Real
	+ AddAssign + SubAssign + MulAssign + DivAssign
{
	fn two() -> Self {
		Self::from(2).unwrap()
	}

	fn three() -> Self {
		Self::from(3).unwrap()
	}

	fn four() -> Self {
		Self::from(4).unwrap()
	}

	fn factorial(a: u64) -> Self {
		let o: u64 = (1..=a).product();
		Self::from(o).unwrap()
	}
}
impl<T> SeaCoefficient for T 
	where 
		T: 
			Clone + Copy + Display + Debug + FromStr
			+ One + Zero + Real
			+ AddAssign + SubAssign + MulAssign + DivAssign
{}