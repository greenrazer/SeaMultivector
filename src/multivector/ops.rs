use super::*;

use std::ops::{
    Add, AddAssign, BitAnd, BitOr, BitXor, Div, DivAssign, Mul, MulAssign, Not, Shl, Shr, Sub,
    SubAssign,
};

use crate::macros::{binary_operator, forward_assign};
use crate::traits::{
    GeometricProduct, InnerProduct, LeftContraction, OuterProduct, RegressiveProduct,
    RightContraction,
};

binary_operator!([T: Coefficient, M: Metric] Mul::mul(a: Multivector<T, M>, b: Multivector<T, M>) -> Multivector<T, M> {
    a.geometric_product(b)
});
binary_operator!([T: Coefficient, M: Metric] BitXor::bitxor(a: Multivector<T, M>, b: Multivector<T, M>) -> Multivector<T, M> {
    a.outer_product(b)
});
binary_operator!([T: Coefficient, M: Metric] BitOr::bitor(a: Multivector<T, M>, b: Multivector<T, M>) -> Multivector<T, M> {
    a.inner_product(b)
});
binary_operator!([T: Coefficient, M: Metric] Shl::shl(a: Multivector<T, M>, b: Multivector<T, M>) -> Multivector<T, M> {
    a.left_contraction(b)
});
binary_operator!([T: Coefficient, M: Metric] Shr::shr(a: Multivector<T, M>, b: Multivector<T, M>) -> Multivector<T, M> {
    a.right_contraction(b)
});
binary_operator!([T: Coefficient, M: Metric] BitAnd::bitand(a: Multivector<T, M>, b: Multivector<T, M>) -> Multivector<T, M> {
    a.regressive_product(b)
});

forward_assign!([T: Coefficient, M: Metric] AddAssign::add_assign(Multivector<T, M>, Multivector<T, M>));
forward_assign!([T: Coefficient, M: Metric] SubAssign::sub_assign(Multivector<T, M>, Multivector<T, M>));

/// `+` and `-` between multivectors reuse an owned left operand's storage.
macro_rules! sum_operator {
    ($Op:ident::$op:ident, $OpAssign:ident::$op_assign:ident) => {
        impl<T: Coefficient, M: Metric> $Op<&Multivector<T, M>> for &Multivector<T, M> {
            type Output = Multivector<T, M>;

            #[inline]
            fn $op(self, rhs: &Multivector<T, M>) -> Multivector<T, M> {
                let mut out = self.clone();
                $OpAssign::$op_assign(&mut out, rhs);
                out
            }
        }

        impl<T: Coefficient, M: Metric> $Op<Multivector<T, M>> for &Multivector<T, M> {
            type Output = Multivector<T, M>;

            #[inline]
            fn $op(self, rhs: Multivector<T, M>) -> Multivector<T, M> {
                $Op::$op(self, &rhs)
            }
        }

        impl<T: Coefficient, M: Metric> $Op<&Multivector<T, M>> for Multivector<T, M> {
            type Output = Multivector<T, M>;

            #[inline]
            fn $op(mut self, rhs: &Multivector<T, M>) -> Multivector<T, M> {
                $OpAssign::$op_assign(&mut self, rhs);
                self
            }
        }

        impl<T: Coefficient, M: Metric> $Op<Multivector<T, M>> for Multivector<T, M> {
            type Output = Multivector<T, M>;

            #[inline]
            fn $op(mut self, rhs: Multivector<T, M>) -> Multivector<T, M> {
                $OpAssign::$op_assign(&mut self, &rhs);
                self
            }
        }
    };
}

sum_operator!(Add::add, AddAssign::add_assign);
sum_operator!(Sub::sub, SubAssign::sub_assign);

/// `x + s`, `x - s`, `x * s` and `x / s` with a scalar `s`, through the `*Assign` form.
macro_rules! scalar_operator {
    ($Op:ident::$op:ident, $OpAssign:ident::$op_assign:ident) => {
        impl<T: Coefficient, M: Metric> $Op<T> for &Multivector<T, M> {
            type Output = Multivector<T, M>;

            #[inline]
            fn $op(self, rhs: T) -> Multivector<T, M> {
                let mut out = self.clone();
                $OpAssign::$op_assign(&mut out, rhs);
                out
            }
        }

        impl<T: Coefficient, M: Metric> $Op<T> for Multivector<T, M> {
            type Output = Multivector<T, M>;

            #[inline]
            fn $op(mut self, rhs: T) -> Multivector<T, M> {
                $OpAssign::$op_assign(&mut self, rhs);
                self
            }
        }
    };
}

scalar_operator!(Add::add, AddAssign::add_assign);
scalar_operator!(Sub::sub, SubAssign::sub_assign);
scalar_operator!(Mul::mul, MulAssign::mul_assign);
scalar_operator!(Div::div, DivAssign::div_assign);

impl<T: Coefficient, M: Metric> Not for &Multivector<T, M> {
    type Output = Multivector<T, M>;

    /// The [dual](Multivector::dual).
    fn not(self) -> Multivector<T, M> {
        self.dual()
    }
}

impl<T: Coefficient, M: Metric> Not for Multivector<T, M> {
    type Output = Self;

    /// The [dual](Multivector::dual).
    fn not(mut self) -> Self {
        self.dualize();
        self
    }
}
