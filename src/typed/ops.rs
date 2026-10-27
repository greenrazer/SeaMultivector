//! Products and operators between typed values and plain multivectors.

use std::ops::{
    Add, AddAssign, BitAnd, BitOr, BitXor, Div, DivAssign, Mul, MulAssign, Neg, Not, Shl, Shr, Sub,
    SubAssign,
};

use super::*;
use crate::macros::{binary_operator, forward_assign};
use crate::traits::{
    GeometricProduct, InnerProduct, LeftContraction, OuterProduct, RegressiveProduct,
    RightContraction,
};

/// A product whose kind follows [`KindMul`], plus its mixed forms with plain multivectors.
macro_rules! typed_product {
    ($Trait:ident::$method:ident, $Op:ident::$op:ident) => {
        impl<T: Coefficient, M: Metric, K1: KindMul<K2>, K2: Kind> $Trait<Typed<T, M, K2>>
            for Typed<T, M, K1>
        {
            type Output = Typed<T, M, Product<K1, K2>>;

            fn $method(&self, rhs: &Typed<T, M, K2>) -> Self::Output {
                Typed::wrap(self.inner.$method(&rhs.inner))
            }
        }

        impl<T: Coefficient, M: Metric, K: Kind> $Trait<Multivector<T, M>> for Typed<T, M, K> {
            type Output = Multivector<T, M>;

            fn $method(&self, rhs: &Multivector<T, M>) -> Multivector<T, M> {
                self.inner.$method(rhs)
            }
        }

        impl<T: Coefficient, M: Metric, K: Kind> $Trait<Typed<T, M, K>> for Multivector<T, M> {
            type Output = Multivector<T, M>;

            fn $method(&self, rhs: &Typed<T, M, K>) -> Multivector<T, M> {
                self.$method(&rhs.inner)
            }
        }

        binary_operator!([T: Coefficient, M: Metric, K1: KindMul<K2>, K2: Kind] $Op::$op(
            a: Typed<T, M, K1>, b: Typed<T, M, K2>
        ) -> Typed<T, M, Product<K1, K2>> {
            a.$method(b)
        });
        binary_operator!([T: Coefficient, M: Metric, K: Kind] $Op::$op(
            a: Typed<T, M, K>, b: Multivector<T, M>
        ) -> Multivector<T, M> {
            a.$method(b)
        });
        binary_operator!([T: Coefficient, M: Metric, K: Kind] $Op::$op(
            a: Multivector<T, M>, b: Typed<T, M, K>
        ) -> Multivector<T, M> {
            a.$method(b)
        });
    };
}

typed_product!(GeometricProduct::geometric_product, Mul::mul);
typed_product!(OuterProduct::outer_product, BitXor::bitxor);
typed_product!(InnerProduct::inner_product, BitOr::bitor);
typed_product!(LeftContraction::left_contraction, Shl::shl);
typed_product!(RightContraction::right_contraction, Shr::shr);

/// The kind of a regressive product: grade r + s - n has the parity of r + s + n.
type Regressive<K1, K2, M> = Product<Product<K1, K2>, DimensionParity<M>>;

impl<T, M, K1, K2> RegressiveProduct<Typed<T, M, K2>> for Typed<T, M, K1>
where
    T: Coefficient,
    M: Metric,
    K1: KindMul<K2>,
    K2: Kind,
    Product<K1, K2>: KindMul<DimensionParity<M>>,
{
    type Output = Typed<T, M, Regressive<K1, K2, M>>;

    fn regressive_product(&self, rhs: &Typed<T, M, K2>) -> Self::Output {
        Typed::wrap(self.inner.regressive_product(&rhs.inner))
    }
}

impl<T: Coefficient, M: Metric, K: Kind> RegressiveProduct<Multivector<T, M>> for Typed<T, M, K> {
    type Output = Multivector<T, M>;

    fn regressive_product(&self, rhs: &Multivector<T, M>) -> Multivector<T, M> {
        self.inner.regressive_product(rhs)
    }
}

impl<T: Coefficient, M: Metric, K: Kind> RegressiveProduct<Typed<T, M, K>> for Multivector<T, M> {
    type Output = Multivector<T, M>;

    fn regressive_product(&self, rhs: &Typed<T, M, K>) -> Multivector<T, M> {
        self.regressive_product(&rhs.inner)
    }
}

/// `&` between typed values, through [`RegressiveProduct`]; written out because of the bound on
/// the doubly projected kind.
macro_rules! regressive_operator {
    ($(($Lhs:ty, $Rhs:ty, $lhs_ref:tt, $rhs_ref:tt)),*) => {$(
        impl<T, M, K1, K2> BitAnd<$Rhs> for $Lhs
        where
            T: Coefficient,
            M: Metric,
            K1: KindMul<K2>,
            K2: Kind,
            Product<K1, K2>: KindMul<DimensionParity<M>>,
        {
            type Output = Typed<T, M, Regressive<K1, K2, M>>;

            #[inline]
            fn bitand(self, rhs: $Rhs) -> Self::Output {
                regressive_operator!(@call self rhs $lhs_ref $rhs_ref)
            }
        }
    )*};
    (@call $lhs:ident $rhs:ident owned owned) => { $lhs.regressive_product(&$rhs) };
    (@call $lhs:ident $rhs:ident owned borrowed) => { $lhs.regressive_product($rhs) };
    (@call $lhs:ident $rhs:ident borrowed owned) => { $lhs.regressive_product(&$rhs) };
    (@call $lhs:ident $rhs:ident borrowed borrowed) => { $lhs.regressive_product($rhs) };
}

regressive_operator!(
    (Typed<T, M, K1>, Typed<T, M, K2>, owned, owned),
    (Typed<T, M, K1>, &Typed<T, M, K2>, owned, borrowed),
    (&Typed<T, M, K1>, Typed<T, M, K2>, borrowed, owned),
    (&Typed<T, M, K1>, &Typed<T, M, K2>, borrowed, borrowed)
);

binary_operator!([T: Coefficient, M: Metric, K: Kind] BitAnd::bitand(
    a: Typed<T, M, K>, b: Multivector<T, M>
) -> Multivector<T, M> {
    a.regressive_product(b)
});
binary_operator!([T: Coefficient, M: Metric, K: Kind] BitAnd::bitand(
    a: Multivector<T, M>, b: Typed<T, M, K>
) -> Multivector<T, M> {
    a.regressive_product(b)
});

impl<T: Coefficient, M: Metric, K: Kind> ScalarProduct<Typed<T, M, K>> for Multivector<T, M> {
    type Output = T;

    fn scalar_product(&self, rhs: &Typed<T, M, K>) -> T {
        self.scalar_product(&rhs.inner)
    }
}

impl<T: Coefficient, M: Metric, K1: Kind, K2: Kind> ScalarProduct<Typed<T, M, K2>>
    for Typed<T, M, K1>
{
    type Output = T;

    fn scalar_product(&self, rhs: &Typed<T, M, K2>) -> T {
        self.inner.scalar_product(&rhs.inner)
    }
}

impl<T: Coefficient, M: Metric, K: Kind> ScalarProduct<Multivector<T, M>> for Typed<T, M, K> {
    type Output = T;

    fn scalar_product(&self, rhs: &Multivector<T, M>) -> T {
        self.inner.scalar_product(rhs)
    }
}

impl<T: Coefficient, M: Metric, K1: Kind, K2: Kind> SandwichProduct<Typed<T, M, K2>>
    for Typed<T, M, K1>
{
    /// A versor maps a value of either parity to the same parity.
    type Output = Typed<T, M, K2>;

    fn sandwich_product(&self, x: &Typed<T, M, K2>) -> Typed<T, M, K2> {
        Typed::wrap(self.inner.sandwich_product(&x.inner))
    }
}

impl<T: Coefficient, M: Metric, K: Kind> SandwichProduct<Multivector<T, M>> for Typed<T, M, K> {
    type Output = Multivector<T, M>;

    fn sandwich_product(&self, x: &Multivector<T, M>) -> Multivector<T, M> {
        self.inner.sandwich_product(x)
    }
}

binary_operator!([T: Coefficient, M: Metric, K1: Kind, K2: Kind] Div::div(
    m: Typed<T, M, K1>, x: Typed<T, M, K2>
) -> Typed<T, M, K2> {
    m.sandwich_product(x)
});
binary_operator!([T: Coefficient, M: Metric, K: Kind] Div::div(
    m: Typed<T, M, K>, x: Multivector<T, M>
) -> Multivector<T, M> {
    m.sandwich_product(x)
});

// Sums: a sum of two values of the same kind keeps the kind; anything else is a multivector.

impl<T: Coefficient, M: Metric, K: Kind> AddAssign<&Typed<T, M, K>> for Typed<T, M, K> {
    fn add_assign(&mut self, rhs: &Self) {
        self.inner += &rhs.inner;
    }
}

impl<T: Coefficient, M: Metric, K: Kind> SubAssign<&Typed<T, M, K>> for Typed<T, M, K> {
    fn sub_assign(&mut self, rhs: &Self) {
        self.inner -= &rhs.inner;
    }
}

impl<T: Coefficient, M: Metric, K: Kind> AddAssign<&Typed<T, M, K>> for Multivector<T, M> {
    fn add_assign(&mut self, rhs: &Typed<T, M, K>) {
        *self += &rhs.inner;
    }
}

impl<T: Coefficient, M: Metric, K: Kind> SubAssign<&Typed<T, M, K>> for Multivector<T, M> {
    fn sub_assign(&mut self, rhs: &Typed<T, M, K>) {
        *self -= &rhs.inner;
    }
}

forward_assign!([T: Coefficient, M: Metric, K: Kind] AddAssign::add_assign(Typed<T, M, K>, Typed<T, M, K>));
forward_assign!([T: Coefficient, M: Metric, K: Kind] SubAssign::sub_assign(Typed<T, M, K>, Typed<T, M, K>));
forward_assign!([T: Coefficient, M: Metric, K: Kind] AddAssign::add_assign(Multivector<T, M>, Typed<T, M, K>));
forward_assign!([T: Coefficient, M: Metric, K: Kind] SubAssign::sub_assign(Multivector<T, M>, Typed<T, M, K>));

binary_operator!([T: Coefficient, M: Metric, K: Kind] Add::add(a: Typed<T, M, K>, b: Typed<T, M, K>) -> Typed<T, M, K> {
    Typed::wrap(&a.inner + &b.inner)
});
binary_operator!([T: Coefficient, M: Metric, K: Kind] Sub::sub(a: Typed<T, M, K>, b: Typed<T, M, K>) -> Typed<T, M, K> {
    Typed::wrap(&a.inner - &b.inner)
});
binary_operator!([T: Coefficient, M: Metric, K: Kind] Add::add(a: Typed<T, M, K>, b: Multivector<T, M>) -> Multivector<T, M> {
    &a.inner + b
});
binary_operator!([T: Coefficient, M: Metric, K: Kind] Sub::sub(a: Typed<T, M, K>, b: Multivector<T, M>) -> Multivector<T, M> {
    &a.inner - b
});
binary_operator!([T: Coefficient, M: Metric, K: Kind] Add::add(a: Multivector<T, M>, b: Typed<T, M, K>) -> Multivector<T, M> {
    a + &b.inner
});
binary_operator!([T: Coefficient, M: Metric, K: Kind] Sub::sub(a: Multivector<T, M>, b: Typed<T, M, K>) -> Multivector<T, M> {
    a - &b.inner
});

// Scalars: scaling keeps every kind; adding a scalar keeps only the even kind.

macro_rules! scalar_operator {
    ([$($generics:tt)*] $Lhs:ty, $Op:ident::$op:ident, $OpAssign:ident::$op_assign:ident) => {
        impl<$($generics)*> $OpAssign<T> for $Lhs {
            #[inline]
            fn $op_assign(&mut self, rhs: T) {
                $OpAssign::$op_assign(&mut self.inner, rhs);
            }
        }

        impl<$($generics)*> $Op<T> for &$Lhs {
            type Output = $Lhs;

            #[inline]
            fn $op(self, rhs: T) -> $Lhs {
                Typed::wrap($Op::$op(&self.inner, rhs))
            }
        }

        impl<$($generics)*> $Op<T> for $Lhs {
            type Output = $Lhs;

            #[inline]
            fn $op(mut self, rhs: T) -> $Lhs {
                $OpAssign::$op_assign(&mut self, rhs);
                self
            }
        }
    };
}

scalar_operator!([T: Coefficient, M: Metric, K: Kind] Typed<T, M, K>, Mul::mul, MulAssign::mul_assign);
scalar_operator!([T: Coefficient, M: Metric, K: Kind] Typed<T, M, K>, Div::div, DivAssign::div_assign);
scalar_operator!([T: Coefficient, M: Metric] Typed<T, M, Even>, Add::add, AddAssign::add_assign);
scalar_operator!([T: Coefficient, M: Metric] Typed<T, M, Even>, Sub::sub, SubAssign::sub_assign);

impl<T: Coefficient, M: Metric, K: Kind> Neg for &Typed<T, M, K> {
    type Output = Typed<T, M, K>;

    fn neg(self) -> Typed<T, M, K> {
        self.negated()
    }
}

impl<T: Coefficient, M: Metric, K: Kind> Neg for Typed<T, M, K> {
    type Output = Self;

    fn neg(mut self) -> Self {
        self.negate();
        self
    }
}

impl<T: Coefficient, M: Metric, K: KindMul<DimensionParity<M>>> Not for &Typed<T, M, K> {
    type Output = Typed<T, M, Product<K, DimensionParity<M>>>;

    /// The [dual](Typed::dual).
    fn not(self) -> Self::Output {
        self.dual()
    }
}

impl<T: Coefficient, M: Metric, K: KindMul<DimensionParity<M>>> Not for Typed<T, M, K> {
    type Output = Typed<T, M, Product<K, DimensionParity<M>>>;

    /// The [dual](Typed::dual).
    fn not(self) -> Self::Output {
        self.dual()
    }
}
