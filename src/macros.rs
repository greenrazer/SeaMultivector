//! Macros that generate the owned and borrowed forms of the operators.

/// Implements `Op` for the four combinations of owned and borrowed `Lhs` and `Rhs`.
///
/// The body is the `&Lhs op &Rhs` form, with the operands bound to the given names; the other
/// three borrow their operands and forward to it.
macro_rules! binary_operator {
    (
        [$($generics:tt)*] $Op:ident::$op:ident(
            $lhs:ident: $Lhs:ty, $rhs:ident: $Rhs:ty
        ) -> $Out:ty $body:block
    ) => {
        impl<$($generics)*> $Op<&$Rhs> for &$Lhs {
            type Output = $Out;

            #[inline]
            fn $op(self, $rhs: &$Rhs) -> $Out {
                let $lhs = self;
                $body
            }
        }

        impl<$($generics)*> $Op<$Rhs> for &$Lhs {
            type Output = $Out;

            #[inline]
            fn $op(self, rhs: $Rhs) -> $Out {
                $Op::$op(self, &rhs)
            }
        }

        impl<$($generics)*> $Op<&$Rhs> for $Lhs {
            type Output = $Out;

            #[inline]
            fn $op(self, rhs: &$Rhs) -> $Out {
                $Op::$op(&self, rhs)
            }
        }

        impl<$($generics)*> $Op<$Rhs> for $Lhs {
            type Output = $Out;

            #[inline]
            fn $op(self, rhs: $Rhs) -> $Out {
                $Op::$op(&self, &rhs)
            }
        }
    };
}

/// Implements `OpAssign<Rhs>` by forwarding to `OpAssign<&Rhs>`.
macro_rules! forward_assign {
    ([$($generics:tt)*] $OpAssign:ident::$op_assign:ident($Lhs:ty, $Rhs:ty)) => {
        impl<$($generics)*> $OpAssign<$Rhs> for $Lhs {
            #[inline]
            fn $op_assign(&mut self, rhs: $Rhs) {
                $OpAssign::$op_assign(self, &rhs);
            }
        }
    };
}

pub(crate) use binary_operator;
pub(crate) use forward_assign;
