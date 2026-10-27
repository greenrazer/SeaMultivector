//! Algebraic identities checked on random multivectors in several signatures.
//!
//! Every sign in a product comes from a mask table, so a wrong mask bit shows up as a broken
//! identity for some blade pair in some signature. Random values with random sparsity reach those
//! pairs without hand-written golden values.

mod common;

use common::{Dense, DenseArray, assert_close};
use proptest::prelude::*;
use seagrass_math::prelude::*;

/// Products of three values with coefficients up to 2 accumulate rounding of about 1e-14 relative
/// to their largest coefficient; 1e-12 leaves room without hiding a sign error (which is O(1)).
const TOL: f64 = 1e-12;

/// Formats coefficients (index = blade bitmask) in the parse syntax, skipping zeros.
fn to_source(signature: &str, coefficients: &[f64]) -> String {
    let mut out = format!("{signature}>");
    for (blade, c) in coefficients.iter().enumerate() {
        if *c == 0.0 {
            continue;
        }
        let sign = if *c < 0.0 { "-" } else { "+" };
        let vectors: String = (0..6)
            .filter(|k| blade >> k & 1 == 1)
            .map(|k| (b'1' + k as u8) as char)
            .collect();
        let name = if vectors.is_empty() {
            String::new()
        } else {
            format!("e{vectors}")
        };
        out.push_str(&format!(" {sign} {}{name}", c.abs()));
    }
    out
}

/// Random coefficients for every blade of a `dim`-dimensional algebra; about half are zero.
fn sparse(dim: u32) -> impl Strategy<Value = Vec<f64>> {
    prop::collection::vec(prop_oneof![Just(0.0), -2.0..2.0_f64], 1 << dim)
}

/// Random coefficients for the grade-1 blades only.
fn vector(dim: u32) -> impl Strategy<Value = Vec<f64>> {
    prop::collection::vec(-2.0..2.0_f64, dim as usize).prop_map(move |v| {
        let mut out = vec![0.0; 1 << dim];
        for (k, c) in v.into_iter().enumerate() {
            out[1 << k] = c;
        }
        out
    })
}

macro_rules! identities {
    ($name:ident, $ty:ty, $dim:literal, $p:literal, $q:literal, $r:literal) => {
        mod $name {
            use super::*;

            const SIGNATURE: &str = concat!($p, ",", $q, ",", $r);

            fn mv(coefficients: &[f64]) -> $ty {
                <$ty>::from_str(&to_source(SIGNATURE, coefficients)).unwrap()
            }

            fn basis_blade(blade: usize) -> $ty {
                let mut coefficients = vec![0.0; 1 << $dim];
                coefficients[blade] = 1.0;
                mv(&coefficients)
            }

            /// The right complement r of x satisfies r ^ x = I and the left complement l satisfies
            /// x ^ l = I (the two coincide in odd dimensions).
            #[test]
            fn complements_wedge_with_their_blade_to_the_pseudoscalar() {
                let pseudoscalar = basis_blade((1 << $dim) - 1);
                for blade in 0..1 << $dim {
                    let b = basis_blade(blade);
                    assert_eq!(&b.right_complemented() ^ &b, pseudoscalar, "blade {blade}");
                    assert_eq!(&b ^ &b.left_complemented(), pseudoscalar, "blade {blade}");
                }
            }

            /// Without null vectors, undual(x) = x I; with them it is the inverse of the right
            /// complement, which is what `dual` returns.
            #[test]
            fn undual_matches_its_definition() {
                let pseudoscalar = basis_blade((1 << $dim) - 1);
                for blade in 0..1 << $dim {
                    let b = basis_blade(blade);
                    if $r == 0 {
                        assert_close(&b.undual(), &(&b * &pseudoscalar), TOL);
                    } else {
                        assert_eq!(b.dual(), b.right_complemented(), "blade {blade}");
                    }
                }
            }

            proptest! {
                #[test]
                fn geometric_product_is_associative(a in sparse($dim), b in sparse($dim), c in sparse($dim)) {
                    let (a, b, c) = (mv(&a), mv(&b), mv(&c));
                    assert_close(&(&(&a * &b) * &c), &(&a * &(&b * &c)), TOL);
                }

                #[test]
                fn reverse_of_product_is_product_of_reverses_swapped(a in sparse($dim), b in sparse($dim)) {
                    let (a, b) = (mv(&a), mv(&b));
                    let lhs = (&a * &b).reversed();
                    let rhs = &b.reversed() * &a.reversed();
                    assert_close(&lhs, &rhs, TOL);
                }

                #[test]
                fn undual_inverts_dual(a in sparse($dim)) {
                    let a = mv(&a);
                    prop_assert_eq!(a.dual().undual(), a);
                }

                /// Every stored blade is written, zeros included, and floats print exactly.
                #[test]
                fn display_round_trips_through_parsing(a in sparse($dim), zero in 0..1usize << $dim) {
                    let mut a = mv(&a);
                    a += &<$ty>::from_blades(a.metric(), [(BladeIndex::new(zero as u8), 0.0)]).unwrap();
                    prop_assert_eq!(a.to_string().parse::<$ty>().unwrap(), a);
                }

                #[test]
                fn vector_squares_to_quadratic_form(v in vector($dim)) {
                    let square = &mv(&v) * &mv(&v);
                    let form: f64 = (0..$dim)
                        .map(|k| {
                            let c = v[1 << k];
                            if k < $p { c * c } else if k < $p + $q { -c * c } else { 0.0 }
                        })
                        .sum();
                    let mut expected = [0.0; 64];
                    expected[0] = form;
                    assert_close(&DenseArray(square.dense()), &DenseArray(expected), TOL);
                }
            }
        }
    };
}

identities!(euclidean_3, Multivector<f64, Dyn3>, 3, 3, 0, 0);
identities!(pga_2, Multivector<f64, Dyn3>, 3, 2, 0, 1);
identities!(mixed_1_2, Multivector<f64, Dyn3>, 3, 1, 2, 0);
identities!(pga_3, Multivector<f64, Dyn4>, 4, 3, 0, 1);
identities!(spacetime, Multivector<f64, Dyn4>, 4, 1, 3, 0);
identities!(mixed_2_1_1, Multivector<f64, Dyn4>, 4, 2, 1, 1);
identities!(conformal_3, Multivector<f64, Dyn5>, 5, 4, 1, 0);
identities!(mixed_3_1_1, Multivector<f64, Dyn5>, 5, 3, 1, 1);
identities!(two_null_3_0_2, Multivector<f64, Dyn5>, 5, 3, 0, 2);
identities!(split_3_3, Multivector<f64, Dyn6>, 6, 3, 3, 0);
identities!(mixed_4_1_1, Multivector<f64, Dyn6>, 6, 4, 1, 1);
identities!(euclidean_6, Multivector<f64, Dyn6>, 6, 6, 0, 0);
