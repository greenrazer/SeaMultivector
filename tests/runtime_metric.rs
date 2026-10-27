//! Behaviour of runtime-metric multivectors (`Multivector<T, Dyn3>` to `Dyn6`): iteration,
//! duals, exponentials, norms, products, addition, scalar operations, and the panic when two
//! operands have different metrics.

use std::fmt::Debug;

use seagrass_math::prelude::*;

type Mv3 = Multivector<f32, Dyn3>;
type Mv4 = Multivector<f32, Dyn4>;
type Mv5 = Multivector<f32, Dyn5>;
type Mv6 = Multivector<f32, Dyn6>;

/// Blades of the 4D fixture with every blade set; prefix a metric before parsing.
const DENSE_4D: &str = "1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e34 + 13e134 + 14e234 + 15e1234";

/// Blades of the left 6D operand shared by the product and scalar-operation tests.
const A_6D: &str = " + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345 +e6 + 12e126 + 8e46";

/// Blades of the right 6D operand shared by the product tests.
const B_6D: &str =
    "  + 3e13 + 1e4 + 2e124  +4e345 +83e1235 + 70e1245 + 4e1345 + 8e2345 + 12e123456";

/// Parses a multivector fixture, panicking on invalid input.
fn mv<M>(s: &str) -> M
where
    M: FromStr,
    M::Err: Debug,
{
    M::from_str(s).unwrap()
}

/// Stored blades as `(blade index, coefficient)` pairs, in iteration order.
///
/// Accepts anything whose `into_iter()` yields `(BladeIndex, &T)` (`&mv`, `mv.grade_iter(..)`).
macro_rules! blades {
    ($blades:expr) => {
        $blades
            .into_iter()
            .map(|(index, c)| (index.bits(), *c))
            .collect::<Vec<_>>()
    };
}

/// Expected `blades!` output, written as blade strings such as `"3e12"`.
fn expected_blades(blades: &[&str]) -> Vec<(u8, f32)> {
    blades
        .iter()
        .map(|s| {
            let blade = Multivector::<f32, Dyn6>::from_str(&format!("6>{s}")).unwrap();
            blades!(&blade)[0]
        })
        .collect()
}

/// Asserts both blade lists store the same blades, in the same order, and that every pair of
/// coefficients differs by strictly less than `tol`.
fn assert_blades_near(actual: &[(u8, f32)], expected: &[(u8, f32)], tol: f32) {
    let actual_indices: Vec<u8> = actual.iter().map(|&(index, _)| index).collect();
    let expected_indices: Vec<u8> = expected.iter().map(|&(index, _)| index).collect();
    assert_eq!(
        actual_indices, expected_indices,
        "blade sets differ\n  left: {actual:?}\n right: {expected:?}"
    );
    for (&(index, a), &(_, e)) in actual.iter().zip(expected) {
        assert!(
            (a - e).abs() < tol,
            "blade {index}: |{a} - {e}| >= {tol}\n  left: {actual:?}\n right: {expected:?}"
        );
    }
}

/// Tests ported from `src/lib.rs`.
mod lib {
    use super::*;

    #[test]
    fn products_of_4d_multivectors_match_expected() {
        let a: Mv4 = mv("3,1,0>e1 + e12");
        let b: Mv4 = mv("3,1,0>e3 + e123 + e4");
        let cases = [
            ("a * b", &a * &b, "3,1,0>-e3 +e13 +e14 +e23 +e123 +e124"),
            ("a | b", &a | &b, "3,1,0>-e3 +e23"),
            ("a ^ b", &a ^ &b, "3,1,0>e13 +e14 +e123 +e124"),
            ("a << b", &a << &b, "3,1,0>-e3 +e23"),
            ("a >> b", &a >> &b, "3,1,0>"),
        ];
        for (operation, actual, expected) in cases {
            assert_eq!(actual, mv::<Mv4>(expected), "{operation}");
        }
    }
}

/// Tests ported from `src/multivector/low_dimensional/custom_metric.rs`.
mod custom_metric {
    use super::*;

    #[test]
    fn iter_yields_every_stored_blade_in_index_order() {
        let a: Mv4 = mv(&format!("4>{DENSE_4D}"));
        assert_eq!(
            blades!(a.iter()),
            expected_blades(&[
                "1", "1e1", "2e2", "3e12", "4e3", "5e13", "6e23", "7e123", "8e4", "9e14", "10e24",
                "11e124", "12e34", "13e134", "14e234", "15e1234",
            ])
        );
    }

    #[test]
    fn grade_iter_yields_only_blades_of_that_grade() {
        let a: Mv4 = mv(&format!("4>{DENSE_4D}"));
        let cases: [(u32, &[&str]); 5] = [
            (0, &["1"]),
            (1, &["1e1", "2e2", "4e3", "8e4"]),
            (2, &["3e12", "5e13", "6e23", "9e14", "10e24", "12e34"]),
            (3, &["7e123", "11e124", "13e134", "14e234"]),
            (4, &["15e1234"]),
        ];
        for (grade, expected) in cases {
            assert_eq!(
                blades!(a.grade_iter(grade)),
                expected_blades(expected),
                "grade {grade}"
            );
        }
    }

    #[test]
    fn grade_iter_mut_writes_only_blades_of_that_grade() {
        let cases = [
            (
                format!("4>{DENSE_4D}"),
                1.12334,
                "4>1 + e1 + 2e2 + 1.12334e12 + 4e3 + 1.12334e13 + 1.12334e23 + 7e123 + 8e4 + 1.12334e14 + 1.12334e24 + 11e124 + 1.12334e34 + 13e134 + 14e234 + 15e1234",
            ),
            (
                "4> e1 + 3e12 + 4e3 + 13e134 + 15e1234".to_string(),
                4.3342,
                "4> e1 + 4.3342e12 + 4e3 + 13e134 + 15e1234",
            ),
        ];
        for (input, value, expected) in cases {
            let mut a: Mv4 = mv(&input);
            for blade in a.grade_iter_mut(2) {
                *blade.1 = value;
            }
            assert_eq!(a, mv::<Mv4>(expected), "{input}");
        }
    }

    #[test]
    fn not_operator_returns_dual_by_reference_and_by_value() {
        let a: Mv5 = mv(
            "2,3>1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e34 + 13e134 + 14e234 + 15e1234 + e5 + e15 + 2e25 + 3e125 + 4e35 + 5e135 + 6e235 + 7e1235 + 8e45 + 9e145 + 10e245 + 11e1245 + 12e345 + 13e1345 + 14e2345 + 15e12345",
        );
        let out_dual: Mv5 = mv(
            "2,3>15 +14e1 -13e2 -11e3 +7e4 -15e5 -12e12 -10e13 +6e14 -14e15 +9e23 -5e24 +13e25 -3e34 +11e35 -7e45 +8e123 -4e124 +12e125 -2e134 +10e135 -6e145 +e234 -9e235 +5e245 +3e345 +e1234 -8e1235 +4e1245 +2e1345 -e2345 -e12345",
        );
        assert_eq!(!&a, out_dual);
        assert_eq!(!a, out_dual);
    }

    #[test]
    fn exp_of_bivector_matches_expected() {
        let cases = [
            (
                blades!(&mv::<Mv3>("2,1>e12 + 2e13 + 3e23").exp()),
                blades!(&mv::<Mv3>(
                    "2,1>15.989522934 +4.6067419052e12 +9.2134838104e13 +13.8202257156e23"
                )),
                1e-6,
            ),
            (
                blades!(&mv::<Mv5>("4,1>-e15 +4e25+e45 -e35").exp()),
                blades!(&mv::<Mv5>(
                    "4,1>39.0919036865 -8.9653635025e15 +35.86145401e25 -8.9653635025e35 +8.9653635025e45"
                )),
                1e-5,
            ),
            (
                blades!(&mv::<Mv4>("3,0,1>e12 + 2e13 + 3e23").exp()),
                blades!(&mv::<Mv4>(
                    "3,0,1>-0.8252990842 -0.1509213299e12 -0.3018426597e13 -0.4527639747e23"
                )),
                1e-6,
            ),
            (
                blades!(&mv::<Mv4>("3,0,1>4e14 + 5e24 + 6e34").exp()),
                blades!(&mv::<Mv4>("3,0,1>1 +4e14 +5e24 +6e34")),
                1e-6,
            ),
        ];
        for (actual, expected, tol) in cases {
            assert_blades_near(&actual, &expected, tol);
        }
    }

    mod norm {
        use super::*;

        #[test]
        fn scalar_norm_depends_on_metric() {
            let cases = [
                (format!("3,0,1>{DENSE_4D}"), 11.874_342),
                (format!("3,1,0>{DENSE_4D}"), 30.967_726),
            ];
            for (input, expected) in cases {
                assert_eq!(mv::<Mv4>(&input).norm(), expected, "{input}");
            }
        }

        #[test]
        fn scalar_normalized_returns_normalized_copy() {
            let a: Mv4 = mv(&format!("3,0,1>{DENSE_4D}"));
            let out: Mv4 = mv(
                "3,0,1>0.084215194 +0.67372155e4 +0.084215194e1 +0.168430388e2 +0.3368607759e3 +0.7579367e14 +0.84215194e24 +1.0105823e34 +0.252645582e12 +0.4210759699e13 +0.5052911639e23 +0.9263671e124 +1.0947975e134 +1.1790127e234 +0.5895063e123 +1.2632279e1234",
            );
            assert_blades_near(&blades!(&a.normalized()), &blades!(&out), 1e-7);
        }

        #[test]
        fn scalar_normalize_normalizes_in_place() {
            let mut a: Mv4 = mv(&format!("3,1,0>{DENSE_4D}"));
            let out: Mv4 = mv(
                "3,1,0>0.0322916843 +0.0322916843e1 +0.0645833686e2 +0.1291667372e3 +0.2583334744e4 +0.0968750492e12 +0.1614584178e13 +0.290625155e14 +0.1937500983e23 +0.3229168355e24 +0.3875001967e34 +0.22604178e123 +0.3552085161e124 +0.41979188e134 +0.4520835876e234 +0.48437524e1234",
            );
            a.normalize();
            assert_blades_near(&blades!(&a), &blades!(&out), 1e-7);
        }
    }
}

/// Tests ported from `src/multivector/low_dimensional/custom_metric/linear_operations/`.
mod linear_operations {
    use super::*;

    mod geometric_product {
        use super::*;

        #[test]
        fn method_and_operator_forms_match_expected() {
            let a: Mv6 = mv(&format!("3,3>{A_6D}"));
            let b: Mv6 = mv(&format!("3,3>{B_6D}"));
            let out: Mv6 = mv(
                "3,3>115+115e1-27e2+1053e3+1240e4+189e5+188e6-21e12+983e13+1159e14+202e15+168e16-919e23-1075e24-215e25-156e26+17e34-351e35-301e45+23e46-846e123-991e124-228e125-128e126+17e134-265e135+3e136-228e145-17e234+143e235-36e236+135e245+303e345-1160e356-924e456-15e1234+230e1235+124e1245-14e1246-560e1256+442e1345+24e1346-152e1356-72e1456-393e2345+44e2356+60e2456-40e3456-360e12345+179e12356+118e12456+76e13456-28e23456-700e123456",
            );
            let results = [
                ("a.geometric_product(&b)", a.geometric_product(&b)),
                ("&a * &b", &a * &b),
                ("a * &b", a.clone() * &b),
                ("&a * b", &a * b.clone()),
                ("a * b", a * b),
            ];
            for (form, result) in results {
                assert_eq!(result, out, "{form}");
            }
        }

        #[test]
        #[should_panic(expected = "assertion `left == right` failed")]
        fn panics_on_mismatched_metrics() {
            let a: Mv6 = mv(&format!("1,2,3>{A_6D}"));
            let b: Mv6 = mv(&format!("3,2,1>{B_6D}"));
            let _ = a.geometric_product(&b);
        }
    }

    mod inner_product {
        use super::*;

        #[test]
        fn method_and_operator_forms_match_expected() {
            let a: Mv6 = mv(&format!("3,3>{A_6D}"));
            let b: Mv6 = mv(&format!("3,3>{B_6D}"));
            // The trailing zero-width space is in the original fixture; the parser ignores it.
            let out: Mv6 = mv(
                "3,3>115+115e1 -27e2 +1053e3 +1240e4 +189e5 +188e6 -87e12 +3e13 -4e14 +202e15 +168e16 +2e24 -215e25 -156e26 -351e35 -301e45 +2e124 -228e125 -144e126 -185e135 -156e145 +33e235 +83e245 -120e345 -132e356 -84e456 +194e1235 +70e1245 +4e1345 -120e1356 -72e1456 +8e2345 +108e2356 +60e2456 -36e3456 +12e12345 +96e12356 +48e12456 -24e13456 +12e23456 +12e123456",
            );
            let results = [
                ("a.inner_product(&b)", a.inner_product(&b)),
                ("&a | &b", &a | &b),
                ("a | &b", a.clone() | &b),
                ("&a | b", &a | b.clone()),
                ("a | b", a | b),
            ];
            for (form, result) in results {
                assert_eq!(result, out, "{form}");
            }
        }

        #[test]
        #[should_panic(expected = "assertion `left == right` failed")]
        fn panics_on_mismatched_metrics() {
            let a: Mv6 = mv(&format!("1,2,3>{A_6D}"));
            let b: Mv6 = mv(&format!("3,2,1>{B_6D}"));
            let _ = a.inner_product(&b);
        }
    }

    mod left_contraction {
        use super::*;

        #[test]
        fn method_and_operator_forms_match_expected() {
            let a: Mv6 = mv(&format!("3,3>{A_6D}"));
            let b: Mv6 = mv(&format!("3,3>{B_6D}"));
            let out: Mv6 = mv(
                "3,3>115+56e1+78e2+3e3-5e4+189e5+180e6-16e12+3e13-4e14+202e15+168e16+2e24-215e25-156e26-333e35-262e45+2e124-228e125-144e126-198e135-156e145+19e235+38e245-120e345-132e356-84e456+179e1235+70e1245+4e1345-120e1356-72e1456+8e2345+108e2356+60e2456-36e3456+12e12345+96e12356+48e12456-24e13456+12e23456+12e123456",
            );
            let results = [
                ("a.left_contraction(&b)", a.left_contraction(&b)),
                ("&a << &b", &a << &b),
                ("a << &b", a.clone() << &b),
                ("&a << b", &a << b.clone()),
                ("a << b", a << b),
            ];
            for (form, result) in results {
                assert_eq!(result, out, "{form}");
            }
        }

        #[test]
        #[should_panic(expected = "assertion `left == right` failed")]
        fn panics_on_mismatched_metrics() {
            let a: Mv6 = mv(&format!("1,2,3>{A_6D}"));
            let b: Mv6 = mv(&format!("3,2,1>{B_6D}"));
            let _ = a.left_contraction(&b);
        }
    }

    mod outer_product {
        use super::*;

        #[test]
        fn method_and_operator_forms_match_expected() {
            let a: Mv6 = mv(&format!("3,3>{A_6D}"));
            let b: Mv6 = mv(&format!("3,3>{B_6D}"));
            let out: Mv6 = mv(
                "3,3>e4 +3e13 +e14 +2e24 +4e34 -e46 -6e123 +5e124 +29e134 +3e136 +6e234 +4e345 -15e1234 +83e1235 +70e1245 -14e1246 +8e1345 +24e1346 +16e2345 -4e3456 -372e12345 +83e12356 +70e12456 +4e13456 +8e23456 -700e123456",
            );
            let results = [
                ("a.outer_product(&b)", a.outer_product(&b)),
                ("&a ^ &b", &a ^ &b),
                ("a ^ &b", a.clone() ^ &b),
                ("&a ^ b", &a ^ b.clone()),
                ("a ^ b", a ^ b),
            ];
            for (form, result) in results {
                assert_eq!(result, out, "{form}");
            }
        }
    }

    mod regressive_product {
        use super::*;

        #[test]
        fn method_and_operator_forms_match_expected() {
            let a: Mv5 = mv(
                "3,2> + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345",
            );
            let b: Mv5 = mv("3,2>  + 3e13 + 1e4 + 2e124  +4e345");
            let out: Mv5 = mv("3,2> 12 +70e3 +83e4 +45e13 +26e14 +28e24 +30e124 +60e345");
            let results = [
                ("a.regressive_product(&b)", a.regressive_product(&b)),
                ("&a & &b", &a & &b),
                ("a & &b", a.clone() & &b),
                ("&a & b", &a & b.clone()),
                ("a & b", a & b),
            ];
            for (form, result) in results {
                assert_eq!(result, out, "{form}");
            }
        }
    }

    mod right_contraction {
        use super::*;

        #[test]
        fn method_and_operator_forms_match_expected() {
            let a: Mv6 = mv(&format!("3,3>{A_6D}"));
            let b: Mv6 = mv(&format!("3,3>{B_6D}"));
            let out: Mv6 = mv(
                "3,3>115+59e1-105e2+1050e3+1245e4+8e6-71e12-18e35-39e45+13e135+14e235+45e245+15e1235",
            );
            let results = [
                ("a.right_contraction(&b)", a.right_contraction(&b)),
                ("&a >> &b", &a >> &b),
                ("&a >> b", &a >> b.clone()),
                ("a >> &b", a.clone() >> &b),
                ("a >> b", a >> b),
            ];
            for (form, result) in results {
                assert_eq!(result, out, "{form}");
            }
        }

        #[test]
        #[should_panic(expected = "assertion `left == right` failed")]
        fn panics_on_mismatched_metrics() {
            let a: Mv6 = mv(&format!("1,2,3>{A_6D}"));
            let b: Mv6 = mv(&format!("3,2,1>{B_6D}"));
            let _ = a.right_contraction(&b);
        }
    }

    mod scalar_product {
        use super::*;

        #[test]
        fn scalar_product_of_degenerate_5d_multivectors() {
            let a: Mv5 = mv(
                "3,1,1> + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345",
            );
            let b: Mv5 = mv("3,1,1> + 3e13 + 1e4 + 2e124  +4e345");
            assert_eq!(a.scalar_product(&b), -1.0);
        }
    }
}

/// Tests ported from
/// `src/multivector/low_dimensional/custom_metric/non_geometric_linear_operations/addition.rs`.
mod non_geometric_linear_operations {
    use super::*;

    mod addition {
        use super::*;

        #[test]
        fn add_assign_basis_blade_inserts_missing_blades() {
            let expected: Mv3 = mv("3>-1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123");
            let cases: [(&str, &[&str]); 2] = [
                ("3>-1 + e1 + 2e2 + 3e12 + 5e13 + 6e23 + 7e123", &["4e3"]),
                ("3>-1 + e1 + 2e2 + 3e12 + 5e13", &["4e3", "7e123", "6e23"]),
            ];
            for (input, added) in cases {
                let mut a: Mv3 = mv(input);
                for blade in added {
                    a += &mv::<Mv3>(&format!("3>{blade}"));
                }
                assert_eq!(a, expected, "{input} += {added:?}");
            }
        }

        #[test]
        fn adding_multivectors_sums_coefficients() {
            let mut a: Mv3 = mv("3>-1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123");
            let b: Mv3 = mv("3>-1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123");
            let out: Mv3 = mv("3>-2 + 2e1 + 4e2 + 6e12 + 8e3 + 10e13 + 12e23 + 14e123");
            assert_eq!(&a + &b, out);
            a += &b;
            assert_eq!(a, out);
        }
    }
}

/// Tests ported from `src/multivector/low_dimensional/custom_metric/scalar_operations/`.
mod scalar_operations {
    use super::*;

    mod addition {
        use super::*;

        #[test]
        fn adding_scalar_changes_only_scalar_blade() {
            let mut a: Mv6 = mv(&format!("3,3>{A_6D}"));
            let out: Mv6 = mv(
                "3,3> + 3 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345 +e6 + 12e126 + 8e46",
            );
            assert_eq!(&a + 2.0, out);
            a += 2.0;
            assert_eq!(a, out);
        }
    }

    mod division {
        use super::*;

        #[test]
        fn dividing_by_scalar_scales_every_blade() {
            let mut a: Mv6 = mv(&format!("3,3>{A_6D}"));
            let out: Mv6 = mv(
                "3,3> +0.5 + 0.5e1 + 1e2 + 1.5e12 + 2e3 + 2.5e13 + 3e23 + 3.5e123 + 4e4 + 4.5e14 + 5e24 + 5.5e124 + 6e345 + 6.5e1345 + 7e2345 + 7.5e12345 +0.5e6 + 6e126 + 4e46",
            );
            assert_eq!(&a / 2.0, out);
            a /= 2.0;
            assert_eq!(a, out);
        }
    }

    mod multiplication {
        use super::*;

        #[test]
        fn multiplying_by_scalar_scales_every_blade() {
            let mut a: Mv6 = mv(&format!("3,3>{A_6D}"));
            let out: Mv6 = mv(
                "3,3> +2 + 2e1 + 4e2 + 6e12 + 8e3 + 10e13 + 12e23 + 14e123 + 16e4 + 18e14 + 20e24 + 22e124 + 24e345 + 26e1345 + 28e2345 + 30e12345 +2e6 + 24e126 + 16e46",
            );
            assert_eq!(&a * 2.0, out);
            a *= 2.0;
            assert_eq!(a, out);
        }
    }

    mod subtraction {
        use super::*;

        #[test]
        fn subtracting_scalar_changes_only_scalar_blade() {
            let mut a: Mv6 = mv(&format!("3,3>{A_6D}"));
            let out: Mv6 = mv(
                "3,3> -1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345 +e6 + 12e126 + 8e46",
            );
            assert_eq!(&a - 2.0, out);
            a -= 2.0;
            assert_eq!(a, out);
        }
    }
}
