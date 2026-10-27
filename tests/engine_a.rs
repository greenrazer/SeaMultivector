//! Golden values from the 0.1 engine's unit tests: parsing, conjugations, duals and
//! complements, iteration, cleanup, and the geometric, inner, left-contraction and outer
//! products in metrics from 3D to 6D. Each module is named after the 0.1 file its cases came
//! from. Products compare with `assert_eq!`, so the stored blades of every result are checked
//! along with the coefficients.

use std::fmt::Debug;

use seagrass_math::prelude::*;

type Mv3 = Multivector<f32, Dyn3>;
type Mv4 = Multivector<f32, Dyn4>;
type Mv5 = Multivector<f32, Dyn5>;
type Mv6 = Multivector<f32, Dyn6>;

/// Operands shared verbatim by several of the original product tests.
const A_5D: &str = " + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345";
const B_5D: &str = "  + 3e13 + 1e4 + 2e124  +4e345";
const A_6D: &str = " + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345 +e6 + 12e126 + 8e46";
const B_6D: &str =
    "  + 3e13 + 1e4 + 2e124  +4e345 +83e1235 + 70e1245 + 4e1345 + 8e2345 + 12e123456";

/// Parses `terms` in the algebra given by `metric` (`"p"`, `"p,q"` or `"p,q,r"`).
fn parse<M>(metric: &str, terms: &str) -> M
where
    M: FromStr,
    M::Err: Debug,
{
    let source = format!("{metric}>{terms}");
    M::from_str(&source).unwrap_or_else(|e| panic!("failed to parse {source:?}: {e:?}"))
}

/// Asserts structural equality (same stored blades, same coefficients) of `op(input)`
/// and `expected`.
fn assert_unary<M>(metric: &str, input: &str, expected: &str, op: impl Fn(&M) -> M)
where
    M: FromStr + PartialEq + Debug,
    M::Err: Debug,
{
    let actual = op(&parse::<M>(metric, input));
    assert_eq!(
        actual,
        parse::<M>(metric, expected),
        "metric {metric}, input {input:?}"
    );
}

/// Asserts structural equality of `op(lhs, rhs)` and `expected`.
fn assert_binary<M>(metric: &str, lhs: &str, rhs: &str, expected: &str, op: impl Fn(&M, &M) -> M)
where
    M: FromStr + PartialEq + Debug,
    M::Err: Debug,
{
    let actual = op(&parse::<M>(metric, lhs), &parse::<M>(metric, rhs));
    assert_eq!(
        actual,
        parse::<M>(metric, expected),
        "metric {metric}, lhs {lhs:?}, rhs {rhs:?}"
    );
}

/// Asserts structural equality of `input` after the in-place `op` and `expected`.
fn assert_in_place<M>(metric: &str, input: &str, expected: &str, op: impl FnOnce(&mut M))
where
    M: FromStr + PartialEq + Debug,
    M::Err: Debug,
{
    let mut actual = parse::<M>(metric, input);
    op(&mut actual);
    assert_eq!(
        actual,
        parse::<M>(metric, expected),
        "metric {metric}, input {input:?}"
    );
}

/// Collects `(blade index, coefficient)` pairs in iteration order. The blade index is the
/// bitmask of its basis vectors: e1 = 1, e2 = 2, e3 = 4, e4 = 8, so e12 = 3, e134 = 13.
fn blades<'a>(items: impl IntoIterator<Item = (BladeIndex, &'a f32)>) -> Vec<(u8, f32)> {
    items
        .into_iter()
        .map(|(index, c)| (index.bits(), *c))
        .collect()
}

/// Parses single-blade strings such as `"3e12"` into `(blade index, coefficient)` pairs.
fn expected_blades(terms: &[&str]) -> Vec<(u8, f32)> {
    terms
        .iter()
        .map(|term| {
            let blade = Multivector::<f32, Dyn6>::from_str(&format!("6>{term}")).unwrap();
            blades(&blade)[0]
        })
        .collect()
}

/// Ported from `src/multivector/low_dimensional/core.rs`.
mod core_multivector {
    use super::*;

    const FULL_3D: &str = "1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123";
    const FULL_4D: &str = "1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e34 + 13e134 + 14e234 + 15e1234";
    const SPARSE_4D: &str = " e1 + 3e12 + 4e3 + 13e134 + 15e1234";

    #[test]
    fn parses_terms_into_canonical_blade_order() {
        // basis_blade_set 255: all eight blades of the 3D algebra are stored, in index order.
        let a: Mv3 = parse("3", "1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123");
        let coefficients = [1_f32, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0];
        assert_eq!(blades(&a), (0..8).zip(coefficients).collect::<Vec<_>>());

        // basis_blade_set 65535: all sixteen blades of the 4D algebra are stored, in index order.
        let a: Mv4 = parse(
            "4",
            "1 - e1 + 2e2 - 3e12 + 4e3 - 5e13 + 6e23 - 7e123 + 8e4 - 9e14 + 10e24 - 11e124 + 12e34 - 13e134 + 14e234 - 15e1234",
        );
        let coefficients = [
            1_f32, -1.0, 2.0, -3.0, 4.0, -5.0, 6.0, -7.0, 8.0, -9.0, 10.0, -11.0, 12.0, -13.0,
            14.0, -15.0,
        ];
        assert_eq!(blades(&a), (0..16).zip(coefficients).collect::<Vec<_>>());
    }

    #[test]
    fn reports_specific_parse_errors() {
        let cases = [
            ("e4", ParseError::InvalidBasisBladeForAlgebra),
            ("e21", ParseError::InvalidBasisVectorOrdering),
            ("e11", ParseError::DuplicateBasisVector),
            ("e1 + e1", ParseError::DuplicateBasisBlade),
        ];
        for (terms, error) in cases {
            assert_eq!(
                Mv3::from_str(&format!("3>{terms}")).unwrap_err(),
                error,
                "terms {terms:?}"
            );
        }
    }

    #[test]
    fn conjugations_flip_signs_by_grade() {
        assert_unary::<Mv3>(
            "3",
            FULL_3D,
            "1 + e1 + 2e2 - 3e12 + 4e3 - 5e13 - 6e23 - 7e123",
            |a| a.reversed(),
        );
        assert_unary::<Mv3>(
            "3",
            FULL_3D,
            "1 - e1 - 2e2 + 3e12 - 4e3 + 5e13 + 6e23 - 7e123",
            |a| a.involuted(),
        );
        assert_unary::<Mv3>(
            "3",
            FULL_3D,
            "1 - e1 - 2e2 - 3e12 - 4e3 - 5e13 - 6e23 + 7e123",
            |a| a.conjugated(),
        );
    }

    #[test]
    fn poincare_dual_and_complements_match_golden_values() {
        assert_unary::<Mv3>(
            "3",
            FULL_3D,
            "7 + 6e1 + 5e2 + 4e12 + 3e3 + 2e13 + 1e23 + 1e123",
            |a| a.poincare_dual(),
        );
        assert_unary::<Mv3>(
            "3",
            FULL_3D,
            "7 + 6e1 - 5e2 + 4e12 + 3e3 - 2e13 + 1e23 + 1e123",
            |a| a.right_complemented(),
        );
        assert_unary::<Mv4>(
            "4",
            FULL_4D,
            "15 + 14e1 + 13e2 + 12e12 + 11e3 + 10e13 + 9e23 + 8e123 + 7e4 + 6e14 + 5e24 + 4e124 + 3e34 + 2e134 + 1e234 + 1e1234",
            |a| a.poincare_dual(),
        );
        assert_unary::<Mv4>(
            "4",
            FULL_4D,
            "15 + 14e1 - 13e2 + 12e12 + 11e3 - 10e13 + 9e23 + 8e123 - 7e4 + 6e14 - 5e24 - 4e124 + 3e34 + 2e134 - 1e234 + 1e1234",
            |a| a.right_complemented(),
        );
        assert_unary::<Mv4>(
            "4",
            FULL_4D,
            "15 - 14e1 + 13e2 + 12e12 - 11e3 - 10e13 + 9e23 - 8e123 + 7e4 + 6e14 - 5e24 + 4e124 + 3e34 - 2e134 + 1e234 + 1e1234",
            |a| a.left_complemented(),
        );
    }

    #[test]
    fn in_place_conjugations_flip_signs_by_grade() {
        assert_in_place::<Mv3>(
            "3",
            FULL_3D,
            "1 + e1 + 2e2 - 3e12 + 4e3 - 5e13 - 6e23 - 7e123",
            |a| a.reverse(),
        );
        assert_in_place::<Mv3>(
            "3",
            FULL_3D,
            "1 - e1 - 2e2 + 3e12 - 4e3 + 5e13 + 6e23 - 7e123",
            |a| a.involute(),
        );
        assert_in_place::<Mv3>(
            "3",
            FULL_3D,
            "1 - e1 - 2e2 - 3e12 - 4e3 - 5e13 - 6e23 + 7e123",
            |a| a.conjugate(),
        );
    }

    #[test]
    fn in_place_poincare_dual_and_complements_match_golden_values() {
        assert_in_place::<Mv3>(
            "3",
            FULL_3D,
            "7 + 6e1 + 5e2 + 4e12 + 3e3 + 2e13 + 1e23 + 1e123",
            |a| a.poincare_dualize(),
        );
        assert_in_place::<Mv3>(
            "3",
            FULL_3D,
            "7 + 6e1 - 5e2 + 4e12 + 3e3 - 2e13 + 1e23 + 1e123",
            |a| a.right_complement(),
        );
        assert_in_place::<Mv3>(
            "3",
            FULL_3D,
            "7 + 6e1 - 5e2 + 4e12 + 3e3 - 2e13 + 1e23 + 1e123",
            |a| a.left_complement(),
        );
        assert_in_place::<Mv4>(
            "4",
            FULL_4D,
            "15 + 14e1 - 13e2 + 12e12 + 11e3 - 10e13 + 9e23 + 8e123 - 7e4 + 6e14 - 5e24 - 4e124 + 3e34 + 2e134 - 1e234 + 1e1234",
            |a| a.right_complement(),
        );
        assert_in_place::<Mv4>(
            "4",
            FULL_4D,
            "15 - 14e1 + 13e2 + 12e12 - 11e3 - 10e13 + 9e23 - 8e123 + 7e4 + 6e14 - 5e24 + 4e124 + 3e34 - 2e134 + 1e234 + 1e1234",
            |a| a.left_complement(),
        );
    }

    #[test]
    fn iter_yields_stored_blades_in_index_order() {
        let a: Mv4 = parse("4", FULL_4D);
        assert_eq!(
            blades(a.iter()),
            expected_blades(&[
                "1", "1e1", "2e2", "3e12", "4e3", "5e13", "6e23", "7e123", "8e4", "9e14", "10e24",
                "11e124", "12e34", "13e134", "14e234", "15e1234",
            ])
        );

        let a: Mv4 = parse("4", SPARSE_4D);
        assert_eq!(
            blades(a.iter()),
            expected_blades(&["1e1", "3e12", "4e3", "13e134", "15e1234"])
        );
    }

    #[test]
    fn iter_mut_writes_every_stored_coefficient() {
        let mut a: Mv4 = parse("4", FULL_4D);
        for blade in a.iter_mut() {
            *blade.1 = 1.12334;
        }
        assert_eq!(
            a,
            parse(
                "4",
                "1.12334 + 1.12334e1 + 1.12334e2 + 1.12334e12 + 1.12334e3 + 1.12334e13 + 1.12334e23 + 1.12334e123 + 1.12334e4 + 1.12334e14 + 1.12334e24 + 1.12334e124 + 1.12334e34 + 1.12334e134 + 1.12334e234 + 1.12334e1234"
            )
        );

        let mut a: Mv4 = parse("4", SPARSE_4D);
        for blade in a.iter_mut() {
            *blade.1 = 4.3342;
        }
        assert_eq!(
            a,
            parse(
                "4",
                " 4.3342e1 + 4.3342e12 + 4.3342e3 + 4.3342e134 + 4.3342e1234"
            )
        );
    }

    #[test]
    fn grade_iter_yields_only_blades_of_that_grade() {
        let a: Mv4 = parse("4", FULL_4D);
        let cases: [(u32, &[&str]); 5] = [
            (0, &["1"]),
            (1, &["1e1", "2e2", "4e3", "8e4"]),
            (2, &["3e12", "5e13", "6e23", "9e14", "10e24", "12e34"]),
            (3, &["7e123", "11e124", "13e134", "14e234"]),
            (4, &["15e1234"]),
        ];
        for (grade, expected) in cases {
            assert_eq!(
                blades(a.grade_iter(grade)),
                expected_blades(expected),
                "grade {grade}"
            );
        }
    }

    #[test]
    fn grade_iter_mut_writes_only_blades_of_that_grade() {
        let mut a: Mv4 = parse("4", FULL_4D);
        for blade in a.grade_iter_mut(2) {
            *blade.1 = 1.12334;
        }
        assert_eq!(
            a,
            parse(
                "4",
                "1 + e1 + 2e2 + 1.12334e12 + 4e3 + 1.12334e13 + 1.12334e23 + 7e123 + 8e4 + 1.12334e14 + 1.12334e24 + 11e124 + 1.12334e34 + 13e134 + 14e234 + 15e1234"
            )
        );

        let mut a: Mv4 = parse("4", SPARSE_4D);
        for blade in a.grade_iter_mut(2) {
            *blade.1 = 4.3342;
        }
        assert_eq!(a, parse("4", " e1 + 4.3342e12 + 4e3 + 13e134 + 15e1234"));
    }

    #[test]
    fn remove_zeros_drops_exactly_zero_blades() {
        assert_in_place::<Mv4>(
            "4",
            "1 + e1 + 0e2 + 0e12 + 0e3 + 0e13 + 0e23 + 0e123 + 0e4 + 0e14 + 0e24 + 0e124 + 0e34 + 0e134 + 0e234 + 0e1234",
            "1 + e1",
            |a| a.remove_zeros(),
        );
        assert_in_place::<Mv4>(
            "4",
            "0 + 0e1 + 0e2 + 0e12 + 0e3 + 1e13 + 0e23 + 0e123 + 0e4 + 0e14 + 0e24 + 0e124 + 10e34 + 0e134 + 0e234 + 0e1234",
            "1e13 + 10e34",
            |a| a.remove_zeros(),
        );
    }

    #[test]
    fn remove_approx_zeros_drops_blades_below_epsilon() {
        assert_in_place::<Mv4>(
            "4",
            "1 + e1 + 0.0000000000000000001e2 + 0.0000000000000000001e12 + 0.0000000000000000001e3 + 0.0000000000000000001e13 + 0.0000000000000000001e23 + 0.0000000000000000001e123 + 0.0000000000000000001e4 + 0.0000000000000000001e14 + 0.0000000000000000001e24 + 0.0000000000000000001e124 + 0.0000000000000000001e34 + 0.0000000000000000001e134 + 0.0000000000000000001e234 + 0.0000000000000000001e1234",
            "1 + e1",
            |a| a.remove_approx_zeros(1e-15),
        );
        assert_in_place::<Mv4>(
            "4",
            "0.0000000000000000001 + 0.0000000000000000001e1 + 0.0000000000000000001e2 + 0.0000000000000000001e12 + 0.0000000000000000001e3 + 1e13 + 0.0000000000000000001e23 + 0.0000000000000000001e123 + 0.0000000000000000001e4 + 0.0000000000000000001e14 + 0.0000000000000000001e24 + 0.0000000000000000001e124 + 10.0000000000000000001e34 + 0.0000000000000000001e134 + 0.0000000000000000001e234 + 0.0000000000000000001e1234",
            "1e13 + 10e34",
            |a| a.remove_approx_zeros(1e-15),
        );
    }
}

/// Ported from `src/multivector/low_dimensional/core/linear_operations/dual.rs`.
mod dual {
    use super::*;

    const FULL_3D: &str = "1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123";
    const FULL_4D: &str = "1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e34 + 13e134 + 14e234 + 15e1234";
    const FULL_5D: &str = "1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e34 + 13e134 + 14e234 + 15e1234 + e5 + e15 + 2e25 + 3e125 + 4e35 + 5e135 + 6e235 + 7e1235 + 8e45 + 9e145 + 10e245 + 11e1245 + 12e345 + 13e1345 + 14e2345 + 15e12345";
    const FULL_6D: &str = "1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e34 + 13e134 + 14e234 + 15e1234 + e5 + e15 + 2e25 + 3e125 + 4e35 + 5e135 + 6e235 + 7e1235 + 8e45 + 9e145 + 10e245 + 11e1245 + 12e345 + 13e1345 + 14e2345 + 15e12345 + e6 + e16 + 2e26 + 3e126 + 4e36 + 5e136 + 6e236 + 7e1236 + 8e46 + 9e146 + 10e246 + 11e1246 + 12e346 + 13e1346 + 14e2346 + 15e12346 + e56 + e156 + 2e256 + 3e1256 + 4e356 + 5e1356 + 6e2356 + 7e12356 + 8e456 + 9e1456 + 10e2456 + 11e12456 + 12e3456 + 13e13456 + 14e23456 + 15e123456";
    const PARTIAL_5D: &str = " + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345";

    #[test]
    fn wedge_with_own_dual_is_signed_pseudoscalar() {
        assert_unary::<Mv3>("3,0,0", "1e13", "1e123", |a| a ^ &a.dual());
        assert_unary::<Mv4>("4,0,0", "1e13", "-1e1234", |a| a ^ &a.dual());
        assert_unary::<Mv5>("5,0,0", "1e13", "-1e12345", |a| a ^ &a.dual());
        assert_unary::<Mv6>("6,0,0", "1e13", "1e123456", |a| a ^ &a.dual());
        assert_unary::<Mv5>("4,0,1", "1e13", "1e12345", |a| a ^ &a.dual());
        assert_unary::<Mv5>("4,0,1", "1e124", "1e12345", |a| a ^ &a.dual());
    }

    #[test]
    fn matches_golden_values_in_positive_metrics() {
        let m = "3,0,0";
        assert_unary::<Mv3>(
            m,
            FULL_3D,
            "7 + 6e1 - 5e2 - 4e12 + 3e3 + 2e13 - 1e23 - 1e123",
            |a| a.dual(),
        );
        assert_unary::<Mv3>(
            m,
            FULL_3D,
            "-7 - 6e1 + 5e2 + 4e12 - 3e3 - 2e13 + 1e23 + 1e123",
            |a| a.undual(),
        );
        assert_unary::<Mv3>(m, FULL_3D, FULL_3D, |a| a.undual().dual());

        let m = "4,0,0";
        assert_unary::<Mv4>(
            m,
            FULL_4D,
            "15 +14e1 -13e2 +11e3 -7e4 -12e12 +10e13 -6e14 -9e23 +5e24 -3e34 -8e123 +4e124 -2e134 +e234 +e1234",
            |a| a.dual(),
        );
        assert_unary::<Mv4>(
            m,
            FULL_4D,
            "15 +14e1 -13e2 +11e3 -7e4 -12e12 +10e13 -6e14 -9e23 +5e24 -3e34 -8e123 +4e124 -2e134 +e234 +e1234",
            |a| a.undual(),
        );
        assert_unary::<Mv4>(m, FULL_4D, FULL_4D, |a| a.dual().undual());

        let m = "5,0,0";
        assert_unary::<Mv5>(
            m,
            FULL_5D,
            "15 +14e1 -13e2 +11e3 -7e4 +15e5 -12e12 +10e13 -6e14 +14e15 -9e23 +5e24 -13e25 -3e34 +11e35 -7e45 -8e123 +4e124 -12e125 -2e134 +10e135 -6e145 +e234 -9e235 +5e245 -3e345 +e1234 -8e1235 +4e1245 -2e1345 +e2345 +e12345",
            |a| a.dual(),
        );
        assert_unary::<Mv5>(
            m,
            FULL_5D,
            "15 +14e1 -13e2 +11e3 -7e4 +15e5 -12e12 +10e13 -6e14 +14e15 -9e23 +5e24 -13e25 -3e34 +11e35 -7e45 -8e123 +4e124 -12e125 -2e134 +10e135 -6e145 +e234 -9e235 +5e245 -3e345 +e1234 -8e1235 +4e1245 -2e1345 +e2345 +e12345",
            |a| a.undual(),
        );
        assert_unary::<Mv5>(m, FULL_5D, FULL_5D, |a| a.dual().undual());

        let m = "6,0,0";
        assert_unary::<Mv6>(
            m,
            FULL_6D,
            "15 +14e1 -13e2 +11e3 -7e4 +15e5 -15e6 -12e12 +10e13 -6e14 +14e15 -14e16 -9e23 +5e24 -13e25 +13e26 -3e34 +11e35 -11e36 -7e45 +7e46 -15e56 -8e123 +4e124 -12e125 +12e126 -2e134 +10e135 -10e136 -6e145 +6e146 -14e156 +e234 -9e235 +9e236 +5e245 -5e246 +13e256 -3e345 +3e346 -11e356 +7e456 +e1234 -8e1235 +8e1236 +4e1245 -4e1246 +12e1256 -2e1345 +2e1346 -10e1356 +6e1456 +e2345 -e2346 +9e2356 -5e2456 +3e3456 +e12345 -e12346 +8e12356 -4e12456 +2e13456 -e23456 -e123456",
            |a| a.dual(),
        );
        assert_unary::<Mv6>(
            m,
            FULL_6D,
            "-15 -14e1 +13e2 -11e3 +7e4 -15e5 +15e6 +12e12 -10e13 +6e14 -14e15 +14e16 +9e23 -5e24 +13e25 -13e26 +3e34 -11e35 +11e36 +7e45 -7e46 +15e56 +8e123 -4e124 +12e125 -12e126 +2e134 -10e135 +10e136 +6e145 -6e146 +14e156 -e234 +9e235 -9e236 -5e245 +5e246 -13e256 +3e345 -3e346 +11e356 -7e456 -e1234 +8e1235 -8e1236 -4e1245 +4e1246 -12e1256 +2e1345 -2e1346 +10e1356 -6e1456 -e2345 +e2346 -9e2356 +5e2456 -3e3456 -e12345 +e12346 -8e12356 +4e12456 -2e13456 +e23456 +e123456",
            |a| a.undual(),
        );
        assert_unary::<Mv6>(m, FULL_6D, FULL_6D, |a| a.dual().undual());
    }

    #[test]
    fn matches_golden_values_in_metrics_with_negative_vectors() {
        let m = "3,1,0";
        assert_unary::<Mv4>(
            m,
            FULL_4D,
            "15 +14e1 -13e2 +11e3 +7e4 -12e12 +10e13 +6e14 -9e23 -5e24 +3e34 -8e123 -4e124 +2e134 -e234 -e1234",
            |a| a.dual(),
        );
        assert_unary::<Mv4>(
            m,
            FULL_4D,
            " -15 -14e1 +13e2 -11e3 -7e4 +12e12 -10e13 -6e14 +9e23 +5e24 -3e34 +8e123 +4e124 -2e134 +e234 +e1234",
            |a| a.undual(),
        );
        assert_unary::<Mv4>(m, FULL_4D, FULL_4D, |a| a.undual().dual());

        let m = "2,2,0";
        assert_unary::<Mv4>(
            m,
            FULL_4D,
            "15 +14e1 -13e2 -11e3 +7e4 -12e12 -10e13 +6e14 +9e23 -5e24 -3e34 +8e123 -4e124 -2e134 +e234 +e1234",
            |a| a.dual(),
        );
        assert_unary::<Mv4>(
            m,
            FULL_4D,
            "15 +14e1 -13e2 -11e3 +7e4 -12e12 -10e13 +6e14 +9e23 -5e24 -3e34 +8e123 -4e124 -2e134 +e234 +e1234",
            |a| a.undual(),
        );
        assert_unary::<Mv4>(m, FULL_4D, FULL_4D, |a| a.dual().undual());

        let m = "2,3,0";
        assert_unary::<Mv5>(
            m,
            FULL_5D,
            "15 +14e1 -13e2 -11e3 +7e4 -15e5 -12e12 -10e13 +6e14 -14e15 +9e23 -5e24 +13e25 -3e34 +11e35 -7e45 +8e123 -4e124 +12e125 -2e134 +10e135 -6e145 +e234 -9e235 +5e245 +3e345 +e1234 -8e1235 +4e1245 +2e1345 -e2345 -e12345",
            |a| a.dual(),
        );
        assert_unary::<Mv5>(
            m,
            FULL_5D,
            " -15 -14e1 +13e2 +11e3 -7e4 +15e5 +12e12 +10e13 -6e14 +14e15 -9e23 +5e24 -13e25 +3e34 -11e35 +7e45 -8e123 +4e124 -12e125 +2e134 -10e135 +6e145 -e234 +9e235 -5e245 -3e345 -e1234 +8e1235 -4e1245 -2e1345 +e2345 +e12345",
            |a| a.undual(),
        );
        assert_unary::<Mv5>(m, FULL_5D, FULL_5D, |a| a.dual().undual());
    }

    #[test]
    fn matches_golden_values_in_degenerate_metrics() {
        let m = "2,0,1";
        assert_unary::<Mv3>(m, FULL_3D, "7 +3e3 +6e1 -5e2 -2e13 +e23 +4e12 +e123", |a| {
            a.dual()
        });
        assert_unary::<Mv3>(m, FULL_3D, FULL_3D, |a| a.dual().undual());

        let m = "1,0,2";
        assert_unary::<Mv3>(m, FULL_3D, "7 -5e2 +3e3 +6e1 +e23 +4e12 -2e13 +e123", |a| {
            a.dual()
        });
        assert_unary::<Mv3>(m, FULL_3D, FULL_3D, |a| a.dual().undual());

        let m = "3,0,1";
        assert_unary::<Mv4>(
            m,
            FULL_4D,
            "15 + 14e1 - 13e2 + 12e12 + 11e3 - 10e13 + 9e23 + 8e123 - 7e4 + 6e14 - 5e24 - 4e124 + 3e34 + 2e134 - 1e234 + 1e1234",
            |a| a.dual(),
        );
        assert_unary::<Mv4>(m, FULL_4D, FULL_4D, |a| a.dual().undual());

        let m = "4,0,1";
        assert_unary::<Mv5>(
            m,
            FULL_5D,
            "15 +15e5 +14e1 -13e2 +11e3 -7e4 -14e15 +13e25 -11e35 +7e45 +12e12 -10e13 +6e14 +9e23 -5e24 +3e34 +12e125 -10e135 +6e145 +9e235 -5e245 +3e345 +8e123 -4e124 +2e134 -e234 -8e1235 +4e1245 -2e1345 +e2345 +e1234 +e12345",
            |a| a.dual(),
        );
        assert_unary::<Mv5>(m, FULL_5D, FULL_5D, |a| a.dual().undual());

        // (input, expected dual), in the original order.
        let m = "3,0,2";
        let cases = [
            (
                PARTIAL_5D,
                "15 +14e1 -13e2 +7e45 -11e35 +12e12 +6e145 -5e245 +3e345 -10e135 +9e235 +4e1245 -2e1345 +e2345 -8e1235 +e12345",
            ),
            (
                "  + 3e13 + 1e4 + 2e124  +4e345",
                "+4e12 -2e35 -3e245 -e1235",
            ),
            (
                PARTIAL_5D,
                "15 +14e1 -13e2 +7e45 -11e35 +12e12 +6e145 -5e245 +3e345 -10e135 +9e235 +4e1245 -2e1345 +e2345 -8e1235 +e12345",
            ),
            (
                "  + 3e12 + 1e4 + 2e124  +4e345",
                "+4e12 -2e35 -e1235 + 3e345",
            ),
        ];
        for (input, expected) in cases {
            assert_unary::<Mv5>(m, input, expected, |a| a.dual());
            assert_unary::<Mv5>(m, input, input, |a| a.dual().undual());
        }
    }
}

/// Ported from `src/multivector/low_dimensional/core/linear_operations/geometric_product.rs`.
mod geometric_product {
    use super::*;

    #[test]
    fn matches_golden_values_in_euclidean_metrics() {
        let cases_3d = [
            ("1", "3e12", "3e12"),
            ("e1 + e2", "e2 + e3", "1 + e12 + e13 + e23"),
            (
                "1 + e1 + e3",
                "     e1  + 3e12",
                "1 + e1 + 3e2 +3e12 -e13 + 3e123",
            ),
            (
                "1 + 2e1 + 1e12 + 2e3",
                "  + 1e2 + 2e13 + 1e23",
                "-3e1 - e2 + 2e12 + 4e3 + 3e13 -3e23 + 2e123",
            ),
        ];
        for (lhs, rhs, expected) in cases_3d {
            assert_binary::<Mv3>("3,0,0", lhs, rhs, expected, |a, b| a * b);
        }

        assert_binary::<Mv4>(
            "4,0,0",
            " + 1 + e1 + 2e12 + e13 + e234 ",
            "  + 1e4 + 2e124 + 3e13 +4e34",
            "-3 -4e2 + 3e3 -3e4 +e13 +5e14 -5e23 + 2e24 + 4e34 + 7e124 + 5e134 + 2e234 + 8e1234",
            |a, b| a * b,
        );
    }

    #[test]
    fn matches_golden_values_in_pga3d() {
        let cases = [
            (
                "1 + 2e1 + 1e12 + 2e3",
                "  + 1e2 + 2e13 + 1e23",
                "-3e1 - e2 + 2e12 + 4e3 + 3e13 -3e23 + 2e123",
            ),
            (" + 2e3 + 2e14 + 1e124 ", "  + 1e24 + 2e134", "-4e14 -2e234"),
            (
                " + 1 + e1 + 2e12 + e13 + e234 ",
                "  + 3e13 + 1e4 + 2e124  +4e34",
                "-3 +3e3 + 3e13 -6e23 -3e4  + 5e14 +2e24 + 7e124 + 4e34  +5e134 + 2e234 + 8e1234",
            ),
            (
                "  + e1 + 2e2 - e234 + e1234 ",
                "  - 3e13 + 1e4 - 2e124  +4e34",
                "-3e3 +6e123 + 5e14 - 3e24 +3e124 + 4e134 + 8e234 ",
            ),
        ];
        for (lhs, rhs, expected) in cases {
            assert_binary::<Mv4>("3,0,1", lhs, rhs, expected, |a, b| a * b);
        }
    }

    #[test]
    fn matches_golden_values_in_hyperbolic_metric() {
        let cases = [
            (
                "1 + 2e1 + 1e12 + 2e3",
                "  + 1e2 + 2e13 + 1e23",
                "-3e1 - e2 + 2e12 + 4e3 + 3e13 -3e23 + 2e123",
            ),
            (
                " + e1 + 2e2 - e234 + e1234  ",
                "  - 3e13 + 1e4 - 2e124  +4e34",
                "-4e2 -5e3 + 4e12 + 2e13 + 5e14 + e23 -3e24 + 5e123 + 3e124 + 4e134 + 8e234",
            ),
            (
                " + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e34 + 13e134 + 14e234 + 15e1234",
                "  + 3e13 + 1e4 + 2e124  +4e34",
                "47 +51e1 + 49e2 + 53e3 - 28e4 + 51e12 + 54e13 - 19e14 - 9e23 + 73e24 + 21e34 - e123 + 75e124 + 21e134 -9e234 -3e1234",
            ),
        ];
        for (lhs, rhs, expected) in cases {
            assert_binary::<Mv4>("3,1,0", lhs, rhs, expected, |a, b| a * b);
        }
    }

    #[test]
    fn matches_golden_values_in_mixed_5d_metrics() {
        let cases = [
            (
                "3,2,0",
                "-49 -53e1 -63e2 +3e3 -5e4 -69e12 + 3e13 - 3e14 - 9e23 + 4e24 + 17e34 + 14e35 -23e45 -6e123 + 5e124 +17e134 + 21e135 - 16e145 -17e234 + 80e235 + 69e245 + 4e345 - 15e1234 +83e1235 + 70e1245 + 4e1345 + 8e2345 + 12e12345",
            ),
            (
                "3,1,1",
                "-1 -e1 -7e2 +3e3 - 5e4 +14e35 -23e45 -9e12 +3e13 - 3e14 -9e23 +4e24 +17e34 +21e135 - 16e145 +80e235 +69e245 + 4e345 -6e123 + 5e124 + 17e134 - 17e234 +83e1235 +70e1245 +4e1345 + 8e2345 -15e1234 +12e12345",
            ),
        ];
        for (metric, expected) in cases {
            assert_binary::<Mv5>(metric, A_5D, B_5D, expected, |a, b| a * b);
        }
    }

    #[test]
    fn matches_golden_values_in_6d_metric() {
        let cases = [
            (
                A_6D,
                B_6D,
                "115+115e1-27e2+1053e3+1240e4+189e5+188e6-21e12+983e13+1159e14+202e15+168e16-919e23-1075e24-215e25-156e26+17e34-351e35-301e45+23e46-846e123-991e124-228e125-128e126+17e134-265e135+3e136-228e145-17e234+143e235-36e236+135e245+303e345-1160e356-924e456-15e1234+230e1235+124e1245-14e1246-560e1256+442e1345+24e1346-152e1356-72e1456-393e2345+44e2356+60e2456-40e3456-360e12345+179e12356+118e12456+76e13456-28e23456-700e123456",
            ),
            (
                "1 + 2e1 + 1e12 + 2e3",
                "  + 1e2 + 2e13 + 1e23",
                "-3e1 - e2 + 2e12 + 4e3 + 3e13 -3e23 + 2e123",
            ),
        ];
        for (lhs, rhs, expected) in cases {
            assert_binary::<Mv6>("3,3,0", lhs, rhs, expected, |a, b| a * b);
        }
    }
}

/// Ported from `src/multivector/low_dimensional/core/linear_operations/inner_product.rs`.
mod inner_product {
    use super::*;

    #[test]
    fn matches_golden_values_in_euclidean_metrics() {
        let cases_3d = [
            ("1", "3e12", "3e12"),
            ("e1 + e2", "e2 + e3", "1"),
            ("1 + e1 + e3", "     e1  + 3e12", "1 + e1 + 3e2 + 3e12"),
            ("e12", "e2", "e1"),
            (
                "1 + 2e1 + 1e12 + 2e3",
                "  + 1e2 + 2e13 + 1e23",
                "-3e1 -e2 +4e3 +2e13 + e23",
            ),
        ];
        for (lhs, rhs, expected) in cases_3d {
            assert_binary::<Mv3>("3,0,0", lhs, rhs, expected, |a, b| a | b);
        }

        assert_binary::<Mv4>(
            "4,0,0",
            " + 1 + e1 + 2e12 + e13 + e234 ",
            "  + 1e4 + 2e124 + 3e13 +4e34",
            "-3 -4e2 +3e3 -3e4 +3e13 +e23 +2e24 +4e34 +2e124",
            |a, b| a | b,
        );
    }

    #[test]
    fn matches_golden_values_in_pga3d() {
        let cases = [
            (
                "  + e1 + 2e2 - e234 + e1234 ",
                "  - 3e13 + 1e4 - 2e124  +4e34",
                "-3e3 +4e14 -5e24",
            ),
            (
                " + 1 + e1 + 2e12 + e13 + e234 ",
                "  + 3e13 + 1e4 + 2e124  +4e34",
                "-3 -3e4 +3e3 +2e24 +4e34 +3e13 +2e124",
            ),
        ];
        for (lhs, rhs, expected) in cases {
            assert_binary::<Mv4>("3,0,1", lhs, rhs, expected, |a, b| a | b);
        }
    }

    #[test]
    fn matches_golden_values_in_hyperbolic_metric() {
        let cases = [
            (
                "1 + 2e1 + 1e12 + 2e3",
                "  + 1e2 + 2e13 + 1e23",
                "-3e1 -e2 +4e3 +2e13 +e23",
            ),
            (
                " + e1 + 2e2 - e234 + e1234  ",
                "  - 3e13 + 1e4 - 2e124  +4e34",
                "-4e2 -5e3 +4e12 +4e14 +e23 -5e24 -e123",
            ),
            (
                " + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e34 + 13e134 + 14e234 + 15e1234",
                "  + 3e13 + 1e4 + 2e124  +4e34",
                "47 +51e1 +49e2 +53e3 -28e4 +33e12 -10e13 -4e14 -14e23 +47e24 +4e34 -15e123 +2e124",
            ),
        ];
        for (lhs, rhs, expected) in cases {
            assert_binary::<Mv4>("3,1,0", lhs, rhs, expected, |a, b| a | b);
        }
    }

    #[test]
    fn matches_golden_values_in_mixed_5d_metrics() {
        let cases = [
            (
                "3,2,0",
                B_5D,
                " -49 -53e1 -63e2 +3e3 -5e4 -87e12 +3e13 -4e14 +2e24 +14e35 -23e45 +2e124 +13e135 +14e235 +45e245 +4e345 +15e1235",
            ),
            (
                "3,1,1",
                "  + 3e13 + 1e4 + 2e124  + 4e345",
                " -1 -e1 -7e2 +3e3 -5e4 +14e35 -23e45 -27e12 +3e13 -4e14 +2e24 +13e135 +14e235 +45e245 +4e345 +2e124 +15e1235",
            ),
        ];
        for (metric, rhs, expected) in cases {
            assert_binary::<Mv5>(metric, A_5D, rhs, expected, |a, b| a | b);
        }
    }

    #[test]
    fn matches_golden_values_in_6d_metric() {
        let cases = [
            (
                A_6D,
                B_6D,
                "115+115e1 -27e2 +1053e3 +1240e4 +189e5 +188e6 -87e12 +3e13 -4e14 +202e15 +168e16 +2e24 -215e25 -156e26 -351e35 -301e45 +2e124 -228e125 -144e126 -185e135 -156e145 +33e235 +83e245 -120e345 -132e356 -84e456 +194e1235 +70e1245 +4e1345 -120e1356 -72e1456 +8e2345 +108e2356 +60e2456 -36e3456 +12e12345 +96e12356 +48e12456 -24e13456 +12e23456 +12e123456",
            ),
            (
                "1 + 2e1 + 1e12 + 2e3",
                "  + 1e2 + 2e13 + 1e23",
                "-3e1 -e2 +4e3 +2e13 +e23",
            ),
        ];
        for (lhs, rhs, expected) in cases {
            assert_binary::<Mv6>("3,3,0", lhs, rhs, expected, |a, b| a | b);
        }
    }
}

/// Ported from `src/multivector/low_dimensional/core/linear_operations/left_contraction.rs`.
mod left_contraction {
    use super::*;

    #[test]
    fn matches_golden_values_in_euclidean_metrics() {
        let cases_3d = [
            ("1", "3e12", "3e12"),
            ("e1 + e2", "e2 + e3", "1"),
            ("1 + e1 + e3", "     e1  + 3e12", "1 + e1 + 3e2 +3e12"),
            (
                "1 + 2e1 + 1e12 + 2e3",
                "  + 1e2 + 2e13 + 1e23",
                "-4e1 - e2 + 4e3 + 2e13 +1e23",
            ),
        ];
        for (lhs, rhs, expected) in cases_3d {
            assert_binary::<Mv3>("3,0,0", lhs, rhs, expected, |a, b| a << b);
        }

        assert_binary::<Mv4>(
            "4,0,0",
            " + 1 + e1 + 2e12 + e13 + e234 ",
            "  + 1e4 + 2e124 + 3e13 +4e34",
            "-3 + 3e3 -3e4 + 3e13 + 2e24 + 4e34 + 2e124",
            |a, b| a << b,
        );
    }

    #[test]
    fn matches_golden_values_in_pga3d() {
        let cases = [
            (
                "  + e1 + 2e2 - e234 + e1234 ",
                "  - 3e13 + 1e4 - 2e124  +4e34",
                "-3e3 +4e14 -2e24",
            ),
            (
                " + 1 + e1 + 2e12 + e13 + e234 ",
                "  + 3e13 + 1e4 + 2e124  +4e34",
                "-3 -3e4 + 3e3 +2e24 + 4e34 + 3e13 + 2e124",
            ),
        ];
        for (lhs, rhs, expected) in cases {
            assert_binary::<Mv4>("3,0,1", lhs, rhs, expected, |a, b| a << b);
        }
    }

    #[test]
    fn matches_golden_values_in_hyperbolic_metric() {
        let cases = [
            (
                "1 + 2e1 + 1e12 + 2e3",
                "  + 1e2 + 2e13 + 1e23",
                "-4e1 -e2 +4e3 +2e13 +e23",
            ),
            (
                " + e1 + 2e2 - e234 + e1234  ",
                "  - 3e13 + 1e4 - 2e124  +4e34",
                "-3e3 + 4e14 -2e24",
            ),
            (
                " + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e34 + 13e134 + 14e234 + 15e1234",
                "  + 3e13 + 1e4 + 2e124  +4e34",
                "47 +8e1 -18e2 +35e3 + 11e4 -16e12 + 3e13 -4e14 +2e24 +4e34 + 2e124",
            ),
        ];
        for (lhs, rhs, expected) in cases {
            assert_binary::<Mv4>("3,1,0", lhs, rhs, expected, |a, b| a << b);
        }
    }

    #[test]
    fn matches_golden_values_in_mixed_5d_metrics() {
        let cases = [
            (
                "3,2,0",
                "-49 +8e1 -18e2 +3e3 -5e4 -16e12 +3e13 -4e14 +2e24 +32e35 +16e45 +2e124 + 4e345",
            ),
            (
                "3,1,1",
                "-1 +8e1 -18e2 + 3e3 -5e4 +32e35 +16e45 -16e12 + 3e13 -4e14 +2e24 +4e345 +2e124",
            ),
        ];
        for (metric, expected) in cases {
            assert_binary::<Mv5>(metric, A_5D, B_5D, expected, |a, b| a << b);
        }
    }

    #[test]
    fn matches_golden_values_in_6d_metric() {
        let cases = [
            (
                A_6D,
                B_6D,
                "115+56e1+78e2+3e3-5e4+189e5+180e6-16e12+3e13-4e14+202e15+168e16+2e24-215e25-156e26-333e35-262e45+2e124-228e125-144e126-198e135-156e145+19e235+38e245-120e345-132e356-84e456+179e1235+70e1245+4e1345-120e1356-72e1456+8e2345+108e2356+60e2456-36e3456+12e12345+96e12356+48e12456-24e13456+12e23456+12e123456",
            ),
            (
                "1 + 2e1 + 1e12 + 2e3",
                "  + 1e2 + 2e13 + 1e23",
                "-4e1 -e2 +4e3 +2e13 +e23",
            ),
        ];
        for (lhs, rhs, expected) in cases {
            assert_binary::<Mv6>("3,3,0", lhs, rhs, expected, |a, b| a << b);
        }
    }
}

/// Ported from `src/multivector/low_dimensional/core/linear_operations/outer_product.rs`.
/// The outer product does not depend on the metric, so every case uses the Euclidean one.
mod outer_product {
    use super::*;

    #[test]
    fn matches_golden_values_in_3d_to_6d() {
        let cases_3d = [
            ("1", "3e12", "3e12"),
            ("e1 + e2", "e2 + e3", "e12 + e13 + e23"),
            ("1 + e1 + e3", "     e1  + 3e12", "e1 + 3e12 - e13 + 3e123"),
            ("e12", "e2", ""),
            (
                "1 + 2e1 + 1e12 + 2e3",
                "  + 1e2 + 2e13 + 1e23",
                "e2 +2e12 +2e13 -e23 +2e123",
            ),
        ];
        for (lhs, rhs, expected) in cases_3d {
            assert_binary::<Mv3>("3", lhs, rhs, expected, |a, b| a ^ b);
        }

        assert_binary::<Mv4>(
            "4",
            " + 1 + e1 + 2e12 + e13 + e234 ",
            "  + 1e4 + 2e124 + 3e13 +4e34",
            "e4 +3e13 +e14 +4e34 +4e124 +5e134 +8e1234",
            |a, b| a ^ b,
        );

        let cases_5d = [
            (
                A_5D,
                B_5D,
                "e4 +3e13 +e14 +2e24 +4e34 -6e123 +5e124 +29e134 +6e234 +4e345 -15e1234 +4e1345 +8e2345 +12e12345",
            ),
            (
                "15 +14e1 -13e2 +7e45 -11e35 +12e12 +6e145 -5e245 +3e345 -10e135 +9e235 +4e1245 -2e1345 +e2345 -8e1235 +e12345",
                "-2e35 +4e12 +3e345 -e1235",
                "-30e35 +60e12 +45e345 -28e135 +26e235 +28e1245 +42e1345 -39e2345 -83e1235 +48e12345",
            ),
        ];
        for (lhs, rhs, expected) in cases_5d {
            assert_binary::<Mv5>("5", lhs, rhs, expected, |a, b| a ^ b);
        }

        assert_binary::<Mv6>(
            "6",
            A_6D,
            B_6D,
            "e4 +3e13 +e14 +2e24 +4e34 -e46 -6e123 +5e124 +29e134 +3e136 +6e234 +4e345 -15e1234 +83e1235 +70e1245 -14e1246 +8e1345 +24e1346 +16e2345 -4e3456 -372e12345 +83e12356 +70e12456 +4e13456 +8e23456 -700e123456",
            |a, b| a ^ b,
        );
    }
}
