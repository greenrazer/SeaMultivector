//! Golden values from the 0.1 engine's unit tests: the regressive product, right contraction,
//! sandwich product, scalar product, addition, exponential and normalization, in metrics from
//! 3D to 6D. Each module is named after the 0.1 file its cases came from.

use std::fmt::Debug;

use seagrass_math::prelude::*;

/// Parses `body` under the metric prefix `metric`, e.g. `parse("3,0,1", "e1 + 2e12")`.
fn parse<M>(metric: &str, body: &str) -> M
where
    M: FromStr<Err = ParseError>,
{
    let s = format!("{metric}>{body}");
    s.parse()
        .unwrap_or_else(|e| panic!("cannot parse {s:?}: {e:?}"))
}

/// Asserts `op(lhs, rhs) == expected` for every `(lhs, rhs, expected)` row, each parsed under
/// `metric`. The comparison is structural: the stored blades must match, not only the values.
#[track_caller]
fn assert_cases<M, F>(metric: &str, op: F, cases: &[(&str, &str, &str)])
where
    M: FromStr<Err = ParseError> + PartialEq + Debug,
    F: Fn(&M, &M) -> M,
{
    for &(lhs, rhs, expected) in cases {
        let actual = op(&parse(metric, lhs), &parse(metric, rhs));
        assert_eq!(
            actual,
            parse(metric, expected),
            "metric {metric}, lhs {lhs:?}, rhs {rhs:?}"
        );
    }
}

/// Asserts that `actual` stores exactly the blades of `expected` and that every coefficient
/// differs from the expected one by less than `tol`. A macro rather than a function because
/// the blade type yielded by `iter()` is not nameable through the prelude.
macro_rules! assert_approx_eq {
    ($actual:expr, $expected:expr, $tol:expr) => {{
        let (actual, expected, tol) = (&$actual, &$expected, $tol);
        let actual_blades: Vec<(u8, f32)> =
            actual.iter().map(|(index, c)| (index.bits(), *c)).collect();
        let expected_blades: Vec<(u8, f32)> =
            expected.iter().map(|(index, c)| (index.bits(), *c)).collect();
        let indices = |blades: &[(u8, f32)]| blades.iter().map(|b| b.0).collect::<Vec<_>>();
        assert_eq!(
            indices(&actual_blades),
            indices(&expected_blades),
            "stored blades differ\n  actual:   {actual:?}\n  expected: {expected:?}"
        );
        for (&(index, a), &(_, e)) in actual_blades.iter().zip(&expected_blades) {
            assert!(
                (a - e).abs() < tol,
                "blade {index:#b}: |{a} - {e}| >= {tol}\n  actual:   {actual:?}\n  expected: {expected:?}"
            );
        }
    }};
}

mod regressive_product {
    use super::*;

    #[test]
    fn matches_golden_values_in_euclidean_space() {
        assert_cases::<Multivector<f32, Dyn3>, _>(
            "3,0,0",
            |a, b| a & b,
            &[
                ("1", "3e12", ""),
                ("e1 + e2", "e2 + e3", ""),
                ("1 + e1 + e3", "     e1  + 3e12", "3"),
                (
                    "1 + 2e1 + 1e12 + 2e3",
                    "  + 1e2 + 2e13 + 1e23",
                    "2 -2e1 -e2",
                ),
            ],
        );
        assert_cases::<Multivector<f32, Dyn4>, _>(
            "4,0,0",
            |a, b| a & b,
            &[(
                " + 1 + e1 + 2e12 + e13 + e234 ",
                "  + 1e4 + 2e124 + 3e13 +4e34",
                "8 -2e1 +3e3 +2e24",
            )],
        );
    }

    #[test]
    fn matches_golden_values_with_null_vectors() {
        assert_cases::<Multivector<f32, Dyn4>, _>(
            "3,0,1",
            |a, b| a & b,
            &[
                ("1 + 2e1 + 1e12 + 2e3", "  + 1e2 + 2e13 + 1e23", ""),
                (" + 2e3 + 2e14 + 1e124 ", "  + 1e24 + 2e134", "2e14"),
                (
                    " + 1 + e1 + 2e12 + e13 + e234 ",
                    "  + 3e13 + 1e4 + 2e124  +4e34",
                    "+8 -2e1 +3e3 -2e24",
                ),
                (
                    "  + e1 + 2e2 - e234 + e1234 ",
                    "  - 3e13 + 1e4 - 2e124  +4e34",
                    "+e4 +3e3 -2e24 +4e34 -3e13 -2e124",
                ),
            ],
        );
        assert_cases::<Multivector<f32, Dyn5>, _>(
            "3,0,2",
            |a, b| a & b,
            &[(
                " + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345",
                "  + 3e12 + 1e4 + 2e124  +4e345",
                "48 +83e4 -39e1 -42e2 +28e3 +26e14 +28e24 +45e12 +60e345 +30e124",
            )],
        );
    }

    #[test]
    fn matches_golden_values_with_negative_vectors() {
        assert_cases::<Multivector<f32, Dyn4>, _>(
            "3,1,0",
            |a, b| a & b,
            &[
                ("1 + 2e1 + 1e12 + 2e3", "  + 1e2 + 2e13 + 1e23", ""),
                (
                    " + e1 + 2e2 - e234 + e1234  ",
                    "  - 3e13 + 1e4 - 2e124  +4e34",
                    "3e3 +e4 -3e13 +2e24 +4e34 -2e124",
                ),
                (
                    " + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e34 + 13e134 + 14e234 + 15e1234",
                    "  + 3e13 + 1e4 + 2e124  +4e34",
                    " -33 -43e1 -12e2 +70e3 +83e4 -14e12 +45e13 +26e14 +28e24 +60e34 +30e124",
                ),
            ],
        );
    }

    #[test]
    fn matches_golden_values_with_mixed_signatures() {
        assert_cases::<Multivector<f32, Dyn5>, _>(
            "3,2,0",
            |a, b| a & b,
            &[(
                " + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345",
                "  + 3e13 + 1e4 + 2e124  +4e345",
                "12 +70e3 +83e4 +45e13 +26e14 +28e24 +30e124 +60e345",
            )],
        );
        assert_cases::<Multivector<f32, Dyn5>, _>(
            "3,1,1",
            |a, b| a & b,
            &[(
                " + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345",
                "  + 3e13 + 1e4 + 2e124  +4e345",
                "12 -14e3 +83e4 +45e13 +26e14 +28e24 +60e345 +30e124",
            )],
        );
    }

    #[test]
    fn matches_golden_values_in_six_dimensions() {
        assert_cases::<Multivector<f32, Dyn6>, _>(
            "3,3,0",
            |a, b| a & b,
            &[
                (
                    " + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345 +e6 + 12e126 + 8e46",
                    "  + 3e13 + 1e4 + 2e124  +4e345 +83e1235 + 70e1245 + 4e1345 + 8e2345 + 12e123456",
                    " -604 -36e1 -72e2 +48e3 +96e4 +12e6 +36e12 +60e13 +108e14 +72e23 +120e24 +96e46 +84e123 +132e124 +144e126 +144e345 +156e1345 +168e2345 +180e12345",
                ),
                ("1 + 2e1 + 1e12 + 2e3", "  + 1e2 + 2e13 + 1e23", ""),
            ],
        );
    }
}

mod right_contraction {
    use super::*;

    #[test]
    fn matches_golden_values_in_euclidean_space() {
        assert_cases::<Multivector<f32, Dyn3>, _>(
            "3,0,0",
            |a, b| a >> b,
            &[
                ("1", "3e12", ""),
                ("e1 + e2", "e2 + e3", "1"),
                ("1 + e1 + e3", "     e1  + 3e12", "1"),
                ("e12", "e2", "e1"),
                ("1 + 2e1 + 1e12 + 2e3", "  + 1e2 + 2e13 + 1e23", "e1"),
            ],
        );
        assert_cases::<Multivector<f32, Dyn4>, _>(
            "4,0,0",
            |a, b| a >> b,
            &[(
                " + 1 + e1 + 2e12 + e13 + e234 ",
                "  + 1e4 + 2e124 + 3e13 +4e34",
                "-3 -4e2 + e23",
            )],
        );
    }

    #[test]
    fn matches_golden_values_with_null_vectors() {
        assert_cases::<Multivector<f32, Dyn4>, _>(
            "3,0,1",
            |a, b| a >> b,
            &[
                (
                    "  + e1 + 2e2 - e234 + e1234 ",
                    "  - 3e13 + 1e4 - 2e124  +4e34",
                    "-3e24",
                ),
                (
                    " + 1 + e1 + 2e12 + e13 + e234 ",
                    "  + 3e13 + 1e4 + 2e124  +4e34",
                    "-3",
                ),
            ],
        );
    }

    #[test]
    fn matches_golden_values_with_negative_vectors() {
        assert_cases::<Multivector<f32, Dyn4>, _>(
            "3,1,0",
            |a, b| a >> b,
            &[
                ("1 + 2e1 + 1e12 + 2e3", "  + 1e2 + 2e13 + 1e23", "e1"),
                (
                    " + e1 + 2e2 - e234 + e1234  ",
                    "  - 3e13 + 1e4 - 2e124  +4e34",
                    "-4e2 - 2e3 +4e12 + e23 - 3e24 - e123",
                ),
                (
                    " + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e34 + 13e134 + 14e234 + 15e1234",
                    "  + 3e13 + 1e4 + 2e124  +4e34",
                    "47+43e1+67e2+18e3-39e4+49e12-13e13-14e23+45e24-15e123",
                ),
            ],
        );
    }

    #[test]
    fn matches_golden_values_with_mixed_signatures() {
        assert_cases::<Multivector<f32, Dyn5>, _>(
            "3,2,0",
            |a, b| a >> b,
            &[(
                " + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345",
                "  + 3e13 + 1e4 + 2e124  +4e345",
                "-49-61e1-45e2-71e12-18e35-39e45+13e135+14e235+45e245+15e1235",
            )],
        );
        assert_cases::<Multivector<f32, Dyn5>, _>(
            "3,1,1",
            |a, b| a >> b,
            &[(
                " + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345",
                "  + 3e13 + 1e4 + 2e124  +4e345",
                "-1 -9e1 +11e2 -18e35 -39e45 -11e12 +13e135 +14e235 +45e245 +15e1235",
            )],
        );
    }

    #[test]
    fn matches_golden_values_in_six_dimensions() {
        assert_cases::<Multivector<f32, Dyn6>, _>(
            "3,3,0",
            |a, b| a >> b,
            &[
                (
                    " + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345 +e6 + 12e126 + 8e46",
                    "  + 3e13 + 1e4 + 2e124  +4e345 +83e1235 + 70e1245 + 4e1345 + 8e2345 + 12e123456",
                    "115+59e1-105e2+1050e3+1245e4+8e6-71e12-18e35-39e45+13e135+14e235+45e245+15e1235",
                ),
                ("1 + 2e1 + 1e12 + 2e3", "  + 1e2 + 2e13 + 1e23", "e1"),
            ],
        );
    }
}

mod sandwich_product {
    use super::*;

    #[test]
    fn matches_golden_values_in_euclidean_space() {
        assert_cases::<Multivector<f32, Dyn3>, _>(
            "3,0,0",
            |a, b| a.sandwich(b),
            &[
                ("1", "3e12", "3e12"),
                ("e1 + e2", "e2 + e3", "2e1 - 2e3"),
                (
                    "1 + e1 + e3",
                    "     e1  + 3e12",
                    "2 + e1 + 2e3 + 3e12 + 6e23 + 6e123",
                ),
                (
                    "1 + 2e1 + 1e12 + 2e3",
                    "  + 1e2 + 2e13 + 1e23",
                    "4 +2e1  -8e2 +8e12 -14e13 -4e23 -4e123",
                ),
            ],
        );
        assert_cases::<Multivector<f32, Dyn4>, _>(
            "4,0,0",
            |a, b| a.sandwich(b),
            &[(
                " + 1 + e1 + 2e12 + e13 + e234 ",
                "  + 1e4 + 2e124 + 3e13 +4e34",
                "6e4 + 20e12 - 16e13 + 16e14 - 4e23 - 6e24 + 24e34 - 20e123 + 14e124 + 16e134 - 6e234",
            )],
        );
    }

    #[test]
    fn matches_golden_values_with_null_vectors() {
        assert_cases::<Multivector<f32, Dyn4>, _>(
            "3,0,1",
            |a, b| a.sandwich(b),
            &[
                // Originally run on the 3D engine type under this 4D metric. The runtime-metric
                // types reject a metric whose dimension differs from theirs, so it runs in 4D.
                (
                    "1 + 2e1 + 1e12 + 2e3",
                    "  + 1e2 + 2e13 + 1e23",
                    "4 +2e1  -8e2 +8e12 -14e13 -4e23 -4e123",
                ),
                (
                    " + 2e3 + 2e14 + 1e124 ",
                    "  + 1e24 + 2e134",
                    "+4e24 + 8e134",
                ),
                (
                    " + 1 + e1 + 2e12 + e13 + e234 ",
                    "  + 3e13 + 1e4 + 2e124  +4e34",
                    "5e4 + 16e14 - 6e24 + 20e34 + 12e12 - 9e13 - 12e23 + 16e124 + 16e134 - 6e234 - 12e123",
                ),
                (
                    "  + e1 + 2e2 - e234 + e1234 ",
                    "  - 3e13 + 1e4 - 2e124  +4e34",
                    "- 5e4 - 12e14 + 6e24 + 20e34 - 9e13 + 12e23 - 16e124 ",
                ),
            ],
        );
    }

    #[test]
    fn matches_golden_values_with_negative_vectors() {
        assert_cases::<Multivector<f32, Dyn4>, _>(
            "3,1,0",
            |a, b| a.sandwich(b),
            &[
                (
                    "1 + 2e1 + 1e12 + 2e3",
                    "  + 1e2 + 2e13 + 1e23",
                    "4 + 2e1 - 8e2 + 8e12 - 14e13 - 4e23 - 4e123",
                ),
                (
                    " + e1 + 2e2 - e234 + e1234  ",
                    "  - 3e13 + 1e4 - 2e124  +4e34",
                    " - 4e3 - 5e4 + 8e12 - 5e13 - 12e14 + 20e23 + 10e24 + 12e34 - 8e123 - 20e124 + 8e134",
                ),
                (
                    " + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e34 + 13e134 + 14e234 + 15e1234",
                    "  + 3e13 + 1e4 + 2e124  +4e34",
                    " - 744 - 738e1 + 100e2 - 40e3 + e4 - 164e12 - 109e13 + 320e14 - 2454e23 + 2078e24 - 2610e34 - 2452e123 + 2102e124 - 2608e134 + 146e234 + 14e1234",
                ),
            ],
        );
    }
}

mod scalar_product {
    use super::*;

    /// Asserts `lhs.scalar_product(&rhs) == expected` exactly for every row, each operand
    /// parsed under `metric`.
    #[track_caller]
    fn assert_scalar_cases<M>(metric: &str, cases: &[(&str, &str, f32)])
    where
        M: FromStr<Err = ParseError> + ScalarProduct<M, Output = f32>,
    {
        for &(lhs, rhs, expected) in cases {
            let lhs_mv: M = parse(metric, lhs);
            assert_eq!(
                lhs_mv.scalar_product(&parse::<M>(metric, rhs)),
                expected,
                "metric {metric}, lhs {lhs:?}, rhs {rhs:?}"
            );
        }
    }

    #[test]
    fn matches_golden_values_in_euclidean_space() {
        assert_scalar_cases::<Multivector<f32, Dyn3>>(
            "3,0,0",
            &[
                ("1", "3e12", 0.0),
                ("e1 + e2", "e2 + e3", 1.0),
                ("1 + e1 + e3", "     e1  + 3e12", 1.0),
                ("1 + 2e1 + 1e12 + 2e3", "  + 1e2 + 2e13 + 1e23", 0.0),
            ],
        );
        assert_scalar_cases::<Multivector<f32, Dyn4>>(
            "4,0,0",
            &[(
                " + 1 + e1 + 2e12 + e13 + e234 ",
                " + 3e13 + 1e4 + 2e124 +4e34",
                -3.0,
            )],
        );
    }

    #[test]
    fn matches_golden_values_with_null_vectors() {
        assert_scalar_cases::<Multivector<f32, Dyn4>>(
            "3,0,1",
            &[
                ("1 + 2e1 + 1e12 + 2e3", "  + 1e2 + 2e13 + 1e23", 0.0),
                (" + 2e3 + 2e14 + 1e124 ", "  + 1e24 + 2e134", 0.0),
                (
                    " + 1 + e1 + 2e12 + e13 + e234 ",
                    "  + 3e13 + 1e4 + 2e124  +4e34",
                    -3.0,
                ),
                (
                    "  + e1 + 2e2 - e234 + e1234 ",
                    "  - 3e13 + 1e4 - 2e124  +4e34",
                    0.0,
                ),
            ],
        );
    }

    #[test]
    fn matches_golden_values_with_negative_vectors() {
        assert_scalar_cases::<Multivector<f32, Dyn4>>(
            "3,1,0",
            &[
                ("1 + 2e1 + 1e12 + 2e3", "  + 1e2 + 2e13 + 1e23", 0.0),
                (
                    " + e1 + 2e2 - e234 + e1234  ",
                    "  - 3e13 + 1e4 - 2e124  +4e34",
                    0.0,
                ),
                (
                    " + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e34 + 13e134 + 14e234 + 15e1234",
                    "  + 3e13 + 1e4 + 2e124  +4e34",
                    47.0,
                ),
            ],
        );
    }

    #[test]
    fn matches_golden_values_with_mixed_signatures() {
        assert_scalar_cases::<Multivector<f32, Dyn5>>(
            "3,2,0",
            &[(
                " + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345",
                "  + 3e13 + 1e4 + 2e124  +4e345",
                -49.0,
            )],
        );
        assert_scalar_cases::<Multivector<f32, Dyn5>>(
            "3,1,1",
            &[(
                " + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345",
                "  + 3e13 + 1e4 + 2e124  +4e345",
                -1.0,
            )],
        );
    }

    #[test]
    fn matches_golden_values_in_six_dimensions() {
        assert_scalar_cases::<Multivector<f32, Dyn6>>(
            "3,3,0",
            &[
                (
                    " + 1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e345 + 13e1345 + 14e2345 + 15e12345 +e6 + 12e126 + 8e46",
                    "  + 3e13 + 1e4 + 2e124  +4e345 +83e1235 + 70e1245 + 4e1345 + 8e2345 + 12e123456",
                    115.0,
                ),
                ("1 + 2e1 + 1e12 + 2e3", "  + 1e2 + 2e13 + 1e23", 0.0),
            ],
        );
    }
}

mod addition {
    use super::*;

    #[test]
    fn sums_coefficients_of_matching_blades() {
        let mut a = parse::<Multivector<f32, Dyn3>>(
            "3",
            "-1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123",
        );
        let b = parse::<Multivector<f32, Dyn3>>(
            "3",
            "-1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123",
        );
        let sum = parse::<Multivector<f32, Dyn3>>(
            "3",
            "-2 + 2e1 + 4e2 + 6e12 + 8e3 + 10e13 + 12e23 + 14e123",
        );

        assert_eq!(&a + &b, sum);

        a += &b;
        assert_eq!(a, sum);
    }

    #[test]
    fn inserts_added_basis_blades() {
        let mut a =
            parse::<Multivector<f32, Dyn3>>("3", "-1 + e1 + 2e2 + 3e12 + 5e13 + 6e23 + 7e123");
        a += &parse::<Multivector<f32, Dyn3>>("3", "4e3");
        assert_eq!(
            a,
            parse("3", "-1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123")
        );

        let mut a = parse::<Multivector<f32, Dyn3>>("3", "-1 + e1 + 2e2 + 3e12 + 5e13");
        a += &parse::<Multivector<f32, Dyn3>>("3", "4e3");
        a += &parse::<Multivector<f32, Dyn3>>("3", "7e123");
        a += &parse::<Multivector<f32, Dyn3>>("3", "6e23");
        assert_eq!(
            a,
            parse("3", "-1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123")
        );
    }
}

mod exponentiation {
    use super::*;

    #[test]
    fn matches_golden_values_when_square_is_scalar() {
        for (input, expected) in [
            (
                "3e12 + 2e13 + 3e23",
                "-0.021971453 -0.6394477487e12 -0.4262984991e13 -0.6394477487e23",
            ),
            ("e1 + e2", "2.1781836 + 1.3682989e1 + 1.3682989e2"),
        ] {
            let actual = parse::<Multivector<f32, Dyn3>>("3,0,0", input).exp();
            assert_approx_eq!(
                actual,
                parse::<Multivector<f32, Dyn3>>("3,0,0", expected),
                1e-6
            );
        }

        for (input, expected) in [
            ("3e12", "-0.9899924994 + 0.1411200017e12"),
            ("3e34", "1 + 3e34"),
            ("e1 + e4", "1.5430807 + 1.1752012e1 + 1.1752012e4"),
            (
                "e12 + 2e13 + 3e23",
                "-0.8252990842 -0.1509213299e12 -0.3018426597e13 -0.4527639747e23",
            ),
            ("4e14 + 5e24 + 6e34", "1 +4e14 +5e24 +6e34"),
        ] {
            let actual = parse::<Multivector<f32, Dyn4>>("3,0,1", input).exp();
            assert_approx_eq!(
                actual,
                parse::<Multivector<f32, Dyn4>>("3,0,1", expected),
                1e-6
            );
        }
    }

    #[test]
    fn matches_golden_values_when_square_is_not_scalar() {
        for (input, expected, tol) in [(
            "1 + e1 + e2 + e12 + e3 + e13 + e23 + e123",
            " -1.0374662 -2.0374662e1 +3.3594248e2 -2.0374662e3 +3.3594248e12 -2.0374662e13 +3.3594248e23 +3.3594248e123",
            1e-4,
        )] {
            let actual = parse::<Multivector<f32, Dyn3>>("3,0,0", input).exp();
            assert_approx_eq!(
                actual,
                parse::<Multivector<f32, Dyn3>>("3,0,0", expected),
                tol
            );
        }

        for (input, expected, tol) in [
            (
                "e12 + 2e13 + 3e23 + 4e14 + 5e24 + 6e34",
                "-0.82529906 -1.7597614e14 +0.016110775e24 -1.2908867e34 -0.15092133e12 -0.30184265e13 -0.45276398e23 -1.2073706e1234",
                1e-4,
            ),
            (
                "1 + e1 + e2 + e12 + e3 + e13 + e23 + e123 + e4 + e14 + e24 + e124 + e34 + e134 + e234 + e1234",
                "-1.0374662 -2.0374662e1 +3.3594248e2 +3.3594248e12 -2.0374662e3 -2.0374662e13 +3.3594248e23 +3.3594248e123 -7.4343572e4 -7.4343572e14 +3.3594248e24 +3.3594248e124 +3.3594248e34 +3.3594248e134 +3.3594248e234 +3.3594248e1234",
                1e-4,
            ),
            // Coefficients sum to 91, so exp squares 8 times; f32 relative error grows to about
            // 2^8 * eps = 3e-5, i.e. 2e-3 on the largest coefficient (66). f64 matches the reference to 1e-13.
            (
                "1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + e123 + e4 + e14 + e24 + 11e124 + 12e34 + 13e134 + 14e234 + 15e1234",
                "-0.55194051 -1.1741264e1 +2.4492923e2 +2.6862086e12 +1.2298737e3 +2.0312071e13 +3.3969575e23 +4.1677590e123 -65.721195e4 -58.232090e14 +11.593442e24 +8.7586579e124 +5.9700175e34 +15.633456e134 +31.371276e234 +35.317598e1234",
                1e-2,
            ),
        ] {
            let actual = parse::<Multivector<f32, Dyn4>>("3,0,1", input).exp();
            assert_approx_eq!(
                actual,
                parse::<Multivector<f32, Dyn4>>("3,0,1", expected),
                tol
            );
        }
    }
}

mod normalization {
    use super::*;

    #[test]
    fn scalar_norm_matches_golden_value() {
        let a = parse::<Multivector<f32, Dyn3>>(
            "3,0,0",
            "1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123",
        );
        assert_eq!(a.norm(), 11.874342);
    }

    #[test]
    fn scalar_normalized_matches_golden_values() {
        for (input, expected) in [
            (
                "1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123",
                "0.084215194 + 0.084215194e1 + 0.16843039e2 + 0.25264558e12 + 0.33686078e3 + 0.42107597e13 + 0.50529116e23 + 0.5895064e123",
            ),
            (
                "1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 1e123",
                "0.10369517 +0.103695172e1 +0.20739034e2 +0.41478068e3 +0.31108552e12 +0.51847583e13 +0.62217104e23 +0.10369517e123",
            ),
        ] {
            let actual = parse::<Multivector<f32, Dyn3>>("3,0,0", input).normalized();
            assert_approx_eq!(
                actual,
                parse::<Multivector<f32, Dyn3>>("3,0,0", expected),
                1e-7
            );
        }

        let a = parse::<Multivector<f32, Dyn4>>(
            "4,0,0",
            "1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e34 + 13e134 + 14e234 + 15e1234",
        );
        assert_eq!(
            a.normalized(),
            parse(
                "4,0,0",
                "0.028386649 + 0.028386649e1 + 0.056773297e2 + 0.08515994e12 + 0.113546595e3 + 0.14193325e13 + 0.17031989e23 + 0.19870654e123 + 0.22709319e4 + 0.25547984e14 + 0.2838665e24 + 0.31225315e124 + 0.34063977e34 + 0.36902642e134 + 0.39741307e234 + 0.42579973e1234",
            )
        );

        let a = parse::<Multivector<f32, Dyn4>>(
            "3,0,1",
            "1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e34 + 13e134 + 14e234 + 15e1234",
        );
        assert_approx_eq!(
            a.normalized(),
            parse::<Multivector<f32, Dyn4>>(
                "3,0,1",
                "0.084215194 +0.67372155e4 +0.084215194e1 +0.168430388e2 +0.3368607759e3 +0.7579367e14 +0.84215194e24 +1.0105823e34 +0.252645582e12 +0.4210759699e13 +0.5052911639e23 +0.92636716e124 +1.0947975e134 +1.1790127e234 +0.5895063281e123 +1.2632279e1234",
            ),
            1e-7
        );
    }
}
