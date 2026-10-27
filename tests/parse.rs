//! Parsing tests: metric strings, basis blade strings and runtime-metric multivector strings,
//! including the exact error variant for every rejected input.

use seagrass_math::prelude::*;

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

/// Expected `blades!` output for a multivector that stores every blade of its algebra
/// (a basis blade set with all bits set): blade `i` holds `coefficients[i]`.
fn dense(coefficients: &[f32]) -> Vec<(u8, f32)> {
    (0..).zip(coefficients.iter().copied()).collect()
}

/// Tests ported from `src/lib.rs`.
mod lib {
    use super::*;

    #[test]
    fn parses_through_4d_alias() {
        let mv = Multivector::<f32, Dyn4>::from_str("3,1,0>e1 + e12");
        assert!(mv.is_ok());
    }
}

/// Tests ported from `src/multivector/metric.rs`.
mod metric {
    use super::*;

    /// Missing counts are zero. Signatures must match the dimension of the bitset type, so the
    /// cases use 3D and 6D runtime metrics.
    #[test]
    fn parses_positive_negative_and_zero_counts() {
        let cases = [("3", (3, 0, 0)), ("1, 2", (1, 2, 0))];
        for (input, (p, q, r)) in cases {
            assert_eq!(
                Dyn3::from_str(input).unwrap(),
                Dyn3::new(p, q, r).unwrap(),
                "{input:?}"
            );
        }
        assert_eq!(
            Dyn6::from_str("1, 2, 3").unwrap(),
            Dyn6::new(1, 2, 3).unwrap()
        );
    }

    #[test]
    fn rejects_malformed_metric_strings() {
        let cases = [
            ("", MetricError::MissingNumPositive),
            ("a", MetricError::InvalidNumPositive),
            ("1,a", MetricError::InvalidNumNegative),
            ("1,2,c", MetricError::InvalidNumZero),
            ("1,2,3,4", MetricError::InvalidMetricLength),
        ];
        for (input, expected) in cases {
            assert_eq!(Dyn3::from_str(input).unwrap_err(), expected, "{input:?}");
        }
    }
}

/// Tests ported from `src/multivector/basis_blade.rs`.
mod basis_blade {
    use super::*;

    #[test]
    fn parses_coefficient_and_blade_index() {
        let cases = [("3e12", (3, 3_f32)), ("3e5", (16, 3_f32))];
        for (input, expected) in cases {
            let a = Multivector::<f32, Dyn6>::from_str(&format!("6>{input}")).unwrap();
            assert_eq!(blades!(&a), [expected], "{input:?}");
        }
    }

    /// The parser reports positions within the whole input, which here starts with `6>`.
    #[test]
    fn rejects_malformed_blades() {
        let cases = [
            ("e", ParseError::InvalidBasisBladeIndex),
            ("ae3", ParseError::UnexpectedToken { position: 2 }),
            ("ea", ParseError::InvalidBasisVector),
            ("e0", ParseError::InvalidBasisVectorZero),
            ("e21", ParseError::InvalidBasisVectorOrdering),
            ("e8", ParseError::InvalidBasisVectorOutOfRange),
            ("e11", ParseError::DuplicateBasisVector),
        ];
        for (input, expected) in cases {
            assert_eq!(
                Multivector::<f32, Dyn6>::from_str(&format!("6>{input}")).unwrap_err(),
                expected,
                "{input:?}"
            );
        }
    }
}

/// Tests ported from `src/multivector/low_dimensional/custom_metric.rs`.
mod custom_metric {
    use super::*;

    #[test]
    fn parses_metric_and_every_blade() {
        let a =
            Multivector::<f32, Dyn3>::from_str("3>1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123")
                .unwrap();
        assert_eq!(a.metric(), Dyn3::new(3, 0, 0).unwrap());
        assert_eq!(
            blades!(&a),
            dense(&[1_f32, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0])
        );

        let cases = [
            (
                "2,2>1 - e1 + 2e2 - 3e12 + 4e3 - 5e13 + 6e23 - 7e123 + 8e4 - 9e14 + 10e24 - 11e124 + 12e34 - 13e134 + 14e234 - 15e1234",
                Dyn4::new(2, 2, 0).unwrap(),
            ),
            (
                "0,0,4>1 - e1 + 2e2 - 3e12 + 4e3 - 5e13 + 6e23 - 7e123 + 8e4 - 9e14 + 10e24 - 11e124 + 12e34 - 13e134 + 14e234 - 15e1234",
                Dyn4::new(0, 0, 4).unwrap(),
            ),
        ];
        for (input, metric) in cases {
            let a = Multivector::<f32, Dyn4>::from_str(input).unwrap();
            assert_eq!(a.metric(), metric, "{input:?}");
            assert_eq!(
                blades!(&a),
                dense(&[
                    1_f32, -1.0, 2.0, -3.0, 4.0, -5.0, 6.0, -7.0, 8.0, -9.0, 10.0, -11.0, 12.0,
                    -13.0, 14.0, -15.0
                ]),
                "{input:?}"
            );
        }
    }

    #[test]
    fn rejects_missing_or_invalid_metric() {
        let cases = [
            ("e3", ParseError::MissingMetric),
            (
                "a>e3",
                ParseError::InvalidMetric(MetricError::InvalidNumPositive),
            ),
            ("3,4,5>e1", ParseError::InvalidMetricDimensions),
        ];
        for (input, expected) in cases {
            assert_eq!(
                Multivector::<f32, Dyn3>::from_str(input).unwrap_err(),
                expected,
                "{input:?}"
            );
        }
    }
}
