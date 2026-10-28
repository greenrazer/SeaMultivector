//! Golden values in every signature (p, q, r) from 3 to 6 dimensions.
//!
//! `tests/fixtures/golden/p_q_r.txt` holds random inputs `a` and `b` and the expected result of
//! every operation, computed once by an independent implementation. Each test reads its
//! signature's file, runs every operation in Seagrass, once with the compile-time
//! `Signature<P, Q, R>` and once with the runtime `Dynamic`, and compares the results.
//!
//! The first case of every file multiplies two values that store every blade, so each operation
//! is checked on every pair of blades. Coefficients are small integers, so the expected results
//! are exact and any difference is an error, not rounding.
//!
//! The stored values follow Seagrass's definitions: `a | b` is the left plus the right contraction
//! with the scalar part counted once, the duals in degenerate signatures are the right and left
//! complements, the right complement r satisfies r ^ x = I, and `a & b` is the undual of the
//! duals' outer product.

mod common;

use common::{DenseArray, TOL_F64, assert_close_in, every_signature};
use seagrass_math::prelude::*;

/// One case of a fixture file: its number and its `name = value` lines.
struct Case<'a> {
    number: usize,
    lines: Vec<(&'a str, &'a str)>,
}

impl<'a> Case<'a> {
    fn value(&self, name: &str) -> &'a str {
        self.lines
            .iter()
            .find(|(n, _)| *n == name)
            .map(|(_, value)| *value)
            .unwrap_or_else(|| panic!("case {} has no `{name}` line", self.number))
    }
}

/// The cases of a fixture file, skipping blank lines and `#` comments.
fn cases(fixture: &str) -> Vec<Case<'_>> {
    let mut cases: Vec<Case<'_>> = Vec::new();
    for line in fixture.lines().map(str::trim) {
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        if let Some(number) = line.strip_prefix("case ") {
            let number = number.parse().expect("a case number");
            cases.push(Case {
                number,
                lines: Vec::new(),
            });
            continue;
        }
        let (name, value) = line.split_once('=').expect("a `name = value` line");
        cases
            .last_mut()
            .expect("a `case` line before the first value")
            .lines
            .push((name.trim(), value.trim()));
    }
    assert!(!cases.is_empty(), "the fixture has no cases");
    cases
}

type Binary<M> = fn(&Multivector<f64, M>, &Multivector<f64, M>) -> Multivector<f64, M>;
type Unary<M> = fn(&Multivector<f64, M>) -> Multivector<f64, M>;

/// Runs every case of a signature's fixture in `metric` and compares each result with the stored
/// one.
fn check<M: Metric>(metric: M, fixture: &str) {
    let prefix = format!(
        "{},{},{}>",
        metric.positive(),
        metric.negative(),
        metric.zero()
    );
    let parse = |value: &str| -> Multivector<f64, M> {
        format!("{prefix}{value}")
            .parse()
            .unwrap_or_else(|e| panic!("cannot parse `{value}`: {e:?}"))
    };
    let binary: [(&str, Binary<M>); 7] = [
        ("a * b", |a, b| a * b),
        ("a ^ b", |a, b| a ^ b),
        ("a << b", |a, b| a << b),
        ("a >> b", |a, b| a >> b),
        ("a | b", |a, b| a | b),
        ("a & b", |a, b| a & b),
        ("a b ã", |a, b| a.sandwich_product(b)),
    ];
    let unary: [(&str, Unary<M>); 7] = [
        ("reverse", Multivector::reversed),
        ("involute", Multivector::involuted),
        ("conjugate", Multivector::conjugated),
        ("dual", Multivector::dual),
        ("undual", Multivector::undual),
        ("left complement", Multivector::left_complemented),
        ("right complement", Multivector::right_complemented),
    ];
    let kind = std::any::type_name::<M>();
    for case in cases(fixture) {
        let (a, b) = (parse(case.value("a")), parse(case.value("b")));
        let context = |name: &str| {
            format!(
                "{name} in {kind}, case {}\n  a = {a}\n  b = {b}",
                case.number
            )
        };
        for (name, op) in binary {
            assert_close_in(
                &op(&a, &b),
                &parse(case.value(name)),
                TOL_F64,
                &context(name),
            );
        }
        for (name, op) in unary {
            assert_close_in(&op(&a), &parse(case.value(name)), TOL_F64, &context(name));
        }
        let mut scalar = [0.0; 64];
        scalar[0] = a.scalar_product(&b);
        let name = "scalar product";
        assert_close_in(
            &DenseArray(scalar),
            &parse(case.value(name)),
            TOL_F64,
            &context(name),
        );
    }
}

/// One test per signature, checking both its compile-time and its runtime form.
macro_rules! test_signature {
    ($name:ident, $bits:ty, $p:literal, $q:literal, $r:literal) => {
        #[test]
        fn $name() {
            let fixture = include_str!(concat!(
                env!("CARGO_MANIFEST_DIR"),
                "/tests/fixtures/golden/",
                $p,
                "_",
                $q,
                "_",
                $r,
                ".txt"
            ));
            check(Signature::<$p, $q, $r>, fixture);
            check(Dynamic::<$bits>::new($p, $q, $r).unwrap(), fixture);
        }
    };
}

every_signature!(test_signature);
