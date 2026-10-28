//! Every signature (p, q, r) from 3 to 6 dimensions, 74 in all, checked against products worked
//! out by hand.
//!
//! The reference here shares no code with the library. It multiplies two basis blades by writing
//! out their vectors, sorting them into ascending order (each swap of two different vectors flips
//! the sign) and replacing each repeated vector by its square. Every other operation is then
//! defined from that product: the outer product, the contractions and the inner product as grade
//! parts, the duals and complements from their definitions, the regressive product through the
//! duals, and the sandwich as m x m̃. Each signature runs as a compile-time `Signature<P, Q, R>`
//! and as a runtime `Dynamic`, on random values from dense to a single blade, so a wrong sign or
//! metric factor for any pair of blades shows up as a wrong coefficient.
//!
//! Coefficients are small integers, so both sides compute exact integers whatever their summation
//! order, and a wrong sign or factor moves a coefficient far outside the shared tolerance.
//! `tests/golden.rs` makes the same comparison against stored golden values.

mod common;

use common::{DenseArray, TOL_F64, assert_close_in, every_signature};
use seagrass_math::prelude::*;

/// A signature and the product of every pair of its basis blades, worked out by hand.
struct Algebra {
    dim: u32,
    /// Whether the signature has null vectors, which changes how the duals are defined.
    degenerate: bool,
    /// `factors[a * blades + b]` is f in e_a e_b = f e_(a xor b), for blades named by bitmasks.
    factors: Vec<f64>,
}

impl Algebra {
    /// Basis vectors are numbered positive, then negative, then null, as in the library.
    fn new(p: u8, q: u8, r: u8) -> Self {
        let squares: Vec<f64> = (0..p + q + r)
            .map(|k| match k {
                k if k < p => 1.0,
                k if k < p + q => -1.0,
                _ => 0.0,
            })
            .collect();
        let blades = 1 << squares.len();
        let factors = (0..blades * blades)
            .map(|i| multiply_blades(&squares, i / blades, i % blades))
            .collect();
        Self {
            dim: squares.len() as u32,
            degenerate: r > 0,
            factors,
        }
    }

    fn blades(&self) -> usize {
        1 << self.dim
    }

    fn factor(&self, a: usize, b: usize) -> f64 {
        self.factors[a * self.blades() + b]
    }

    /// The bilinear product that keeps the part of e_a e_b for which `keep(grade a, grade b,
    /// grade of e_a e_b)` holds.
    fn product(&self, x: &[f64], y: &[f64], keep: impl Fn(u32, u32, u32) -> bool) -> Vec<f64> {
        let mut out = vec![0.0; self.blades()];
        for a in 0..self.blades() {
            for b in 0..self.blades() {
                let grades = (a.count_ones(), b.count_ones(), (a ^ b).count_ones());
                if keep(grades.0, grades.1, grades.2) {
                    out[a ^ b] += self.factor(a, b) * x[a] * y[b];
                }
            }
        }
        out
    }

    fn geometric(&self, x: &[f64], y: &[f64]) -> Vec<f64> {
        self.product(x, y, |_, _, _| true)
    }

    fn outer(&self, x: &[f64], y: &[f64]) -> Vec<f64> {
        self.product(x, y, |a, b, out| out == a + b)
    }

    fn left_contraction(&self, x: &[f64], y: &[f64]) -> Vec<f64> {
        self.product(x, y, |a, b, out| b >= a && out == b - a)
    }

    fn right_contraction(&self, x: &[f64], y: &[f64]) -> Vec<f64> {
        self.product(x, y, |a, b, out| a >= b && out == a - b)
    }

    /// The grade |a - b| part, which counts the scalar part of two equal grades once.
    fn inner(&self, x: &[f64], y: &[f64]) -> Vec<f64> {
        self.product(x, y, |a, b, out| out == a.abs_diff(b))
    }

    /// Sends each blade e_a to `f(a) e_(target(a))`.
    fn map(
        &self,
        x: &[f64],
        target: impl Fn(usize) -> usize,
        f: impl Fn(usize) -> f64,
    ) -> Vec<f64> {
        let mut out = vec![0.0; self.blades()];
        for a in 0..self.blades() {
            out[target(a)] += f(a) * x[a];
        }
        out
    }

    /// Flips the sign of the blades of grade k for which `flip(k)` holds.
    fn by_grade(&self, x: &[f64], flip: impl Fn(u32) -> bool) -> Vec<f64> {
        self.map(x, |a| a, |a| if flip(a.count_ones()) { -1.0 } else { 1.0 })
    }

    fn reverse(&self, x: &[f64]) -> Vec<f64> {
        self.by_grade(x, |k| k * k.saturating_sub(1) / 2 % 2 == 1)
    }

    fn involute(&self, x: &[f64]) -> Vec<f64> {
        self.by_grade(x, |k| k % 2 == 1)
    }

    fn conjugate(&self, x: &[f64]) -> Vec<f64> {
        self.by_grade(x, |k| k * (k + 1) / 2 % 2 == 1)
    }

    fn pseudoscalar(&self) -> usize {
        self.blades() - 1
    }

    /// Each blade moved to its complement with no sign change.
    fn poincare_dual(&self, x: &[f64]) -> Vec<f64> {
        self.map(x, |a| self.pseudoscalar() ^ a, |_| 1.0)
    }

    /// r with r ^ e_a = I: the complement e_c times the sign that orders e_c e_a into I.
    fn right_complement(&self, x: &[f64]) -> Vec<f64> {
        let full = self.pseudoscalar();
        self.map(x, |a| full ^ a, |a| self.factor(full ^ a, a))
    }

    /// l with e_a ^ l = I.
    fn left_complement(&self, x: &[f64]) -> Vec<f64> {
        let full = self.pseudoscalar();
        self.map(x, |a| full ^ a, |a| self.factor(a, full ^ a))
    }

    /// x I⁻¹, or the right complement when null vectors leave I without an inverse.
    fn dual(&self, x: &[f64]) -> Vec<f64> {
        if self.degenerate {
            return self.right_complement(x);
        }
        // I I is ±1, so I⁻¹ = (I I) I.
        let full = self.pseudoscalar();
        let mut inverse = vec![0.0; self.blades()];
        inverse[full] = self.factor(full, full);
        self.geometric(x, &inverse)
    }

    /// x I, or the left complement when null vectors are present.
    fn undual(&self, x: &[f64]) -> Vec<f64> {
        if self.degenerate {
            return self.left_complement(x);
        }
        let mut pseudoscalar = vec![0.0; self.blades()];
        pseudoscalar[self.pseudoscalar()] = 1.0;
        self.geometric(x, &pseudoscalar)
    }

    fn regressive(&self, x: &[f64], y: &[f64]) -> Vec<f64> {
        self.undual(&self.outer(&self.dual(x), &self.dual(y)))
    }

    fn sandwich(&self, m: &[f64], x: &[f64]) -> Vec<f64> {
        self.geometric(&self.geometric(m, x), &self.reverse(m))
    }
}

/// f in e_a e_b = f e_(a xor b), by hand: write out the vectors of e_a then e_b, bubble-sort them
/// (swapping two different vectors flips the sign), then replace each vector that now appears
/// twice in a row by its square.
fn multiply_blades(squares: &[f64], a: usize, b: usize) -> f64 {
    let vectors = |blade: usize| (0..squares.len()).filter(move |k| blade >> k & 1 == 1);
    let mut word: Vec<usize> = vectors(a).chain(vectors(b)).collect();
    let mut factor = 1.0;
    for end in (1..word.len()).rev() {
        for i in 0..end {
            if word[i] > word[i + 1] {
                word.swap(i, i + 1);
                factor = -factor;
            }
        }
    }
    for pair in word.windows(2) {
        if pair[0] == pair[1] {
            factor *= squares[pair[0]];
        }
    }
    factor
}

/// xorshift64, seeded per signature so a failure reproduces.
struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }

    fn below(&mut self, n: usize) -> usize {
        (self.next() % n as u64) as usize
    }

    /// A non-zero integer from -4 to 4.
    fn coefficient(&mut self) -> f64 {
        let c = self.below(4) as f64 + 1.0;
        if self.next() & 1 == 0 { c } else { -c }
    }

    /// Coefficients for one of five shapes: every blade, about half of them, one grade, the even
    /// grades (the shape of a rotor or motor), or one to three blades.
    fn value(&mut self, algebra: &Algebra) -> Vec<f64> {
        let shape = self.below(5);
        let grade = self.below(algebra.dim as usize + 1) as u32;
        let few = [0; 3].map(|_| self.below(algebra.blades()));
        let count = 1 + self.below(3);
        (0..algebra.blades())
            .map(|blade| {
                let stored = match shape {
                    0 => true,
                    1 => self.next() & 1 == 0,
                    2 => blade.count_ones() == grade,
                    3 => blade.count_ones() % 2 == 0,
                    _ => few[..count].contains(&blade),
                };
                if stored { self.coefficient() } else { 0.0 }
            })
            .collect()
    }
}

fn value<M: Metric>(metric: M, coefficients: &[f64]) -> Multivector<f64, M> {
    let blades = coefficients
        .iter()
        .enumerate()
        .filter(|(_, c)| **c != 0.0)
        .map(|(blade, c)| (BladeIndex::new(blade as u8), *c));
    Multivector::from_blades(metric, blades).unwrap()
}

/// Asserts that `actual` agrees with the reference coefficients `expected`, counting blades it
/// does not store as zero.
#[track_caller]
fn assert_matches<M: Metric>(actual: &Multivector<f64, M>, expected: &[f64], what: &str) {
    let mut dense = [0.0; 64];
    dense[..expected.len()].copy_from_slice(expected);
    assert_close_in(actual, &DenseArray(dense), TOL_F64, what);
}

/// Runs every operation on random values in `metric` against the hand-worked reference.
fn check<M: Metric>(metric: M) {
    let (p, q, r) = (metric.positive(), metric.negative(), metric.zero());
    let algebra = Algebra::new(p, q, r);
    // Fewer cases in higher dimensions, where each value already holds up to 64 blades.
    let cases = [0, 0, 0, 60, 40, 20, 10][algebra.dim as usize];
    let mut rng =
        Rng(0x9E37_79B9_7F4A_7C15 ^ (u64::from(p) << 16 | u64::from(q) << 8 | u64::from(r)));
    for case in 0..cases {
        let (x, y) = (rng.value(&algebra), rng.value(&algebra));
        let (a, b) = (value(metric, &x), value(metric, &y));
        let what = |op: &str| {
            format!(
                "{op} in {}, case {case}\n  a = {a}\n  b = {b}",
                std::any::type_name::<M>()
            )
        };

        assert_matches(&(&a * &b), &algebra.geometric(&x, &y), &what("a * b"));
        assert_matches(&(&a ^ &b), &algebra.outer(&x, &y), &what("a ^ b"));
        assert_matches(
            &(&a << &b),
            &algebra.left_contraction(&x, &y),
            &what("a << b"),
        );
        assert_matches(
            &(&a >> &b),
            &algebra.right_contraction(&x, &y),
            &what("a >> b"),
        );
        assert_matches(&(&a | &b), &algebra.inner(&x, &y), &what("a | b"));
        assert_matches(&(&a & &b), &algebra.regressive(&x, &y), &what("a & b"));
        assert_matches(
            &a.sandwich_product(&b),
            &algebra.sandwich(&x, &y),
            &what("a b ã"),
        );
        let (mut actual, mut expected) = ([0.0; 64], [0.0; 64]);
        actual[0] = a.scalar_product(&b);
        expected[0] = algebra.geometric(&x, &y)[0];
        assert_close_in(
            &DenseArray(actual),
            &DenseArray(expected),
            TOL_F64,
            &what("scalar product"),
        );

        assert_matches(&a.reversed(), &algebra.reverse(&x), &what("reverse"));
        assert_matches(&a.involuted(), &algebra.involute(&x), &what("involution"));
        assert_matches(&a.conjugated(), &algebra.conjugate(&x), &what("conjugate"));
        assert_matches(&a.dual(), &algebra.dual(&x), &what("dual"));
        assert_matches(&a.undual(), &algebra.undual(&x), &what("undual"));
        assert_matches(
            &a.poincare_dual(),
            &algebra.poincare_dual(&x),
            &what("Poincaré dual"),
        );
        assert_matches(
            &a.right_complemented(),
            &algebra.right_complement(&x),
            &what("right complement"),
        );
        assert_matches(
            &a.left_complemented(),
            &algebra.left_complement(&x),
            &what("left complement"),
        );
    }
}

/// Two dense 6D values give a sandwich of about 64³ / 2 terms, more than a 16-bit offset counts;
/// such plans once indexed past their terms.
#[test]
fn sandwich_of_two_dense_6d_values() {
    let metric = Signature::<3, 3, 0>;
    let algebra = Algebra::new(3, 3, 0);
    let m: Vec<f64> = (0..64)
        .map(|i| [1.0, -2.0, 3.0, -1.0, 2.0][i % 5])
        .collect();
    let x: Vec<f64> = (0..64).map(|i| [2.0, -1.0, 1.0][i % 3]).collect();
    assert_matches(
        &value(metric, &m).sandwich_product(&value(metric, &x)),
        &algebra.sandwich(&m, &x),
        "a b ã for dense a and b in (3, 3, 0)",
    );
}

/// One test per signature, checking both its compile-time and its runtime form.
macro_rules! test_signature {
    ($name:ident, $bits:ty, $p:literal, $q:literal, $r:literal) => {
        #[test]
        fn $name() {
            check(Signature::<$p, $q, $r>);
            check(Dynamic::<$bits>::new($p, $q, $r).unwrap());
        }
    };
}

every_signature!(test_signature);
