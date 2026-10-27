//! The per-call product loops that plans replaced, kept as the reference the plans are tested
//! against. Each product finds its pairs and signs with the bitset tricks on every call.

use super::*;

use crate::plan::{SandwichPlan, with_sandwich_plan};

/// The sum over the left blades x in `pairs` of `±a[x] * b[x XOR output]`, negated where
/// `signs` holds x.
fn dot<T: Coefficient, B: BladeSet>(
    (a_blades, a): (B, &[T]),
    (b_blades, b): (B, &[T]),
    pairs: B,
    signs: B,
    output: BladeIndex,
) -> T {
    let mut out = T::zero();
    for x in pairs.iter() {
        let mut c =
            a[a_blades.count_below(x) as usize] * b[b_blades.count_below(x ^ output) as usize];
        if signs.has(x) {
            c = -c;
        }
        out += c;
    }
    out
}

/// Runs one product the per-call way: `structure` gives the output blades and `pairs` the pair
/// masks and signs of each output blade.
fn per_call<T: Coefficient, M: Metric>(
    a: &Multivector<T, M>,
    b: &Multivector<T, M>,
    structure: impl Fn(M::Bits, M::Bits, ProductConstants) -> M::Bits,
    pairs: impl Fn(M::Bits, M::Bits, BladeIndex, ProductConstants) -> Vec<(M::Bits, M::Bits)>,
) -> Multivector<T, M> {
    let constants = a.constants();
    let output = structure(a.blades, b.blades, constants);
    let coefficients = output
        .iter()
        .map(|o| {
            let mut c = T::zero();
            for (n, (pairs, signs)) in pairs(a.blades, b.blades, o, constants)
                .into_iter()
                .enumerate()
            {
                let term = dot(
                    (a.blades, &a.coefficients),
                    (b.blades, &b.coefficients),
                    pairs,
                    signs,
                    o,
                );
                c = if n == 0 { term } else { c + term };
            }
            c
        })
        .collect();
    Multivector::from_raw(a.metric, output, coefficients)
}

pub(crate) fn geometric<T: Coefficient, M: Metric>(
    a: &Multivector<T, M>,
    b: &Multivector<T, M>,
) -> Multivector<T, M> {
    per_call(
        a,
        b,
        |a, b, c| a.geometric_structure(b, c.non_null_count),
        |a, b, o, c| {
            vec![(
                a & b.geometric_pairs(o, c.non_null),
                M::Bits::product_signs(o, c.negative),
            )]
        },
    )
}

pub(crate) fn outer<T: Coefficient, M: Metric>(
    a: &Multivector<T, M>,
    b: &Multivector<T, M>,
) -> Multivector<T, M> {
    per_call(
        a,
        b,
        |a, b, _| a.outer_structure(b),
        |a, b, o, _| {
            vec![(
                a & b.reverse_bits().outer_pairs(o),
                M::Bits::reorder_signs(o),
            )]
        },
    )
}

pub(crate) fn left_contraction<T: Coefficient, M: Metric>(
    a: &Multivector<T, M>,
    b: &Multivector<T, M>,
) -> Multivector<T, M> {
    per_call(
        a,
        b,
        |a, b, c| a.left_contraction_structure(b, c.non_null_count),
        |a, b, o, c| {
            vec![(
                a & b.left_contraction_pairs(o, c.non_null),
                M::Bits::product_signs(o, c.negative),
            )]
        },
    )
}

pub(crate) fn right_contraction<T: Coefficient, M: Metric>(
    a: &Multivector<T, M>,
    b: &Multivector<T, M>,
) -> Multivector<T, M> {
    per_call(
        a,
        b,
        |a, b, c| a.right_contraction_structure(b, c.non_null_count),
        |a, b, o, c| {
            vec![(
                a & b.right_contraction_pairs(o, c.non_null),
                M::Bits::product_signs(o, c.negative),
            )]
        },
    )
}

pub(crate) fn inner<T: Coefficient, M: Metric>(
    a: &Multivector<T, M>,
    b: &Multivector<T, M>,
) -> Multivector<T, M> {
    per_call(
        a,
        b,
        |a, b, c| a.inner_structure(b, c.non_null_count),
        |a, b, o, c| {
            let signs = M::Bits::product_signs(o, c.negative);
            let mut out = vec![(a & b.left_contraction_pairs(o, c.non_null), signs)];
            if o != BladeIndex::SCALAR {
                out.push((a & b.right_contraction_pairs(o, c.non_null), signs));
            }
            out
        },
    )
}

/// The two-stage sandwich, with the reverse of `m` folded into the second product's signs and
/// numerically zero blades removed.
pub(crate) fn sandwich<T: Coefficient, M: Metric>(
    m: &Multivector<T, M>,
    x: &Multivector<T, M>,
) -> Multivector<T, M> {
    let first = geometric(m, x);
    let constants = m.constants();
    let output = first
        .blades
        .geometric_structure(m.blades, constants.non_null_count);
    let coefficients = output
        .iter()
        .map(|o| {
            let pairs = first.blades & m.blades.geometric_pairs(o, constants.non_null);
            let signs = M::Bits::product_signs(o, constants.negative)
                ^ M::Bits::REVERSE_SIGNS.geometric_pairs(o, constants.non_null);
            dot(
                (first.blades, &first.coefficients),
                (m.blades, &m.coefficients),
                pairs,
                signs,
                o,
            )
        })
        .collect();
    let mut out = Multivector::from_raw(m.metric, output, coefficients);
    out.remove_zeros();
    out
}

/// The flattened sandwich plan's result, without the removal of numeric zeros.
pub(crate) fn sandwich_plan_only<T: Coefficient, M: Metric>(
    m: &Multivector<T, M>,
    x: &Multivector<T, M>,
) -> Multivector<T, M> {
    with_sandwich_plan(&m.metric, m.blades, x.blades, |plan: &SandwichPlan| {
        Multivector::from_raw(
            m.metric,
            plan.output(),
            plan.run(&m.coefficients, &x.coefficients),
        )
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::metric::Dynamic;
    use crate::traits::{
        GeometricProduct, InnerProduct, LeftContraction, OuterProduct, RightContraction,
        SandwichProduct,
    };

    /// A small deterministic generator, so failures reproduce.
    struct XorShift(u64);

    impl XorShift {
        fn next(&mut self) -> u64 {
            self.0 ^= self.0 << 13;
            self.0 ^= self.0 >> 7;
            self.0 ^= self.0 << 17;
            self.0
        }

        /// A coefficient in [-2, 2) with a few significant bits, so products are rarely exact.
        fn coefficient(&mut self) -> f64 {
            (self.next() % 4_000_001) as f64 / 1_000_000.0 - 2.0
        }

        /// A random multivector; the density of stored blades varies from case to case.
        fn multivector<B: BladeSet>(&mut self, metric: Dynamic<B>) -> Multivector<f64, Dynamic<B>> {
            let mut blades = self.next();
            for _ in 0..self.next() % 3 {
                blades &= self.next();
            }
            let blades = B::from_u64(blades);
            let coefficients = blades.iter().map(|_| self.coefficient()).collect();
            Multivector::from_raw(metric, blades, coefficients)
        }
    }

    fn assert_identical<B: BladeSet>(
        plan: &Multivector<f64, Dynamic<B>>,
        reference: &Multivector<f64, Dynamic<B>>,
        what: &str,
    ) {
        assert_eq!(plan.blades, reference.blades, "{what}: stored blades");
        for (p, r) in plan.coefficients.iter().zip(&reference.coefficients) {
            assert_eq!(p.to_bits(), r.to_bits(), "{what}: {p} != {r}");
        }
    }

    /// Values agree blade by blade within `tol` times the largest coefficient; blades only one
    /// side stores must be within `tol` of zero.
    fn assert_close<B: BladeSet>(
        plan: &Multivector<f64, Dynamic<B>>,
        reference: &Multivector<f64, Dynamic<B>>,
        tol: f64,
        what: &str,
    ) {
        let scale = reference
            .coefficients
            .iter()
            .fold(1.0_f64, |m, c| m.max(c.abs()));
        for index in (plan.blades | reference.blades).iter() {
            let (p, r) = (plan.get(index), reference.get(index));
            assert!(
                (p - r).abs() <= tol * scale,
                "{what}: blade {index}: {p} != {r}"
            );
        }
    }

    fn check_every_signature<B: BladeSet>(cases: usize) {
        let mut rng = XorShift(0x9E37_79B9_7F4A_7C15 ^ u64::from(B::DIM));
        let dim = B::DIM;
        for p in 0..=dim {
            for q in 0..=dim - p {
                let metric = Dynamic::<B>::new(p, q, dim - p - q).unwrap();
                for case in 0..cases {
                    let (a, b) = (rng.multivector(metric), rng.multivector(metric));
                    let what = |op: &str| format!("{op} in ({p},{q},{}) case {case}", dim - p - q);
                    assert_identical(
                        &a.geometric_product(&b),
                        &geometric(&a, &b),
                        &what("geometric"),
                    );
                    assert_identical(&a.outer_product(&b), &outer(&a, &b), &what("outer"));
                    assert_identical(
                        &a.left_contraction(&b),
                        &left_contraction(&a, &b),
                        &what("left"),
                    );
                    assert_identical(
                        &a.right_contraction(&b),
                        &right_contraction(&a, &b),
                        &what("right"),
                    );
                    // The plan sums the left and right contraction pairs in one running sum.
                    let inner_plan = a.inner_product(&b);
                    assert_eq!(inner_plan.blades, inner(&a, &b).blades, "{}", what("inner"));
                    assert_close(&inner_plan, &inner(&a, &b), 1e-14, &what("inner"));
                    // The flattened sandwich multiplies in a different order and drops blades
                    // whose terms cancel symbolically.
                    let sandwich_plan = a.sandwich_product(&b);
                    let reference = sandwich(&a, &b);
                    assert_close(&sandwich_plan, &reference, 1e-13, &what("sandwich"));
                    assert!(
                        (sandwich_plan_only(&a, &b).blades
                            & !geometric(&geometric(&a, &b), &a).blades)
                            .is_empty(),
                        "{}",
                        what("sandwich structure")
                    );
                }
            }
        }
    }

    #[test]
    fn plans_match_the_per_call_products_in_every_3d_signature() {
        check_every_signature::<u8>(200);
    }

    #[test]
    fn plans_match_the_per_call_products_in_every_4d_signature() {
        check_every_signature::<u16>(120);
    }

    #[test]
    fn plans_match_the_per_call_products_in_every_5d_signature() {
        check_every_signature::<u32>(40);
    }

    #[test]
    fn plans_match_the_per_call_products_in_every_6d_signature() {
        check_every_signature::<u64>(12);
    }
}
