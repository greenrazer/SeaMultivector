//! Helpers shared by the integration tests.
//!
//! Each file in `tests/` is its own crate and uses a different subset of these helpers.
#![allow(dead_code)]

use std::fmt::Debug;

use seagrass_math::prelude::*;

/// Tolerance for computed `f64` results.
pub(crate) const TOL_F64: f64 = 1e-12;

/// Tolerance for computed `f32` results.
pub(crate) const TOL_F32: f64 = 1e-6;

/// The coefficient of every basis blade, with zero for blades that are not stored.
pub(crate) trait Dense {
    fn dense(&self) -> [f64; 64];
}

impl<T: Coefficient + Into<f64>, M: Metric> Dense for Multivector<T, M> {
    fn dense(&self) -> [f64; 64] {
        let mut out = [0.0; 64];
        for (index, c) in self.iter() {
            out[index.bits() as usize] = (*c).into();
        }
        out
    }
}

impl<T: Coefficient + Into<f64>, M: Metric, K: Kind> Dense for Typed<T, M, K> {
    fn dense(&self) -> [f64; 64] {
        self.as_multivector().dense()
    }
}

/// A dense array of coefficients, for comparisons against values computed by hand.
#[derive(Debug)]
pub(crate) struct DenseArray(pub(crate) [f64; 64]);

impl Dense for DenseArray {
    fn dense(&self) -> [f64; 64] {
        self.0
    }
}

/// Asserts that two values agree blade by blade, treating blades that are not stored as zero.
///
/// Each coefficient may differ by `tol` times the largest expected coefficient, or by `tol` itself
/// when every expected coefficient is below 1. Stored zeros are ignored, so `1 + 0e12` is close
/// to `1`.
#[track_caller]
pub(crate) fn assert_close<A: Dense + Debug, E: Dense + Debug>(actual: &A, expected: &E, tol: f64) {
    assert_close_in(actual, expected, tol, "");
}

/// [`assert_close`], with `context` (such as the operation and case) heading the failure message.
#[track_caller]
pub(crate) fn assert_close_in<A: Dense + Debug, E: Dense + Debug>(
    actual: &A,
    expected: &E,
    tol: f64,
    context: &str,
) {
    let (a, e) = (actual.dense(), expected.dense());
    let scale = e.iter().fold(1.0_f64, |m, c| m.max(c.abs()));
    for blade in 0..64 {
        let diff = (a[blade] - e[blade]).abs();
        assert!(
            diff <= tol * scale,
            "{context}\nblade {blade}: {} != {} (|diff| = {diff:e}, allowed {:e})\nactual:   {actual:?}\nexpected: {expected:?}",
            a[blade],
            e[blade],
            tol * scale,
        );
    }
}

/// Invokes `$each!(name, bits, p, q, r)` for every signature (p, q, r) from 3 to 6 dimensions,
/// 74 in all, where `bits` is the blade-set type of that dimension and `name` is
/// `signature_p_q_r`.
// Unused in the test files that do not loop over signatures, as is the re-export below.
#[allow(unused_macros)]
macro_rules! every_signature {
    ($each:ident) => {
        $each!(signature_3_0_0, u8, 3, 0, 0);
        $each!(signature_2_1_0, u8, 2, 1, 0);
        $each!(signature_2_0_1, u8, 2, 0, 1);
        $each!(signature_1_2_0, u8, 1, 2, 0);
        $each!(signature_1_1_1, u8, 1, 1, 1);
        $each!(signature_1_0_2, u8, 1, 0, 2);
        $each!(signature_0_3_0, u8, 0, 3, 0);
        $each!(signature_0_2_1, u8, 0, 2, 1);
        $each!(signature_0_1_2, u8, 0, 1, 2);
        $each!(signature_0_0_3, u8, 0, 0, 3);

        $each!(signature_4_0_0, u16, 4, 0, 0);
        $each!(signature_3_1_0, u16, 3, 1, 0);
        $each!(signature_3_0_1, u16, 3, 0, 1);
        $each!(signature_2_2_0, u16, 2, 2, 0);
        $each!(signature_2_1_1, u16, 2, 1, 1);
        $each!(signature_2_0_2, u16, 2, 0, 2);
        $each!(signature_1_3_0, u16, 1, 3, 0);
        $each!(signature_1_2_1, u16, 1, 2, 1);
        $each!(signature_1_1_2, u16, 1, 1, 2);
        $each!(signature_1_0_3, u16, 1, 0, 3);
        $each!(signature_0_4_0, u16, 0, 4, 0);
        $each!(signature_0_3_1, u16, 0, 3, 1);
        $each!(signature_0_2_2, u16, 0, 2, 2);
        $each!(signature_0_1_3, u16, 0, 1, 3);
        $each!(signature_0_0_4, u16, 0, 0, 4);

        $each!(signature_5_0_0, u32, 5, 0, 0);
        $each!(signature_4_1_0, u32, 4, 1, 0);
        $each!(signature_4_0_1, u32, 4, 0, 1);
        $each!(signature_3_2_0, u32, 3, 2, 0);
        $each!(signature_3_1_1, u32, 3, 1, 1);
        $each!(signature_3_0_2, u32, 3, 0, 2);
        $each!(signature_2_3_0, u32, 2, 3, 0);
        $each!(signature_2_2_1, u32, 2, 2, 1);
        $each!(signature_2_1_2, u32, 2, 1, 2);
        $each!(signature_2_0_3, u32, 2, 0, 3);
        $each!(signature_1_4_0, u32, 1, 4, 0);
        $each!(signature_1_3_1, u32, 1, 3, 1);
        $each!(signature_1_2_2, u32, 1, 2, 2);
        $each!(signature_1_1_3, u32, 1, 1, 3);
        $each!(signature_1_0_4, u32, 1, 0, 4);
        $each!(signature_0_5_0, u32, 0, 5, 0);
        $each!(signature_0_4_1, u32, 0, 4, 1);
        $each!(signature_0_3_2, u32, 0, 3, 2);
        $each!(signature_0_2_3, u32, 0, 2, 3);
        $each!(signature_0_1_4, u32, 0, 1, 4);
        $each!(signature_0_0_5, u32, 0, 0, 5);

        $each!(signature_6_0_0, u64, 6, 0, 0);
        $each!(signature_5_1_0, u64, 5, 1, 0);
        $each!(signature_5_0_1, u64, 5, 0, 1);
        $each!(signature_4_2_0, u64, 4, 2, 0);
        $each!(signature_4_1_1, u64, 4, 1, 1);
        $each!(signature_4_0_2, u64, 4, 0, 2);
        $each!(signature_3_3_0, u64, 3, 3, 0);
        $each!(signature_3_2_1, u64, 3, 2, 1);
        $each!(signature_3_1_2, u64, 3, 1, 2);
        $each!(signature_3_0_3, u64, 3, 0, 3);
        $each!(signature_2_4_0, u64, 2, 4, 0);
        $each!(signature_2_3_1, u64, 2, 3, 1);
        $each!(signature_2_2_2, u64, 2, 2, 2);
        $each!(signature_2_1_3, u64, 2, 1, 3);
        $each!(signature_2_0_4, u64, 2, 0, 4);
        $each!(signature_1_5_0, u64, 1, 5, 0);
        $each!(signature_1_4_1, u64, 1, 4, 1);
        $each!(signature_1_3_2, u64, 1, 3, 2);
        $each!(signature_1_2_3, u64, 1, 2, 3);
        $each!(signature_1_1_4, u64, 1, 1, 4);
        $each!(signature_1_0_5, u64, 1, 0, 5);
        $each!(signature_0_6_0, u64, 0, 6, 0);
        $each!(signature_0_5_1, u64, 0, 5, 1);
        $each!(signature_0_4_2, u64, 0, 4, 2);
        $each!(signature_0_3_3, u64, 0, 3, 3);
        $each!(signature_0_2_4, u64, 0, 2, 4);
        $each!(signature_0_1_5, u64, 0, 1, 5);
        $each!(signature_0_0_6, u64, 0, 0, 6);
    };
}
#[allow(unused_imports)]
pub(crate) use every_signature;
