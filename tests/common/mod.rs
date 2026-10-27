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
    let (a, e) = (actual.dense(), expected.dense());
    let scale = e.iter().fold(1.0_f64, |m, c| m.max(c.abs()));
    for blade in 0..64 {
        let diff = (a[blade] - e[blade]).abs();
        assert!(
            diff <= tol * scale,
            "blade {blade}: {} != {} (|diff| = {diff:e}, allowed {:e})\nactual:   {actual:?}\nexpected: {expected:?}",
            a[blade],
            e[blade],
            tol * scale,
        );
    }
}
