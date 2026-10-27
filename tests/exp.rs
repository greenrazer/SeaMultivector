//! Regression tests for `exp` on arguments whose square is not a scalar, which go through the
//! scaled Taylor series rather than a closed form.

mod common;

use common::{TOL_F64, assert_close};
use seagrass_math::prelude::*;

fn euclid3(s: &str) -> Multivector<f64, Dyn3> {
    Multivector::<_, Dyn3>::from_str(&format!("3>{s}")).unwrap()
}

/// e1 squares to 1 and commutes with scalars, so exp(x + x e1) = e^x (cosh x + sinh x e1).
/// The square (x + x e1)^2 = 2x^2 + 2x^2 e1 is not a scalar, so this exercises the series.
#[test]
fn exp_of_scalar_plus_vector_matches_closed_form() {
    for x in [0.25_f64, 1.0, 4.0, 10.0, 25.0] {
        let a = euclid3(&format!("{x} + {x}e1"));
        let expected = euclid3(&format!(
            "{} + {}e1",
            x.exp() * x.cosh(),
            x.exp() * x.sinh()
        ));
        assert_close(&a.exp(), &expected, TOL_F64);
    }
}

/// Rotation about the z axis and translation along it commute, so the exponential of their sum
/// equals the composed motors. The screw bivector squares to a scalar plus an e1234 part, so this
/// also goes through the series.
#[test]
fn exp_of_screw_bivector_matches_composed_motors() {
    for (angle, distance) in [(0.5_f64, 1.0_f64), (3.0, 10.0), (-2.5, 40.0)] {
        let rotation = Multivector::<_, Pga3>::new_line(0.0, 0.0, angle);
        let translation = Multivector::<_, Pga3>::new_ideal_line(0.0, 0.0, distance);
        let screw = ((&rotation + &translation) * 0.5).exp();

        let composed = Motor::<_, Pga3>::new_translation_motor(0.0, 0.0, distance)
            * Motor::<_, Pga3>::new_rotation_motor(0.0, 0.0, angle);
        assert_close(&screw, &composed, TOL_F64);
    }
}

/// exp(A) exp(-A) = 1 for a general mixed-grade argument, as long as A commutes with -A.
#[test]
fn exp_of_mixed_grade_argument_has_inverse_exp_of_negation() {
    let a =
        Multivector::<f64, Dyn4>::from_str("3,0,1>1 + e1 - 2e23 + 0.5e4 + 3e124 - e1234").unwrap();
    let product = &a.exp() * &(-&a).exp();
    assert_close(
        &product,
        &Multivector::<f64, Dyn4>::from_str("3,0,1>1").unwrap(),
        TOL_F64,
    );
}
