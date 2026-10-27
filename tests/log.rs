//! The PGA3D motor logarithm inverts the exponential.

mod common;

use std::f64::consts::PI;

use common::{TOL_F64, assert_close};
use seagrass_math::prelude::*;

type Pga = Multivector<f64, Pga3>;

/// Unit motors: rotations by various angles (including half turns), translations and screws.
fn motors() -> Vec<Motor<f64, Pga3>> {
    let mut out = Vec::new();
    for angle in [0.3, 1.0, 2.5, PI - 1e-3, PI] {
        for axis in [(1.0, 0.0, 0.0), (0.0, 0.6, 0.8), (-0.48, 0.6, 0.64)] {
            let rotation = Motor::<f64, Pga3>::new_rotation_motor(
                axis.0 * angle,
                axis.1 * angle,
                axis.2 * angle,
            );
            out.push(rotation.clone());
            out.push(Motor::new_translation_motor(1.0, -2.0, 0.5) * rotation.clone());
            out.push(rotation * Motor::new_translation_motor(-3.0, 0.25, 2.0));
        }
    }
    out.push(Motor::new_translation_motor(1.0, 2.0, 3.0));
    out
}

#[test]
fn exp_of_log_gives_back_the_motor() {
    for motor in motors() {
        let back = motor.log().exp();
        // M and -M are the same motion, and log picks the representative with a non-negative
        // scalar part.
        let target = motor.quaternion_angle_positivized();
        assert_close(&back.quaternion_angle_positivized(), &target, TOL_F64);
    }
}

#[test]
fn log_is_a_bivector() {
    for motor in motors() {
        let log = motor.log();
        for k in [0, 3, 4] {
            assert!(
                log.grade(k).iter().all(|(_, c)| c.abs() < 1e-12),
                "grade {k} of log({motor}) = {log}"
            );
        }
    }
}

#[test]
fn half_turn_log_moves_halfway() {
    // A half turn about the z axis through (1, 0, 0): the halfway motor is a quarter turn.
    let shift = Motor::<f64, Pga3>::new_translation_motor(1.0, 0.0, 0.0);
    let half_turn = &shift * &Motor::new_rotation_motor(0.0, 0.0, PI) * shift.reversed();
    let quarter = (half_turn.log() * 0.5).exp();
    let p = &quarter / Pga::new_point(0.0, 0.0, 0.0);
    assert!(
        (p.point_to_euclidean() - Vector3::new(1.0, -1.0, 0.0)).norm() < 1e-12,
        "{p}"
    );
}
