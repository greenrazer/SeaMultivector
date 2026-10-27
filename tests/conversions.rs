//! Conversions between PGA3D motors and nalgebra types describe the same rigid motion.

mod common;

use common::{TOL_F64, assert_close};
use nalgebra::{Isometry3, Translation3, UnitDualQuaternion, UnitQuaternion};
use seagrass_math::prelude::*;

/// A few rigid motions as homogeneous matrices: rotation, translation, and both.
fn motions() -> Vec<Matrix4<f64>> {
    vec![
        Matrix4::identity(),
        Matrix4::from_euler_angles(0.3, -0.2, 1.1),
        Matrix4::new_translation(&Vector3::new(1.0, -2.0, 3.0)),
        Matrix4::new_translation(&Vector3::new(1.0, 2.0, 3.0))
            * Matrix4::from_euler_angles(0.3, -0.2, 1.1),
        Matrix4::new_translation(&Vector3::new(-4.0, 0.5, 2.0))
            * Matrix4::from_euler_angles(2.9, 1.2, -2.5),
    ]
}

#[test]
fn dual_quaternion_describes_the_same_motion_as_the_matrix() {
    for m in motions() {
        let motor = Motor::<f64, Pga3>::from(m);
        let dq = UnitDualQuaternion::new_unchecked(motor.to_dual_quaternion());
        assert!((dq.to_homogeneous() - m).norm() < 1e-12, "{m}");
    }
}

#[test]
fn dual_quaternion_round_trips_through_nalgebra_isometries() {
    let isometry = Isometry3::from_parts(
        Translation3::new(1.0, 2.0, 3.0),
        UnitQuaternion::from_euler_angles(0.3, -0.2, 1.1),
    );
    let motor = Motor::<f64, Pga3>::from(isometry.to_homogeneous());
    let dq = UnitDualQuaternion::new_unchecked(motor.to_dual_quaternion());
    let expected = UnitDualQuaternion::from_isometry(&isometry);
    // q and -q are the same motion.
    let distance = (dq.as_ref() - expected.as_ref())
        .real
        .norm()
        .min((dq.as_ref() + expected.as_ref()).real.norm());
    assert!(distance < 1e-12);
    assert!((dq.translation().vector - isometry.translation.vector).norm() < 1e-12);
}

#[test]
fn matrix_round_trips_through_a_motor() {
    for m in motions() {
        let motor = Motor::<f64, Pga3>::from(m);
        let (forward, inverse) = motor.to_homogeneous_transformation_matrices();
        assert!((forward - m).norm() < 1e-12, "{m}");
        assert!(
            (forward * inverse - Matrix4::identity()).norm() < 1e-12,
            "{m}"
        );
        let again = Motor::<f64, Pga3>::from(forward);
        assert_close(&again, &motor.quaternion_angle_positivized(), TOL_F64);
    }
}

#[test]
fn applying_the_motor_matches_multiplying_by_its_matrix() {
    for m in motions() {
        let motor = Motor::<f64, Pga3>::from(m);
        for p in [
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, -2.0, 0.5),
            Vector3::new(-3.0, 4.0, 7.0),
        ] {
            let moved = &motor / Multivector::<f64, Pga3>::new_point(p.x, p.y, p.z);
            let expected = (m * p.push(1.0)).xyz();
            assert!(
                (moved.point_to_euclidean() - expected).norm() < 1e-12,
                "{m} {p}"
            );
        }
    }
}
