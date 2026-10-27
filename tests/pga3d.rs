//! PGA3D points, planes, lines, motors and flectors: camera projections, matrices both ways,
//! quaternions, the logarithm and sandwiches. The cases come from the 0.1 unit tests.

use std::f64::consts::TAU;

use nalgebra::{Matrix3x4, Matrix4x3, Vector4};
use seagrass_math::prelude::*;

/// Tolerance of the original `assert_approx_eq_vec3` helper.
const VEC3_TOLERANCE: f64 = 1e-15;
/// Tolerance of the original `assert_approx_eq_vec4` helper.
const VEC4_TOLERANCE: f64 = 1e-14;
/// Tolerance of the original `assert_approx_eq_mat4` helper (it applied the vec4 helper to
/// every row).
const MAT4_TOLERANCE: f64 = 1e-14;
/// Tolerance of the original `assert_approx_eq_mv` helper.
const MV_TOLERANCE: f64 = 1e-14;

fn assert_vec3_approx_eq(left: Vector3<f64>, right: Vector3<f64>) {
    let diff = left - right;
    assert!(
        diff.iter().all(|d| d.abs() < VEC3_TOLERANCE),
        "{left:?} not ~= {right:?}"
    );
}

fn assert_vec4_approx_eq(left: Vector4<f64>, right: Vector4<f64>) {
    let diff = left - right;
    assert!(
        diff.iter().all(|d| d.abs() < VEC4_TOLERANCE),
        "{left:?} not ~= {right:?}"
    );
}

fn assert_mat4_approx_eq(left: Matrix4<f64>, right: Matrix4<f64>) {
    let diff = left - right;
    assert!(
        diff.iter().all(|d| d.abs() < MAT4_TOLERANCE),
        "{left:?} not ~= {right:?}"
    );
}

/// Both multivectors must store exactly the same basis blades, and every coefficient must
/// agree within `MV_TOLERANCE`.
fn assert_mv_approx_eq(left: Multivector<f64, Pga3>, right: Multivector<f64, Pga3>) {
    let (left_blades, left_coefficients) = left.raw();
    let (right_blades, right_coefficients) = right.raw();
    assert_eq!(left_blades, right_blades, "{left:?} != {right:?}");
    for (l, r) in left_coefficients.iter().zip(&right_coefficients) {
        assert!((l - r).abs() < MV_TOLERANCE, "{left:?} != {right:?}");
    }
}

/// `rotation * translation`, where each part is `new_motor` of a line and an amount, as
/// the original camera and sandwich-product tests built it.
fn camera_motor(
    translation_dir: [f64; 3],
    distance: f64,
    rotation_axis: [f64; 3],
    angle: f64,
) -> Motor<f64, Pga3> {
    let [tx, ty, tz] = translation_dir;
    let translation_dir = Multivector::<_, Pga3>::new_ideal_line(tx, ty, tz);
    let camera_pos = Multivector::<_, Pga3>::new_motor(&translation_dir, distance);

    let [rx, ry, rz] = rotation_axis;
    let rotation_axis = Multivector::<_, Pga3>::new_line(rx, ry, rz);
    let camera_rot = Multivector::<_, Pga3>::new_motor(&rotation_axis, angle);

    camera_rot * camera_pos
}

/// Ported from `pga3d.rs` (`tests::camera`).
mod camera {
    use super::*;

    /// Original: `camera::simple_projection` and `camera::translation`.
    #[test]
    fn unrotated_camera_projects_scene_point_onto_image_plane() {
        // (translation direction, distance, rotation axis, angle, expected screen point)
        let cases = [
            // simple_projection
            (
                [1.0, 1.0, 1.0],
                0.0,
                [1.0, 0.0, 0.0],
                0.0,
                Vector3::new(0.0, 1.0, 0.0),
            ),
            // translation
            (
                [1.0, 1.0, 0.0],
                -3.0,
                [1.0, 0.0, 0.0],
                0.0,
                Vector3::new(-1.0, 0.0, 0.0),
            ),
        ];
        for (translation_dir, distance, rotation_axis, angle, expected) in cases {
            let image_plane = Multivector::<_, Pga3>::new_plane(0.0, 0.0, 1.0, 0.0);
            let focal_point = Multivector::<_, Pga3>::new_point(0.0, 0.0, -1.0);

            let camera_trans = camera_motor(translation_dir, distance, rotation_axis, angle);
            let camera_trans_inv = camera_trans.inverse();

            let scene_point = Multivector::<_, Pga3>::new_point(0.0, 3.0, 2.0);
            let scene_point_trans = camera_trans * scene_point * camera_trans_inv;

            let fp_to_scene = focal_point & scene_point_trans;
            let screen_point = fp_to_scene ^ image_plane;

            assert_vec3_approx_eq(screen_point.point_to_euclidean(), expected);
        }
    }

    /// Original: `camera::rotation` (both cases) and `camera::rotate_and_trans`.
    #[test]
    fn rotated_camera_projects_scene_point_onto_image_plane() {
        // (translation direction, distance, rotation axis, angle, simplify, expected)
        let cases = [
            // rotation, first case
            (
                [0.0, 0.0, 0.0],
                -3.0,
                [0.0, 0.0, 1.0],
                -TAU / 4.0,
                true,
                Vector3::new(1.0, 0.0, 0.0),
            ),
            // rotation, second case (the original does not simplify here)
            (
                [0.0, 0.0, 0.0],
                0.0,
                [0.0, 0.0, 1.0],
                TAU / 2.0,
                false,
                Vector3::new(0.0, -1.0, 0.0),
            ),
            // rotate_and_trans
            (
                [1.0, 1.0, 0.0],
                -3.0,
                [0.0, 0.0, 1.0],
                -TAU / 4.0,
                true,
                Vector3::new(0.0, 1.0, 0.0),
            ),
        ];
        for (translation_dir, distance, rotation_axis, angle, simplify, expected) in cases {
            let image_plane = Multivector::<_, Pga3>::new_plane(0.0, 0.0, 1.0, 0.0);
            let focal_point = Multivector::<_, Pga3>::new_point(0.0, 0.0, -1.0);

            let camera_trans = camera_motor(translation_dir, distance, rotation_axis, angle);
            let camera_trans_inv = camera_trans.inverse();

            let scene_point = Multivector::<_, Pga3>::new_point(0.0, 3.0, 2.0);
            let scene_point_trans = camera_trans * scene_point * camera_trans_inv;

            let fp_to_scene = &focal_point & scene_point_trans;
            let mut screen_point = fp_to_scene ^ &image_plane;
            if simplify {
                screen_point.simplify();
            }

            assert_vec3_approx_eq(screen_point.point_to_euclidean(), expected);
        }
    }

    /// Original: `camera::sanity_check`.
    #[test]
    fn identity_camera_projects_point_on_image_plane_to_itself() {
        let image_plane = Multivector::<_, Pga3>::new_plane(0.0, 0.0, 1.0, 0.0);
        let mut focal_point = Multivector::<_, Pga3>::new_point(0.0, 0.0, -1.0);
        focal_point.simplify();

        let camera_trans = camera_motor([0.0, 0.0, 0.0], 0.0, [0.0, 0.0, 0.0], 0.0);
        let camera_trans_inv = camera_trans.inverse();

        let scene_point = Multivector::<_, Pga3>::new_point(-0.08682, 0.4924, 0.00);
        let scene_point_trans = camera_trans * scene_point * camera_trans_inv;

        // The original also evaluates the by-reference join and discards the result.
        let _ = &focal_point & &scene_point_trans;

        let mut screen_point = (focal_point & scene_point_trans) ^ image_plane;
        screen_point.simplify();

        assert_vec3_approx_eq(
            screen_point.point_to_euclidean(),
            Vector3::new(-0.08682, 0.4924, 0.0),
        );
    }
}

/// Ported from `pga3d/multivector.rs`.
mod multivector {
    use super::*;

    /// Original: `create_plane_from_points`.
    #[test]
    fn join_of_three_points_is_their_plane() {
        // (three points, expected plane as (x, y, z, w))
        let cases = [
            (
                [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 0.0]],
                [0.0, 0.0, 1.0, 0.0],
            ),
            (
                [[1.0, 2.0, 3.0], [2.0, 1.0, 4.0], [2.0, 5.0, 7.0]],
                [-7.0, -3.0, 4.0, 1.0],
            ),
            // The original repeats the previous case verbatim.
            (
                [[1.0, 2.0, 3.0], [2.0, 1.0, 4.0], [2.0, 5.0, 7.0]],
                [-7.0, -3.0, 4.0, 1.0],
            ),
        ];
        for ([p1, p2, p3], [x, y, z, w]) in cases {
            let p1 = Multivector::<_, Pga3>::new_point(p1[0], p1[1], p1[2]);
            let p2 = Multivector::<_, Pga3>::new_point(p2[0], p2[1], p2[2]);
            let p3 = Multivector::<_, Pga3>::new_point(p3[0], p3[1], p3[2]);
            let c = p1 & p2 & p3;
            let out = Multivector::<_, Pga3>::new_plane(x, y, z, w);

            assert_eq!(c, out);
        }
    }

    /// Original: `project_point_onto_plane`.
    #[test]
    fn point_projects_onto_plane() {
        // (plane as (x, y, z, w), point, expected projection)
        let cases = [
            (
                [1.0, 0.0, 0.0, 0.0],
                [1.0, 1.0, 1.0],
                Vector3::new(0.0, 1.0, 1.0),
            ),
            (
                [0.0, 1.0, 0.0, 0.0],
                [1.0, 1.0, 1.0],
                Vector3::new(1.0, 0.0, 1.0),
            ),
            (
                [0.0, 0.0, 1.0, 0.0],
                [1.0, 1.0, 1.0],
                Vector3::new(1.0, 1.0, 0.0),
            ),
            (
                [-7.0, -3.0, 4.0, 1.0],
                [0.0, 0.0, 0.0],
                Vector3::new(
                    0.0945945945945946,
                    0.04054054054054054,
                    -0.05405405405405406,
                ),
            ),
            (
                [1.0, 0.0, 0.0, -1.0],
                [100.0, 0.0, 0.0],
                Vector3::new(1.0, 0.0, 0.0),
            ),
        ];
        for ([x, y, z, w], [px, py, pz], expected) in cases {
            let plane = Multivector::<_, Pga3>::new_plane(x, y, z, w);
            let p = Multivector::<_, Pga3>::new_point(px, py, pz);
            let proj = (&plane | p) * plane;

            assert_vec3_approx_eq(proj.point_to_euclidean(), expected);
        }
    }

    /// Original: first halves of `rotor` and `translator`.
    #[test]
    fn motor_times_its_inverse_is_scalar_one() {
        let cases = [
            // rotor
            Multivector::<_, Pga3>::new_motor(
                &Multivector::<_, Pga3>::new_line(5.0, 3.0, 2.0),
                TAU / 1.6,
            ),
            // translator
            Multivector::<_, Pga3>::new_motor(
                &Multivector::<_, Pga3>::new_ideal_line(2.0, 3.0, 5.0),
                TAU / 3.2,
            ),
        ];
        for r in cases {
            let r_inv = r.inverse();
            let mut approx_one = r * r_inv;
            approx_one.simplify();

            assert_eq!(
                Multivector::<_, Pga3>::from(approx_one),
                Multivector::<_, Pga3>::new_scalar(1.0).into(),
            );
        }
    }

    /// Original: second halves of `rotor` and `translator`.
    #[test]
    fn motor_conjugation_moves_point() {
        // (motor, point, expected image)
        let cases = [
            // rotor: half turn about z
            (
                Multivector::<_, Pga3>::new_motor(
                    &Multivector::<_, Pga3>::new_line(0.0, 0.0, 1.0),
                    TAU / 2.0,
                ),
                Multivector::<_, Pga3>::new_point(1.0, 0.0, 0.0),
                Vector3::new(-1.0, 0.0, 0.0),
            ),
            // translator: 4 along x
            (
                Multivector::<_, Pga3>::new_motor(
                    &Multivector::<_, Pga3>::new_ideal_line(1.0, 0.0, 0.0),
                    4.0,
                ),
                Multivector::<_, Pga3>::new_point(0.0, 0.0, 0.0),
                Vector3::new(4.0, 0.0, 0.0),
            ),
        ];
        for (r, p, expected) in cases {
            let r_inv = r.inverse();
            let p_moved = r * p * r_inv;

            assert_vec3_approx_eq(p_moved.point_to_euclidean(), expected);
        }
    }
}

/// Ported from `pga3d/multivector/matrix.rs`.
mod multivector_matrix {
    use super::*;

    /// Original: `point_to_vec`.
    #[test]
    fn point_converts_to_euclidean_coordinates() {
        let p1 = Multivector::<_, Pga3>::new_point(1.0, 2.0, 3.0);
        let v = p1.point_to_euclidean();

        assert_vec3_approx_eq(v, Vector3::new(1.0, 2.0, 3.0));
    }
}

/// Ported from `pga3d/motor/matrix.rs`.
mod motor_matrix {
    use super::*;

    /// Original: `translation`, first two blocks.
    #[test]
    fn translator_matrix_matches_sandwich_and_has_exact_inverse() {
        // (point, ideal line direction, amount)
        let cases = [
            ([1.0, 0.0, 0.0], [0.0, 1.0, 0.0], TAU / 4.0),
            ([5.0, -3.0, 7.7], [1.0, 3.0, 2.0], TAU / 3.2),
        ];
        for ([px, py, pz], [lx, ly, lz], amount) in cases {
            let p = Multivector::<_, Pga3>::new_point(px, py, pz);
            let l = Multivector::<_, Pga3>::new_ideal_line(lx, ly, lz);
            let r = Multivector::<_, Pga3>::new_motor(&l, amount);

            let p_trans = &r * &p * r.inverse();
            let (m, m_inv) = r.to_homogeneous_transformation_matrices();
            let p_homo_trans = m * p.point_to_homogeneous();

            assert_vec4_approx_eq(p_homo_trans, p_trans.point_to_homogeneous());
            assert_eq!(m * m_inv, Matrix4::identity());
            assert_eq!(m * m_inv, m_inv * m);
        }
    }

    /// Original: `translation`, third block.
    #[test]
    fn translator_inverse_matrix_matches_opposite_translator() {
        let p = Multivector::<_, Pga3>::new_point(5.0, -3.0, 7.7);

        let l = Multivector::<_, Pga3>::new_ideal_line(1.0, 3.0, 2.0);
        let r1 = Multivector::<_, Pga3>::new_motor(&l, TAU / 3.2);
        let r2 = Multivector::<_, Pga3>::new_motor(&l, -TAU / 3.2);

        let p_trans1 = &r1 * &p * r1.inverse();
        let p_trans2 = &r2 * &p * r2.inverse();

        let (m1, m1_inv) = r1.to_homogeneous_transformation_matrices();

        let p_homo_trans1 = m1 * p.point_to_homogeneous();
        let p_homo_trans2 = m1_inv * p.point_to_homogeneous();

        assert_vec4_approx_eq(p_homo_trans1, p_trans1.point_to_homogeneous());
        assert_vec4_approx_eq(p_homo_trans2, p_trans2.point_to_homogeneous());
    }

    /// Original: `rotation`, first two blocks.
    #[test]
    fn rotor_matrix_matches_sandwich_and_has_inverse() {
        // (point, line direction, angle)
        let cases = [
            ([1.0, 0.0, 0.0], [0.0, 1.0, 0.0], TAU / 4.0),
            ([5.0, -3.0, 7.7], [1.0, 3.0, 2.0], TAU / 3.2),
        ];
        for ([px, py, pz], [lx, ly, lz], angle) in cases {
            let p = Multivector::<_, Pga3>::new_point(px, py, pz);
            let l = Multivector::<_, Pga3>::new_line(lx, ly, lz);
            let r = Multivector::<_, Pga3>::new_motor(&l, angle);

            let p_trans = &r * &p * r.inverse();
            let (m, m_inv) = r.to_homogeneous_transformation_matrices();
            let p_homo_trans = m * p.point_to_homogeneous();

            assert_vec4_approx_eq(p_homo_trans, p_trans.point_to_homogeneous());
            assert_mat4_approx_eq(m * m_inv, Matrix4::identity());
            assert_mat4_approx_eq(m * m_inv, m_inv * m);
        }
    }

    /// Original: `rotation`, third block.
    #[test]
    fn rotor_inverse_matrix_matches_opposite_rotor() {
        let p = Multivector::<_, Pga3>::new_point(5.0, -3.0, 7.7);

        let l = Multivector::<_, Pga3>::new_line(1.0, 3.0, 2.0);
        let r1 = Multivector::<_, Pga3>::new_motor(&l, TAU / 3.2);
        let r2 = Multivector::<_, Pga3>::new_motor(&l, -TAU / 3.2);

        let p_trans1 = &r1 * &p * r1.inverse();
        let p_trans2 = &r2 * &p * r2.inverse();

        let (m1, m1_inv) = r1.to_homogeneous_transformation_matrices();

        let p_homo_trans1 = m1 * p.point_to_homogeneous();
        let p_homo_trans2 = m1_inv * p.point_to_homogeneous();

        assert_vec4_approx_eq(p_homo_trans1, p_trans1.point_to_homogeneous());
        assert_vec4_approx_eq(p_homo_trans2, p_trans2.point_to_homogeneous());
    }
}

/// Ported from `pga3d/motor/quaternion.rs`.
mod motor_quaternion {
    use super::*;

    /// Converts `Matrix4::from_euler_angles(roll, pitch, yaw)` to a motor and compares it
    /// with `rot_quat`.
    fn assert_euler_matrix_converts_to(
        roll: f64,
        pitch: f64,
        yaw: f64,
        rot_quat: Motor<f64, Pga3>,
    ) {
        let rot_mat = Matrix4::from_euler_angles(roll, pitch, yaw);
        let rot_mat_quat = Motor::<f64, Pga3>::from(rot_mat);

        assert_mv_approx_eq(rot_mat_quat.into(), rot_quat.into());
    }

    /// Original: `roll`, `pitch` and `yaw`.
    #[test]
    fn single_axis_euler_matrix_converts_to_rotation_motor() {
        // (roll, pitch, yaw)
        let cases = [
            // roll
            (TAU / 3.1, 0.0, 0.0),
            (-TAU / 3.1, 0.0, 0.0),
            // pitch
            (0.0, TAU / 4.9, 0.0),
            (0.0, -TAU / 4.9, 0.0),
            // yaw (the original checks the same angle twice)
            (0.0, 0.0, TAU / 2.8),
            (0.0, 0.0, TAU / 2.8),
        ];
        for (roll, pitch, yaw) in cases {
            let rot_quat = Motor::<_, Pga3>::new_rotation_motor(roll, pitch, yaw);

            assert_euler_matrix_converts_to(roll, pitch, yaw, rot_quat);
        }
    }

    /// Original: `roll_pitch`, and `pitch_yaw`, whose body repeats `roll_pitch` verbatim.
    #[test]
    fn roll_pitch_euler_matrix_converts_to_pitch_times_roll_motor() {
        // (roll, pitch)
        let cases = [
            // roll_pitch
            (TAU / 3.1, TAU / 4.9),
            (-TAU / 3.1, TAU / 4.9),
            (TAU / 3.1, -TAU / 4.9),
            (-TAU / 3.1, -TAU / 4.9),
            // pitch_yaw (identical to the rows above in the original)
            (TAU / 3.1, TAU / 4.9),
            (-TAU / 3.1, TAU / 4.9),
            (TAU / 3.1, -TAU / 4.9),
            (-TAU / 3.1, -TAU / 4.9),
        ];
        for (roll, pitch) in cases {
            let rot_quat = Motor::<_, Pga3>::new_rotation_motor(0.0, pitch, 0.0)
                * Motor::<_, Pga3>::new_rotation_motor(roll, 0.0, 0.0);

            assert_euler_matrix_converts_to(roll, pitch, 0.0, rot_quat);
        }
    }

    /// Original: `roll_yaw`.
    #[test]
    fn roll_yaw_euler_matrix_converts_to_yaw_times_roll_motor() {
        // (roll, yaw)
        let cases = [
            (TAU / 3.1, TAU / 2.8),
            (-TAU / 3.1, TAU / 2.8),
            (TAU / 3.1, -TAU / 2.8),
            (-TAU / 3.1, -TAU / 2.8),
        ];
        for (roll, yaw) in cases {
            let rot_quat = Motor::<_, Pga3>::new_rotation_motor(0.0, 0.0, yaw)
                * Motor::<_, Pga3>::new_rotation_motor(roll, 0.0, 0.0);

            assert_euler_matrix_converts_to(roll, 0.0, yaw, rot_quat);
        }
    }

    /// Original: `roll_pitch_yaw`.
    #[test]
    fn roll_pitch_yaw_euler_matrix_converts_to_positivized_composed_motor() {
        // (roll, pitch, yaw)
        let cases = [
            (TAU / 4.0, TAU / 4.0, TAU / 4.0),
            (-TAU / 3.1, TAU / 4.9, TAU / 2.8),
            (TAU / 3.1, -TAU / 4.9, TAU / 2.8),
            (-TAU / 3.1, -TAU / 4.9, TAU / 2.8),
            (TAU / 3.1, TAU / 4.9, -TAU / 2.8),
            (-TAU / 3.1, TAU / 4.9, -TAU / 2.8),
            (TAU / 3.1, -TAU / 4.9, -TAU / 2.8),
            (-TAU / 3.1, -TAU / 4.9, -TAU / 2.8),
        ];
        for (roll, pitch, yaw) in cases {
            let mut rot_quat = Motor::<_, Pga3>::new_rotation_motor(0.0, 0.0, yaw)
                * Motor::<_, Pga3>::new_rotation_motor(0.0, pitch, 0.0)
                * Motor::<_, Pga3>::new_rotation_motor(roll, 0.0, 0.0);
            rot_quat.quaternion_angle_positivize();

            assert_euler_matrix_converts_to(roll, pitch, yaw, rot_quat);
        }
    }
}

/// Ported from `pga3d/motor/from/from_matrix.rs`.
mod motor_from_matrix {
    use super::*;

    /// Original: `from_matrix3` (it converts a `Matrix4` despite its name).
    #[test]
    fn motor_from_rotation_matrix_moves_point_like_matrix() {
        // (roll, pitch, yaw)
        let cases = [
            (TAU / 4.0, 0.0, 0.0),
            (0.0, TAU / 4.0, 0.0),
            (0.0, 0.0, TAU / 4.0),
            (-TAU / 3.1, TAU / 4.9, -TAU / 2.8),
        ];
        for (roll, pitch, yaw) in cases {
            let rot_mat = Matrix4::from_euler_angles(roll, pitch, yaw);
            let rot_quat = Motor::<f64, Pga3>::from(rot_mat);

            let p = Multivector::<_, Pga3>::new_point(1.0, 2.0, 3.0);
            let p_vec = p.point_to_homogeneous();

            let p_mat = rot_mat * p_vec;
            let p_mv = rot_quat / p;

            assert_vec4_approx_eq(p_mv.point_to_homogeneous(), p_mat);
        }
    }

    /// Original: `from_matrix4`.
    #[test]
    fn motor_from_rigid_transform_matrix_moves_point_like_matrix() {
        let rot_mat = Matrix4::from_euler_angles(TAU / 3.1, TAU / 4.9, TAU / 2.8);
        let trans_mat = Matrix4::new_translation(&Vector3::new(1.2, 2.3, 3.4));
        let full_mat = trans_mat * rot_mat;
        let trans_mot = Motor::<f64, Pga3>::from(full_mat);

        let p = Multivector::<_, Pga3>::new_point(1.0, 2.0, 3.0);
        let p_vec = p.point_to_homogeneous();

        let p_mat = full_mat * p_vec;
        let p_mv = trans_mot / p;

        assert_vec4_approx_eq(p_mv.point_to_homogeneous(), p_mat);
    }
}

/// Ported from `pga3d/motor/linear_operations/sandwich_product.rs`.
mod motor_sandwich_product {
    use super::*;

    /// Original: `rotation`, `translation` and `combined`.
    #[test]
    fn motor_sandwich_operator_moves_point() {
        // (translation direction, rotation axis, expected image of (0, 3, 2)); every case
        // uses distance -3 and angle -TAU/4.
        let cases = [
            // rotation
            (
                [0.0, 0.0, 0.0],
                [0.0, 0.0, 1.0],
                Vector3::new(3.0, 0.0, 2.0),
            ),
            // translation
            (
                [0.0, 1.0, 1.0],
                [0.0, 0.0, 0.0],
                Vector3::new(0.0, 0.0, -1.0),
            ),
            // combined
            (
                [0.0, 1.0, 1.0],
                [0.0, 1.0, 0.0],
                Vector3::new(1.0, 0.0, 0.0),
            ),
        ];
        for (translation_dir, rotation_axis, expected) in cases {
            let camera_trans = camera_motor(translation_dir, -3.0, rotation_axis, -TAU / 4.0);

            let scene_point = Multivector::<_, Pga3>::new_point(0.0, 3.0, 2.0);
            let scene_point_trans = camera_trans / scene_point;

            assert_vec3_approx_eq(scene_point_trans.point_to_euclidean(), expected);
        }
    }
}

/// Ported from `pga3d/motor/non_linear_operations/logarithm.rs`.
mod motor_logarithm {
    use super::*;

    /// Original: `log`.
    #[test]
    fn log_of_motor_returns_its_generating_line() {
        let mut normalized_line = Multivector::<_, Pga3>::new_line(TAU / 8.0, TAU / 4.0, TAU / 3.2);
        normalized_line.normalize();

        let cases = [
            Multivector::<_, Pga3>::new_ideal_line(1.0, TAU / 4.0, 10.0),
            Multivector::<_, Pga3>::new_line(0.0, TAU / 4.0, 0.0),
            normalized_line,
        ];
        for l in cases {
            let r = Multivector::<_, Pga3>::new_motor(&l, 2.0);

            let mut lout = r.log();
            lout.simplify();

            assert_mv_approx_eq(l.into(), lout.into());
        }
    }
}

/// Ported from `pga3d/flector/matrix.rs`.
mod flector_matrix {
    use super::*;

    /// Original: `basic_plane`.
    #[test]
    fn plane_matrix_is_axis_reflection_matching_sandwich() {
        let p = Multivector::<_, Pga3>::new_point(1.0, 2.0, 3.0);

        let f = Multivector::<_, Pga3>::new_plane(0.0, 1.0, 0.0, 0.0);

        let (m, m_inv) = f.to_homogeneous_transformation_matrices();

        let p_mat_homo = m * p.point_to_homogeneous();

        let p_trans = &f / &p;
        let p_trans_homo = p_trans.point_to_homogeneous();

        assert_vec4_approx_eq(p_mat_homo, p_trans_homo);
        assert_eq!(
            m,
            Matrix4::new_nonuniform_scaling(&Vector3::new(1.0, -1.0, 1.0))
        );
        assert_eq!(m * m_inv, Matrix4::identity());
        assert_eq!(m * m_inv, m_inv * m);
    }

    /// Original: `different_plane`.
    #[test]
    fn normalized_plane_matrix_matches_sandwich_and_has_inverse() {
        let p = Multivector::<_, Pga3>::new_point(1.0, 2.0, 3.0);

        let mut f = Multivector::<_, Pga3>::new_plane(1.2, 3.2342, 5.534, 2.775);
        f.normalize();

        let (m, m_inv) = f.to_homogeneous_transformation_matrices();

        let p_mat_homo = m * p.point_to_homogeneous();

        let p_trans = &f / &p;
        let p_trans_homo = p_trans.point_to_homogeneous();

        assert_vec4_approx_eq(p_mat_homo, p_trans_homo);
        assert_mat4_approx_eq(m * m_inv, Matrix4::identity());
        assert_eq!(m * m_inv, m_inv * m);
    }

    /// Original: `basic_point` and `different_points`.
    #[test]
    fn point_matrix_reflects_through_point_like_sandwich() {
        // (reflecting point, expected image of (1, 2, 3) in homogeneous coordinates)
        let cases = [
            // basic_point
            ([0.0, 0.0, 0.0], Vector4::new(-1.0, -2.0, -3.0, 1.0)),
            // different_points
            ([1.2, 3.4, 7.5], Vector4::new(1.4, 4.8, 12.0, 1.0)),
            ([0.4, 1.1, 3.1], Vector4::new(-0.2, 0.2, 3.2, 1.0)),
        ];
        for ([fx, fy, fz], expected) in cases {
            let p = Multivector::<_, Pga3>::new_point(1.0, 2.0, 3.0);

            let f = Multivector::<_, Pga3>::new_point(fx, fy, fz);

            let (m, m_inv) = f.to_homogeneous_transformation_matrices();

            let p_mat_homo = m * p.point_to_homogeneous();

            let p_trans = &f / &p;
            let p_trans_homo = p_trans.point_to_homogeneous();

            assert_vec4_approx_eq(p_mat_homo, expected);
            assert_vec4_approx_eq(p_mat_homo, p_trans_homo);
            assert_mat4_approx_eq(m * m_inv, Matrix4::identity());
            assert_eq!(m * m_inv, m_inv * m);
        }
    }
}

/// Ported from `pga3d/flector/from/from_matrix.rs`.
mod flector_from_matrix {
    use super::*;

    /// Original: `from_matrix3`.
    #[test]
    fn flector_from_plane_matrix3_recovers_plane() {
        let f = Multivector::<_, Pga3>::new_plane(-1.0, 2.0, 3.0, 0.0);

        let (m, _) = f.to_homogeneous_transformation_matrices();

        let pm = Matrix3x4::<f64>::identity() * m * Matrix4x3::identity();
        let fout = Flector::<_, Pga3>::from(pm);

        assert_mv_approx_eq(f.into(), fout.into());
    }

    /// Original: `from_matrix4`.
    #[test]
    fn flector_from_matrix4_recovers_flector() {
        let cases = [
            Multivector::<_, Pga3>::new_point(-7.0, 5.0, 2.0),
            Multivector::<_, Pga3>::new_plane(-1.0, 2.0, 3.0, 1e-15),
            Multivector::<_, Pga3>::new_plane(-3.5, 12.3, -2.5, 1e-15),
        ];
        for p in cases {
            let (m, _) = p.to_homogeneous_transformation_matrices();

            let pfout = Flector::<_, Pga3>::from(m);

            assert_mv_approx_eq(p.into(), pfout.into());
        }
    }

    /// Original: `test_matrix_transform`.
    #[test]
    fn normalized_flector_matrix_moves_point_like_sandwich() {
        let cases = [
            // plane
            Multivector::<_, Pga3>::new_plane(-1.0, 2.0, 3.0, 0.0),
            // plane + point
            Multivector::<_, Pga3>::new_plane(-3.5, 12.3, -2.5, 1e-15)
                + Multivector::<_, Pga3>::new_point(-7.0, 5.0, 2.0),
        ];
        for mut p in cases {
            let w = Multivector::<_, Pga3>::new_point(1.0, 2.0, 3.0);
            p.normalize();

            let wm = w.point_to_homogeneous();
            let (m, _) = p.to_homogeneous_transformation_matrices();

            let out_mv = p / w;
            let out_mat = m * wm;

            assert_vec4_approx_eq(out_mat, out_mv.point_to_homogeneous());
        }
    }
}
