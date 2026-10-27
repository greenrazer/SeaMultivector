//! PGA3D motors: rotations, translations and screw motions.

#[cfg(feature = "nalgebra")]
use nalgebra::{DualQuaternion, Matrix3, Matrix3x4, Matrix4, Matrix4x3, Quaternion};

use super::*;
use crate::traits::GeometricProduct;

/// Which part of a combined motor acts first, for
/// [`new_combined_motor`](Motor::new_combined_motor).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum MotorCreationOrder {
    /// Translate, then rotate about the origin.
    TranslationFirst,
    /// Rotate about the origin, then translate.
    RotationFirst,
}

impl<T: Coefficient> Motor<T, Pga3> {
    /// The motor that does nothing, the scalar 1.
    pub fn identity() -> Self {
        Self::from(T::one())
    }

    /// The translation by (x, y, z).
    ///
    /// # Examples
    ///
    /// ```
    /// use seagrass_math::prelude::*;
    ///
    /// let m = Motor::<f64, Pga3>::new_translation_motor(1.0, 2.0, 3.0);
    /// let p = &m / Multivector::<f64, Pga3>::new_point(0.0, 0.0, 0.0);
    /// assert_eq!(p.point_to_euclidean(), Vector3::new(1.0, 2.0, 3.0));
    /// ```
    pub fn new_translation_motor(x: T, y: T, z: T) -> Self {
        let line = Multivector::new_ideal_line(x, y, z);
        Multivector::new_motor(&line, T::one())
    }

    /// The rotation about the axis (rx, ry, rz) through the origin, by the axis length in radians,
    /// counterclockwise when looking down the axis.
    ///
    /// # Examples
    ///
    /// ```
    /// use seagrass_math::prelude::*;
    ///
    /// let quarter = Motor::<f64, Pga3>::new_rotation_motor(0.0, 0.0, std::f64::consts::FRAC_PI_2);
    /// let p = &quarter / Multivector::<f64, Pga3>::new_point(1.0, 0.0, 0.0);
    /// assert!((p.point_to_euclidean() - Vector3::new(0.0, 1.0, 0.0)).norm() < 1e-12);
    /// ```
    pub fn new_rotation_motor(rx: T, ry: T, rz: T) -> Self {
        let line = Multivector::new_line(rx, ry, rz);
        Multivector::new_motor(&line, T::one())
    }

    /// A rotation and a translation combined in the given order.
    ///
    /// # Examples
    ///
    /// ```
    /// use seagrass_math::prelude::*;
    ///
    /// use std::f64::consts::FRAC_PI_2;
    ///
    /// let first_rotate = Motor::<f64, Pga3>::new_combined_motor(
    ///     1.0, 0.0, 0.0, 0.0, 0.0, FRAC_PI_2, MotorCreationOrder::RotationFirst,
    /// );
    /// let p = &first_rotate / Multivector::<f64, Pga3>::new_point(1.0, 0.0, 0.0);
    /// assert!((p.point_to_euclidean() - Vector3::new(1.0, 1.0, 0.0)).norm() < 1e-12);
    /// ```
    #[allow(clippy::too_many_arguments)] // Three translation and three rotation components.
    pub fn new_combined_motor(
        x: T,
        y: T,
        z: T,
        rx: T,
        ry: T,
        rz: T,
        order: MotorCreationOrder,
    ) -> Self {
        let translation = Self::new_translation_motor(x, y, z);
        let rotation = Self::new_rotation_motor(rx, ry, rz);
        match order {
            MotorCreationOrder::TranslationFirst => rotation * translation,
            MotorCreationOrder::RotationFirst => translation * rotation,
        }
    }

    /// The rotation with quaternion coefficients `qw + qx i + qy j + qz k`, stored on
    /// `qw + qz e12 + qy e13 + qx e23`.
    pub fn from_quaternion(qx: T, qy: T, qz: T, qw: T) -> Self {
        pga(ROTATION_BLADES, &[qw, qz, qy, qx])
    }

    /// The translation part `dx e14 + dy e24 + dz e34 + dw e1234`.
    pub fn from_dual_quaternion(dx: T, dy: T, dz: T, dw: T) -> Self {
        pga(TRANSLATION_BLADES, &[dx, dy, dz, dw])
    }

    /// A motor from its eight coefficients, `qw + qx e12 + qy e13 + qz e23` plus
    /// `dx e14 + dy e24 + dz e34 + dw e1234`.
    #[allow(clippy::too_many_arguments)] // The eight coefficients of a motor.
    pub fn from_combined_quaternion(
        dx: T,
        dy: T,
        dz: T,
        dw: T,
        qx: T,
        qy: T,
        qz: T,
        qw: T,
    ) -> Self {
        pga(MOTOR_BLADES, &[qw, qx, qy, qz, dx, dy, dz, dw])
    }

    /// Negates the motor if its scalar part is negative. M and -M are the same motion.
    pub fn quaternion_angle_positivize(&mut self) {
        if self.scalar().is_sign_negative() {
            self.negate();
        }
    }

    /// A copy with a non-negative scalar part; see
    /// [`quaternion_angle_positivize`](Self::quaternion_angle_positivize).
    #[must_use]
    pub fn quaternion_angle_positivized(&self) -> Self {
        let mut out = self.clone();
        out.quaternion_angle_positivize();
        out
    }

    /// The logarithm: the bivector B with `exp(B) = self`, for a unit motor. Half of the
    /// bivector of the motion, so `(m.log() * 0.5).exp()` moves halfway along the same screw.
    ///
    /// Closed form (De Keninck): writing the motor as `s + B + p I`, with `Bₑ` the part of B
    /// without the null vector, the logarithm is `b B + c Bₑ I`, where
    /// `b = atan2(|Bₑ|, s) / |Bₑ|` and `c = p (1 - s b) / |Bₑ|²`. Nothing divides by s, so half
    /// turns (s = 0) are exact; a pure translation (`Bₑ = 0`) gives `B / s`. Scaling the motor
    /// does not change the result.
    ///
    /// # Examples
    ///
    /// ```
    /// use seagrass_math::prelude::*;
    ///
    /// let screw = Motor::<f64, Pga3>::new_translation_motor(0.0, 0.0, 4.0)
    ///     * Motor::new_rotation_motor(0.0, 0.0, std::f64::consts::PI);
    /// let halfway = (screw.log() * 0.5).exp(); // a quarter turn and half the translation
    /// let p = &halfway / Multivector::<f64, Pga3>::new_point(1.0, 0.0, 0.0);
    /// assert!((p.point_to_euclidean() - Vector3::new(0.0, 1.0, 2.0)).norm() < 1e-12);
    /// ```
    pub fn log(&self) -> Self {
        let s = self.scalar();
        let bivector = self.grade(2);
        let euclidean_mask = <u16 as BladeSet>::GRADE[2] & <u16 as BladeSet>::EXCLUDES[4];
        let euclidean = self.as_multivector().masked(euclidean_mask);
        let norm_squared = euclidean
            .coefficients()
            .iter()
            .fold(T::zero(), |sum, c| sum + *c * *c);

        if norm_squared.is_zero() {
            return Typed::wrap(bivector / s);
        }

        let norm = norm_squared.sqrt();
        let b = norm.atan2(s) / norm;
        let c = self.pseudoscalar() * (T::one() - s * b) / norm_squared;
        let pseudoscalar = Multivector::basis(Pga3::default(), BladeIndex::E1234);
        Typed::wrap(bivector * b - euclidean.geometric_product(&pseudoscalar) * c)
    }
}

#[cfg(feature = "nalgebra")]
impl<T: Coefficient> Motor<T, Pga3> {
    /// The rotation quaternion `(rw, rx, ry, rz)` and translation part `(tx, ty, tz, tw)` in the
    /// sign convention of [`to_homogeneous_transformation_matrices`](Self::to_homogeneous_transformation_matrices).
    fn quaternion_parts(&self) -> ([T; 4], [T; 4]) {
        let coefficients = self.as_multivector();
        let get = |index| coefficients.get(index);
        let rotation = [
            get(BladeIndex::SCALAR),
            -get(BladeIndex::E23),
            get(BladeIndex::E13),
            -get(BladeIndex::E12),
        ];
        let translation = [
            get(BladeIndex::E14),
            get(BladeIndex::E24),
            get(BladeIndex::E34),
            get(BladeIndex::E1234),
        ];
        (rotation, translation)
    }

    /// The motor as a nalgebra dual quaternion describing the same rigid motion: the real part
    /// is the rotation quaternion and the dual part is half the translation times it.
    ///
    /// # Examples
    ///
    /// ```
    /// use seagrass_math::prelude::*;
    ///
    /// use nalgebra::UnitDualQuaternion;
    ///
    /// let m = Matrix4::new_translation(&Vector3::new(1.0, 2.0, 3.0)) * Matrix4::from_euler_angles(0.3, -0.2, 1.1);
    /// let dq = UnitDualQuaternion::new_unchecked(Motor::<f64, Pga3>::from(m).to_dual_quaternion());
    /// assert!((dq.to_homogeneous() - m).norm() < 1e-12);
    /// ```
    pub fn to_dual_quaternion(&self) -> DualQuaternion<T> {
        let ([rw, rx, ry, rz], [tx, ty, tz, tw]) = self.quaternion_parts();
        DualQuaternion::from_real_and_dual(
            Quaternion::new(rw, rx, ry, rz),
            Quaternion::new(tw, tx, ty, tz),
        )
    }

    /// The homogeneous matrix of the motor and of its inverse, for a unit motor.
    ///
    /// Reversing a motor negates its bivector part and keeps its scalar and e1234 parts. Each
    /// matrix entry is a sum of terms that keep their sign under reversal (constants and products
    /// of two bivector coefficients), collected in A, and terms that flip it (a scalar or e1234
    /// coefficient times a bivector coefficient), collected in B. So `M = A + B`, and since a
    /// unit motor's inverse is its reverse, `M⁻¹ = A - B`: no matrix inversion runs.
    ///
    /// # Examples
    ///
    /// ```
    /// use seagrass_math::prelude::*;
    ///
    /// let m = Matrix4::new_translation(&Vector3::new(1.0, 2.0, 3.0)) * Matrix4::from_euler_angles(0.3, -0.2, 1.1);
    /// let (forward, inverse) = Motor::<f64, Pga3>::from(m).to_homogeneous_transformation_matrices();
    /// assert!((forward - m).norm() < 1e-12);
    /// assert!((forward * inverse - Matrix4::identity()).norm() < 1e-12);
    /// ```
    pub fn to_homogeneous_transformation_matrices(&self) -> (Matrix4<T>, Matrix4<T>) {
        let ([rw, rx, ry, rz], [tx, ty, tz, tw]) = self.quaternion_parts();
        let (zero, one, two) = (T::zero(), T::one(), constant::<T>(2.0));

        #[rustfmt::skip] // Matrix layout.
        let a = Matrix4::new(
            one - two * (ry * ry + rz * rz), two * rx * ry, two * rx * rz, two * (ry * tz - rz * ty),
            two * rx * ry, one - two * (rx * rx + rz * rz), two * ry * rz, two * (rz * tx - rx * tz),
            two * rx * rz, two * ry * rz, one - two * (rx * rx + ry * ry), two * (rx * ty - ry * tx),
            zero, zero, zero, one,
        );
        #[rustfmt::skip] // Matrix layout.
        let b = Matrix4::new(
            zero, -two * rz * rw, two * ry * rw, two * (rw * tx - rx * tw),
            two * rz * rw, zero, -two * rx * rw, two * (rw * ty - ry * tw),
            -two * ry * rw, two * rx * rw, zero, two * (rw * tz - rz * tw),
            zero, zero, zero, zero,
        );
        (a + b, a - b)
    }
}

#[cfg(feature = "nalgebra")]
impl<T: Coefficient> From<&Matrix3<T>> for Motor<T, Pga3> {
    /// The rotation with this rotation matrix, with a non-negative scalar part.
    fn from(rot: &Matrix3<T>) -> Self {
        let (one, two, four) = (T::one(), constant::<T>(2.0), constant::<T>(4.0));
        let trace = rot.trace();
        let mut quat = [T::zero(); 4];
        if trace > T::zero() {
            let s = one / (two * (trace + one).sqrt());
            quat[0] = (rot[(2, 1)] - rot[(1, 2)]) * s;
            quat[1] = (rot[(0, 2)] - rot[(2, 0)]) * s;
            quat[2] = (rot[(1, 0)] - rot[(0, 1)]) * s;
            quat[3] = one / (four * s);
        } else if rot[(0, 0)] > rot[(1, 1)] && rot[(0, 0)] > rot[(2, 2)] {
            let s = two * (one + rot[(0, 0)] - rot[(1, 1)] - rot[(2, 2)]).sqrt();
            quat[0] = s / four;
            quat[1] = (rot[(0, 1)] + rot[(1, 0)]) / s;
            quat[2] = (rot[(0, 2)] + rot[(2, 0)]) / s;
            quat[3] = (rot[(2, 1)] - rot[(1, 2)]) / s;
        } else if rot[(1, 1)] > rot[(2, 2)] {
            let s = two * (one + rot[(1, 1)] - rot[(0, 0)] - rot[(2, 2)]).sqrt();
            quat[0] = (rot[(0, 1)] + rot[(1, 0)]) / s;
            quat[1] = s / four;
            quat[2] = (rot[(1, 2)] + rot[(2, 1)]) / s;
            quat[3] = (rot[(0, 2)] - rot[(2, 0)]) / s;
        } else {
            let s = two * (one + rot[(2, 2)] - rot[(0, 0)] - rot[(1, 1)]).sqrt();
            quat[0] = (rot[(0, 2)] + rot[(2, 0)]) / s;
            quat[1] = (rot[(1, 2)] + rot[(2, 1)]) / s;
            quat[2] = s / four;
            quat[3] = (rot[(1, 0)] - rot[(0, 1)]) / s;
        }

        let mut motor = Self::from_quaternion(-quat[0], quat[1], -quat[2], quat[3]);
        motor.quaternion_angle_positivize();
        motor
    }
}

#[cfg(feature = "nalgebra")]
impl<T: Coefficient> From<Matrix3<T>> for Motor<T, Pga3> {
    fn from(rot: Matrix3<T>) -> Self {
        Self::from(&rot)
    }
}

#[cfg(feature = "nalgebra")]
impl<T: Coefficient> From<&Matrix4<T>> for Motor<T, Pga3> {
    /// The rigid motion with this homogeneous matrix: its rotation part, then its translation.
    fn from(value: &Matrix4<T>) -> Self {
        let rotation = Matrix3x4::identity() * (value * Matrix4x3::identity());
        let translation = Self::new_translation_motor(value[(0, 3)], value[(1, 3)], value[(2, 3)]);
        translation * Self::from(rotation)
    }
}

#[cfg(feature = "nalgebra")]
impl<T: Coefficient> From<Matrix4<T>> for Motor<T, Pga3> {
    fn from(value: Matrix4<T>) -> Self {
        Self::from(&value)
    }
}
