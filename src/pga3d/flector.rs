//! PGA3D flectors: reflections in planes and point reflections.

#[cfg(feature = "nalgebra")]
use nalgebra::{Matrix3, Matrix3x4, Matrix4, Matrix4x3};

use super::*;

impl<T: Coefficient> Flector<T, Pga3> {
    /// The point (x, y, z); see [`Multivector::new_point`].
    pub fn from_point(x: T, y: T, z: T) -> Self {
        Multivector::new_point(x, y, z)
    }

    /// The point (x/w, y/w, z/w); see [`Multivector::new_homogeneous_point`].
    pub fn from_homogeneous_point(x: T, y: T, z: T, w: T) -> Self {
        Multivector::new_homogeneous_point(x, y, z, w)
    }

    /// The plane `x X + y Y + z Z + d = 0`; see [`Multivector::new_plane`]. As a versor it is
    /// the reflection in that plane.
    ///
    /// # Examples
    ///
    /// ```
    /// use seagrass_math::prelude::*;
    ///
    /// let mirror = Flector::<f64, Pga3>::from_plane(0.0, 0.0, 1.0, -1.0); // z = 1
    /// let p = &mirror / Multivector::<f64, Pga3>::new_point(2.0, 3.0, 0.0);
    /// assert!((p.point_to_euclidean() - Vector3::new(2.0, 3.0, 2.0)).norm() < 1e-12);
    /// ```
    pub fn from_plane(x: T, y: T, z: T, d: T) -> Self {
        Multivector::new_plane(x, y, z, d)
    }

    /// The plane `x X + y Y + z Z = 0` through the origin.
    pub fn from_origin_plane(x: T, y: T, z: T) -> Self {
        Multivector::new_plane(x, y, z, T::zero())
    }
}

#[cfg(feature = "nalgebra")]
impl<T: Coefficient> Flector<T, Pga3> {
    /// The homogeneous matrix of the flector and of its inverse, for a unit flector; the same
    /// A ± B construction as [`Motor::to_homogeneous_transformation_matrices`].
    pub fn to_homogeneous_transformation_matrices(&self) -> (Matrix4<T>, Matrix4<T>) {
        let hx = self.get(BladeIndex::E1);
        let hy = self.get(BladeIndex::E2);
        let hz = self.get(BladeIndex::E3);
        let hw = self.get(BladeIndex::E4);
        let sw = self.get(BladeIndex::E123);
        let sz = -self.get(BladeIndex::E124);
        let sy = self.get(BladeIndex::E134);
        let sx = -self.get(BladeIndex::E234);
        let (zero, one, two) = (T::zero(), T::one(), constant::<T>(2.0));

        #[rustfmt::skip] // Matrix layout.
        let a = Matrix4::new(
            two * (hy * hy + hz * hz) - one, -two * hx * hy, -two * hx * hz, two * (sx * sw - hx * hw),
            -two * hx * hy, two * (hx * hx + hz * hz) - one, -two * hy * hz, two * (sy * sw - hy * hw),
            -two * hx * hz, -two * hy * hz, two * (hx * hx + hy * hy) - one, two * (sz * sw - hz * hw),
            zero, zero, zero, one,
        );
        #[rustfmt::skip] // Matrix layout.
        let b = Matrix4::new(
            zero, two * hz * sw, -two * hy * sw, two * (hy * sz - hz * sy),
            -two * hz * sw, zero, two * hx * sw, two * (hz * sx - hx * sz),
            two * hy * sw, -two * hx * sw, zero, two * (hx * sy - hy * sx),
            zero, zero, zero, zero,
        );
        (a + b, a - b)
    }
}

#[cfg(feature = "nalgebra")]
impl<T: Coefficient> From<&Matrix3<T>> for Flector<T, Pga3> {
    /// The reflection in a plane through the origin with this 3x3 reflection matrix.
    fn from(plane_mat: &Matrix3<T>) -> Self {
        let (zero, one, two, four) = (T::zero(), T::one(), constant::<T>(2.0), constant::<T>(4.0));

        let sum_squares =
            ((plane_mat[(0, 0)] + plane_mat[(1, 1)] + plane_mat[(2, 2)] - one) / four) + one;
        let sum_two_x_squares = (plane_mat[(1, 1)] + plane_mat[(2, 2)] + two) / two;
        let sum_two_y_squares = (plane_mat[(0, 0)] + plane_mat[(2, 2)] + two) / two;
        let two_z_squares = (plane_mat[(0, 0)] + plane_mat[(1, 1)] + two) / two;

        let mut normal = [
            (sum_two_x_squares - sum_squares).sqrt(),
            (sum_two_y_squares - sum_squares).sqrt(),
            (two_z_squares - sum_squares).sqrt(),
        ];
        let x_y = -(plane_mat[(1, 0)] + plane_mat[(0, 1)]) / four;
        let x_z = -(plane_mat[(2, 0)] + plane_mat[(0, 2)]) / four;
        if x_y < zero {
            normal[1] = -normal[1];
        }
        if x_z < zero {
            normal[2] = -normal[2];
        }
        Self::from_origin_plane(-normal[0], -normal[1], -normal[2])
    }
}

#[cfg(feature = "nalgebra")]
impl<T: Coefficient> From<Matrix3<T>> for Flector<T, Pga3> {
    fn from(plane_mat: Matrix3<T>) -> Self {
        Self::from(&plane_mat)
    }
}

#[cfg(feature = "nalgebra")]
impl<T: Coefficient> From<&Matrix4<T>> for Flector<T, Pga3> {
    /// The flector with this homogeneous matrix: a plane reflection, or a plane reflection
    /// combined with a point reflection.
    fn from(flec_mat: &Matrix4<T>) -> Self {
        let (two, four) = (constant::<T>(2.0), constant::<T>(4.0));
        let plane_mat = Matrix3x4::identity() * flec_mat * Matrix4x3::identity();
        let plane = Self::from(&plane_mat);

        let hx = plane.get(BladeIndex::E1);
        let hy = plane.get(BladeIndex::E2);
        let hz = plane.get(BladeIndex::E3);
        let hw = plane.get(BladeIndex::E4);

        let sum = hx + hy + hz;
        let sw = if sum.is_zero() {
            T::one()
        } else {
            ((plane_mat[(0, 1)] - plane_mat[(1, 0)])
                + (plane_mat[(2, 0)] - plane_mat[(0, 2)])
                + (plane_mat[(1, 2)] - plane_mat[(2, 1)]))
                / (four * sum)
        };

        if sw.is_zero() {
            let hw =
                -(flec_mat[(0, 3)] + flec_mat[(1, 3)] + flec_mat[(2, 3)]) / (two * (hx + hy + hz));
            if hw.is_sign_negative() {
                return Self::from_plane(-hx, -hy, -hz, -hw);
            }
            return Self::from_plane(hx, hy, hz, hw);
        }

        let d1 = flec_mat[(0, 3)] / two + hx * hw;
        let d2 = flec_mat[(1, 3)] / two + hy * hw;
        let d3 = flec_mat[(2, 3)] / two + hz * hw;

        let p1 = Self::from_plane(sw, hz, hy, -d1);
        let p2 = Self::from_plane(hz, sw, hx, -d2);
        let p3 = Self::from_plane(hy, hx, sw, -d3);

        plane + (p1 ^ p2 ^ p3)
    }
}

#[cfg(feature = "nalgebra")]
impl<T: Coefficient> From<Matrix4<T>> for Flector<T, Pga3> {
    fn from(value: Matrix4<T>) -> Self {
        Self::from(&value)
    }
}
