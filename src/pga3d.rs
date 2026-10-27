//! Projective geometric algebra of 3D space, signature (3, 0, 1): points, planes, lines and the
//! motors and flectors that move them.
//!
//! Planes are vectors `x e1 + y e2 + z e3 + w e4` (the plane `xX + yY + zZ + w = 0`), lines are
//! bivectors and points are trivectors. `&` joins and `^` meets.

#[cfg(feature = "nalgebra")]
use nalgebra::{Vector3, Vector4};

use crate::blade::{BladeIndex, BladeSet};
use crate::coefficient::{Coefficient, constant};
use crate::metric::Pga3;
use crate::multivector::Multivector;
use crate::plan::Coefficients;
use crate::typed::{Flector, Kind, Motor, Typed};

mod flector;
mod motor;

pub use motor::MotorCreationOrder;

/// e23, e13 and e12: a line through the origin.
const LINE_BLADES: u16 = 0x0068;
/// e14, e24 and e34: a line at infinity.
const IDEAL_LINE_BLADES: u16 = 0x1600;
/// 1, e12, e13 and e23: a rotation about the origin.
const ROTATION_BLADES: u16 = 0x0069;
/// e14, e24, e34 and e1234: a translation.
const TRANSLATION_BLADES: u16 = 0x9600;
/// Every even blade: a general motor.
const MOTOR_BLADES: u16 = 0x9669;

/// A PGA3D value with the given blades and coefficients, with zero coefficients dropped.
fn pga<K: Kind, T: Coefficient>(blades: u16, coefficients: &[T]) -> Typed<T, Pga3, K> {
    let mut value = Multivector::from_raw(
        Pga3::default(),
        blades,
        Coefficients::from_slice(coefficients),
    );
    value.remove_zeros();
    Typed::wrap(value)
}

impl<T: Coefficient> Multivector<T, Pga3> {
    /// The scalar `s`, as a motor.
    pub fn new_scalar(s: T) -> Motor<T, Pga3> {
        pga(<u16 as BladeSet>::SCALAR, &[s])
    }

    /// The point (x, y, z).
    ///
    /// # Examples
    ///
    /// ```
    /// use seagrass_math::prelude::*;
    ///
    /// let p = Multivector::<f64, Pga3>::new_point(1.0, 2.0, 3.0);
    /// assert_eq!(p.to_string(), "3,0,1>e123 - 3e124 + 2e134 - e234");
    /// ```
    pub fn new_point(x: T, y: T, z: T) -> Flector<T, Pga3> {
        Self::new_homogeneous_point(x, y, z, T::one())
    }

    /// The point (x/w, y/w, z/w); w = 0 gives a direction (a point at infinity).
    pub fn new_homogeneous_point(x: T, y: T, z: T, w: T) -> Flector<T, Pga3> {
        pga(<u16 as BladeSet>::GRADE[3], &[w, -z, y, -x])
    }

    /// The plane `x X + y Y + z Z + w = 0`.
    ///
    /// # Examples
    ///
    /// ```
    /// use seagrass_math::prelude::*;
    ///
    /// // The meet of the planes x = 1 and y = 2 and z = 3 is the point (1, 2, 3).
    /// let x = Multivector::<f64, Pga3>::new_plane(1.0, 0.0, 0.0, -1.0);
    /// let y = Multivector::<f64, Pga3>::new_plane(0.0, 1.0, 0.0, -2.0);
    /// let z = Multivector::<f64, Pga3>::new_plane(0.0, 0.0, 1.0, -3.0);
    /// assert_eq!((&x ^ &y ^ &z).point_to_euclidean(), Vector3::new(1.0, 2.0, 3.0));
    /// ```
    pub fn new_plane(x: T, y: T, z: T, w: T) -> Flector<T, Pga3> {
        pga(<u16 as BladeSet>::GRADE[1], &[x, y, z, w])
    }

    /// The line through the origin along (x, y, z), oriented so that a motor built from it
    /// rotates counterclockwise about that direction.
    ///
    /// # Examples
    ///
    /// ```
    /// use seagrass_math::prelude::*;
    ///
    /// // The z axis is the join of the origin and the point (0, 0, 1).
    /// let origin = Multivector::<f64, Pga3>::new_point(0.0, 0.0, 0.0);
    /// let up = Multivector::<f64, Pga3>::new_point(0.0, 0.0, 1.0);
    /// assert_eq!(&origin & &up, Multivector::<f64, Pga3>::new_line(0.0, 0.0, 1.0));
    /// ```
    pub fn new_line(x: T, y: T, z: T) -> Motor<T, Pga3> {
        pga(LINE_BLADES, &[-z, y, -x])
    }

    /// The line at infinity orthogonal to (x, y, z); a motor built from it translates along
    /// (x, y, z).
    pub fn new_ideal_line(x: T, y: T, z: T) -> Motor<T, Pga3> {
        pga(IDEAL_LINE_BLADES, &[x, y, z])
    }

    /// The motor `exp(value / 2 * line)`: a rotation by `value` times the line's weight about a
    /// line, or a translation by `value` times its length along an ideal line.
    ///
    /// # Examples
    ///
    /// ```
    /// use seagrass_math::prelude::*;
    ///
    /// // A quarter turn about the vertical line through (1, 0, 0).
    /// let axis = &Multivector::<f64, Pga3>::new_point(1.0, 0.0, 0.0)
    ///     & &Multivector::<f64, Pga3>::new_point(1.0, 0.0, 1.0);
    /// let turn = Multivector::<f64, Pga3>::new_motor(&axis, std::f64::consts::FRAC_PI_2);
    /// let p = &turn / Multivector::<f64, Pga3>::new_point(0.0, 0.0, 0.0);
    /// assert!((p.point_to_euclidean() - Vector3::new(1.0, -1.0, 0.0)).norm() < 1e-12);
    /// ```
    pub fn new_motor(line: &impl AsRef<Self>, value: T) -> Motor<T, Pga3> {
        Typed::wrap((line.as_ref() * value / constant(2.0)).exp())
    }

    /// The scalar 1.
    pub fn identity() -> Self {
        Self::from(T::one())
    }

    /// The stored blades and the coefficient of every blade, zero where not stored.
    pub fn raw(&self) -> (u16, [T; 16]) {
        let mut dense = [T::zero(); 16];
        for (index, c) in self.iter() {
            dense[index.bits() as usize] = *c;
        }
        (self.blades(), dense)
    }
}

#[cfg(feature = "nalgebra")]
impl<T: Coefficient> Multivector<T, Pga3> {
    /// The homogeneous coordinates (x, y, z, w) of a point.
    pub fn point_to_homogeneous(&self) -> Vector4<T> {
        Vector4::new(
            -self.get(BladeIndex::E234),
            self.get(BladeIndex::E134),
            -self.get(BladeIndex::E124),
            self.get(BladeIndex::E123),
        )
    }

    /// The Euclidean coordinates (x/w, y/w, z/w) of a point.
    pub fn point_to_euclidean(&self) -> Vector3<T> {
        let h = self.point_to_homogeneous();
        Vector3::new(h.x / h.w, h.y / h.w, h.z / h.w)
    }
}

#[cfg(feature = "nalgebra")]
impl<T: Coefficient, K: Kind> Typed<T, Pga3, K> {
    /// The homogeneous coordinates (x, y, z, w) of a point.
    pub fn point_to_homogeneous(&self) -> Vector4<T> {
        self.as_multivector().point_to_homogeneous()
    }

    /// The Euclidean coordinates (x/w, y/w, z/w) of a point.
    pub fn point_to_euclidean(&self) -> Vector3<T> {
        self.as_multivector().point_to_euclidean()
    }
}

impl<T: Coefficient, K: Kind> Typed<T, Pga3, K> {
    /// The stored blades and the coefficient of every blade, zero where not stored.
    pub fn raw(&self) -> (u16, [T; 16]) {
        self.as_multivector().raw()
    }
}
