use std::ops::{Deref, DerefMut, Shl, Shr, BitOr, BitXor, Mul, Not, BitAnd, Index, AddAssign, Add, Neg, Sub, Div, SubAssign, MulAssign, DivAssign};

use nalgebra::{Vector3, Vector4, Matrix4};

use crate::coefficient::SeaCoefficient;
use crate::multivector::metric::SeaMetric;
use crate::multivector::multivector_traits::{Multivector, MultivectorMut, Versor, Motor};
use crate::multivector::multivector_product_traits::{
	LeftContraction, 
	RightContraction, 
	InnerProduct, 
	OuterProduct, 
	GeometricProduct, 
	RegressiveProduct, 
	ScalarProduct,
	SandwichProduct
};
use crate::multivector::multivector_unary_traits::{
	MultivectorConjugation, 
	MultivectorConjugationMut, 
	MultivectorDual, 
	MultivectorDualMut, 
	MultivectorNorm,
	MultivectorNormMut
};
use crate::multivector::grade::SeaGrade;
use crate::multivector::basis_blade_index::SeaBasisBladeIndexLow;
use crate::multivector::settings::MotorCreationOrder;

use super::basis_blade_set::SeaBasisBladeSet;
use super::core::SeaLowDimensionalSparseMultivector;
use super::SeaMultivector4D;

mod multivector;
mod motor;

pub use multivector::SeaMultivectorPGA3D;
pub use motor::SeaMotorPGA3D;

type SeaCoreMultivector4D<T> = SeaLowDimensionalSparseMultivector::<T, u16>;
const METRIC: SeaMetric = SeaMetric{positive:3, negitive:0, zero:1, dimensions:4};

const LINE_BASIS_BLADES: u16 = 0x0068;
const IDEAL_LINE_BASIS_BLADES: u16 = 0x1600;


#[cfg(test)]
mod tests {
	use std::f64::consts::TAU;

	use super::*;

	pub fn assert_approx_eq_vec4<T: SeaCoefficient + 'static>(a:Vector4<T>, b: Vector4<T>) {
        let epsilon = T::from(1e-14).unwrap();
		let q = a - b;
		assert!(q.x.abs() < epsilon, "{:?} not ~= {:?}", a, b);
		assert!(q.y.abs() < epsilon, "{:?} not ~= {:?}", a, b);
		assert!(q.z.abs() < epsilon, "{:?} not ~= {:?}", a, b);
		assert!(q.w.abs() < epsilon, "{:?} not ~= {:?}", a, b);
    }

	pub fn assert_approx_eq_vec3<T: SeaCoefficient + 'static>(a:Vector3<T>, b: Vector3<T>) {
        let epsilon = T::from(1e-15).unwrap();
		let q = a - b;
		assert!(q.x.abs() < epsilon, "{:?} not ~= {:?}", a, b);
		assert!(q.y.abs() < epsilon, "{:?} not ~= {:?}", a, b);
		assert!(q.z.abs() < epsilon, "{:?} not ~= {:?}", a, b);
    }

	pub fn assert_approx_eq_mat4<T: SeaCoefficient + 'static>(a:Matrix4<T>, b: Matrix4<T>) {
		assert_approx_eq_vec4(a.row(0).transpose(),b.row(0).transpose());
		assert_approx_eq_vec4(a.row(1).transpose(),b.row(1).transpose());
		assert_approx_eq_vec4(a.row(2).transpose(),b.row(2).transpose());
		assert_approx_eq_vec4(a.row(3).transpose(),b.row(3).transpose());
    }

	mod camera {
		use super::*;

		#[test]
		fn simple_projection() {
			let image_plane = SeaMultivectorPGA3D::new_plane(0.0,0.0,1.0, 0.0);
			let focal_point = SeaMultivectorPGA3D::new_point(0.0, 0.0, -1.0);
	
			let translation_dir = SeaMultivectorPGA3D::new_ideal_line(1.0,1.0,1.0);
			let camera_pos = SeaMultivectorPGA3D::new_motor(&translation_dir, 0.0);
	
			let rotation_axis = SeaMultivectorPGA3D::new_line(1.0,0.0,0.0);
			let camera_rot = SeaMultivectorPGA3D::new_motor(&rotation_axis, 0.0);
	
			let camera_trans = camera_rot * camera_pos;
			let camera_trans_inv = camera_trans.inverse();
			
			let scene_point = SeaMultivectorPGA3D::new_point(0.0,3.0,2.0);
			let scene_point_trans = camera_trans * scene_point * camera_trans_inv;
	
			let fp_to_scene = focal_point&scene_point_trans;
	
			let screen_point = fp_to_scene^image_plane;
	
			assert_approx_eq_vec3(
				screen_point.point_to_euclidean(),
				Vector3::new(0.0, 1.0, 0.0)
			);
		}

		#[test]
		fn translation() {
			let image_plane = SeaMultivectorPGA3D::new_plane(0.0,0.0,1.0, 0.0);
			let focal_point = SeaMultivectorPGA3D::new_point(0.0, 0.0, -1.0);
	
			let translation_dir = SeaMultivectorPGA3D::new_ideal_line(1.0,1.0,0.0);
			let camera_pos = SeaMultivectorPGA3D::new_motor(&translation_dir, -3.0);
	
			let rotation_axis = SeaMultivectorPGA3D::new_line(1.0,0.0,0.0);
			let camera_rot = SeaMultivectorPGA3D::new_motor(&rotation_axis, 0.0);
	
			let camera_trans = camera_rot * camera_pos;
			let camera_trans_inv = camera_trans.inverse();
			
			let scene_point = SeaMultivectorPGA3D::new_point(0.0,3.0,2.0);
			let scene_point_trans = camera_trans * scene_point * camera_trans_inv;
	
			let fp_to_scene = focal_point&scene_point_trans;
	
			let screen_point = fp_to_scene^image_plane;
	
			assert_approx_eq_vec3(
				screen_point.point_to_euclidean(),
				Vector3::new(-1.0, 0.0, 0.0)
			);
		}

		
		#[test]
		fn rotation() {
			let image_plane = SeaMultivectorPGA3D::new_plane(0.0,0.0,1.0, 0.0);
			let focal_point = SeaMultivectorPGA3D::new_point(0.0, 0.0, -1.0);
	
			let translation_dir = SeaMultivectorPGA3D::new_ideal_line(0.0,0.0,0.0);
			let camera_pos = SeaMultivectorPGA3D::new_motor(&translation_dir, -3.0);
	
			let rotation_axis = SeaMultivectorPGA3D::new_line(0.0,0.0,1.0);
			let camera_rot = SeaMultivectorPGA3D::new_motor(&rotation_axis, -TAU/4.0);
	
			let camera_trans = camera_rot * camera_pos;
			let camera_trans_inv = camera_trans.inverse();
			
			let scene_point = SeaMultivectorPGA3D::new_point(0.0,3.0,2.0);
			let scene_point_trans = camera_trans * scene_point * camera_trans_inv;
	
			let fp_to_scene = &focal_point&scene_point_trans;
	
			let mut screen_point = fp_to_scene^&image_plane;
			
			screen_point.simplify();

			assert_approx_eq_vec3(
				screen_point.point_to_euclidean(),
				Vector3::new(-1.0, 0.0, 0.0)
			);

			let image_plane = SeaMultivectorPGA3D::new_plane(0.0,0.0,1.0, 0.0);
			let focal_point = SeaMultivectorPGA3D::new_point(0.0, 0.0, -1.0);
	
			let translation_dir = SeaMultivectorPGA3D::new_ideal_line(0.0,0.0,0.0);
			let camera_pos = SeaMultivectorPGA3D::new_motor(&translation_dir, 0.0);
	
			let rotation_axis = SeaMultivectorPGA3D::new_line(0.0,0.0,1.0);
			let camera_rot = SeaMultivectorPGA3D::new_motor(&rotation_axis, -TAU/2.0);
	
			let camera_trans = camera_rot * camera_pos;
			let camera_trans_inv = camera_trans.inverse();
			
			let scene_point = SeaMultivectorPGA3D::new_point(0.0,3.0,2.0);
			let scene_point_trans = camera_trans * scene_point * camera_trans_inv;
	
			let fp_to_scene = &focal_point&scene_point_trans;
	
			let screen_point = fp_to_scene^&image_plane;
	
			assert_approx_eq_vec3(
				screen_point.point_to_euclidean(),
				Vector3::new(0.0, -1.0, 0.0)
			);
		}

		
		#[test]
		fn rotate_and_trans() {
			let image_plane = SeaMultivectorPGA3D::new_plane(0.0,0.0,1.0, 0.0);
			let focal_point = SeaMultivectorPGA3D::new_point(0.0, 0.0, -1.0);
	
			let translation_dir = SeaMultivectorPGA3D::new_ideal_line(1.0,1.0,0.0);
			let camera_pos = SeaMultivectorPGA3D::new_motor(&translation_dir, -3.0);
	
			let rotation_axis = SeaMultivectorPGA3D::new_line(0.0,0.0,1.0);
			let camera_rot = SeaMultivectorPGA3D::new_motor(&rotation_axis, -TAU/4.0);
			
			let camera_trans = camera_rot * camera_pos;
			let camera_trans_inv = camera_trans.inverse();
			
			let scene_point = SeaMultivectorPGA3D::new_point(0.0,3.0,2.0);
			let scene_point_trans = camera_trans * scene_point * camera_trans_inv;
	
			let fp_to_scene = &focal_point&scene_point_trans;
	
			let mut screen_point = fp_to_scene^&image_plane;
			
			screen_point.simplify();

			assert_approx_eq_vec3(
				screen_point.point_to_euclidean(),
				Vector3::new(0.0, -1.0, 0.0)
			);
		}
	
		#[test]
		fn precomputed() {
			let image_plane = SeaMultivectorPGA3D::new_plane(0.0,0.0,1.0, 0.0);
			let focal_point = SeaMultivectorPGA3D::new_point(0.0, 0.0, -1.0);
	
			let translation_dir = SeaMultivectorPGA3D::new_ideal_line(1.0,1.0,0.0);
			let camera_pos = SeaMultivectorPGA3D::new_motor(&translation_dir, -3.0);
	
			let rotation_axis = SeaMultivectorPGA3D::new_line(0.0,0.0,1.0);
			let camera_rot = SeaMultivectorPGA3D::new_motor(&rotation_axis, -TAU/4.0);
			
			let camera_trans = camera_rot * camera_pos;
			let camera_trans_inv = camera_trans.inverse();
			
			let scene_point = SeaMultivectorPGA3D::new_point(0.0,3.0,2.0);
			let scene_point_trans = camera_trans * scene_point * camera_trans_inv;

			let mut screen_point = focal_point&scene_point_trans^image_plane;

			screen_point.simplify();

			assert_approx_eq_vec3(
				screen_point.point_to_euclidean(),
				Vector3::new(0.0, -1.0, 0.0)
			);
		}

		#[test]
		fn sanity_check() {
			let image_plane = SeaMultivectorPGA3D::new_plane(0.0,0.0,1.0, 0.0);
			let mut focal_point = SeaMultivectorPGA3D::new_point(0.0, 0.0, -1.0);
			focal_point.simplify();

			let translation_dir = SeaMultivectorPGA3D::new_ideal_line(0.0,0.0,0.0);
			let camera_pos = SeaMultivectorPGA3D::new_motor(&translation_dir, 0.0);
	
			let rotation_axis = SeaMultivectorPGA3D::new_line(0.0,0.0,0.0);
			let camera_rot = SeaMultivectorPGA3D::new_motor(&rotation_axis, 0.0);
			
			let camera_trans = camera_rot * camera_pos;
			let camera_trans_inv = camera_trans.inverse();
			
			let scene_point = SeaMultivectorPGA3D::new_point(-0.08682,0.4924,0.00);
			let scene_point_trans = camera_trans * scene_point * camera_trans_inv;

			let _screen_point2 = &focal_point & &scene_point_trans;

			let mut screen_point = focal_point&scene_point_trans^image_plane;

			screen_point.simplify();

			assert_approx_eq_vec3(
				screen_point.point_to_euclidean(),
				Vector3::new(-0.08682,0.4924, 0.0)
			);
		}
	}
}