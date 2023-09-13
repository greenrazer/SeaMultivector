use std::ops::{Deref, DerefMut, Shl, Shr, BitOr, BitXor, Mul, Not, BitAnd, Index, AddAssign, Add, Neg, Sub, Div, SubAssign, MulAssign, DivAssign};

use crate::coefficient::SeaCoefficient;
use crate::multivector::metric::SeaMetric;
use crate::multivector::multivector_traits::{Multivector, MultivectorMut, Versor, Motor};
use crate::multivector::multivector_product_traits::{ScalarProduct, LeftContraction, RightContraction, InnerProduct, OuterProduct, GeometricProduct, RegressiveProduct};
use crate::multivector::grade::SeaGrade;
use crate::multivector::basis_blade_index::SeaBasisBladeIndexLow;

use super::basis_blade_set::SeaBasisBladeSet;
use super::core::SeaLowDimensionalSparseMultivector;

mod multivector;
mod motor;

pub use multivector::SeaMultivectorPGA3D;
pub use motor::SeaMotorPGA3D;


// impl<T:SeaCoefficient> From<SeaMotorPGA3D<T>> for SeaMultivectorPGA3D<T> {
// 	fn from(value: SeaMotorPGA3D<T>) -> Self {
// 		value.0
// 	}
// }
// impl<T:SeaCoefficient> From<&SeaMotorPGA3D<T>> for SeaMultivectorPGA3D<T> {
// 	fn from(value: &SeaMotorPGA3D<T>) -> Self {
// 		value.0.clone()
// 	}
// }

// impl<T:SeaCoefficient> Neg for &SeaMultivectorPGA3D<T> {
// 	type Output = SeaMultivectorPGA3D<T>;

// 	fn neg(self) -> Self::Output {
// 		SeaMultivectorPGA3D(-&self.0)
// 	}
// }
// impl<T:SeaCoefficient> Neg for SeaMultivectorPGA3D<T> {
// 	type Output = SeaMultivectorPGA3D<T>;

// 	fn neg(self) -> Self::Output {
// 		-&self
// 	}
// }

// impl<T:SeaCoefficient> AddAssign<T> for SeaMultivectorPGA3D<T> {
//     fn add_assign(&mut self, rhs: T){
// 		self.0 += rhs
// 	}
// }
// impl<T:SeaCoefficient> Add<T> for SeaMultivectorPGA3D<T> {
//     type Output = SeaMultivectorPGA3D<T>;

//     fn add(self, rhs: T) -> Self::Output {
// 		&self + rhs
// 	}
// }
// impl<T:SeaCoefficient> Add<T> for &SeaMultivectorPGA3D<T> {
//     type Output = SeaMultivectorPGA3D<T>;

//     fn add(self, rhs: T) -> Self::Output {
// 		SeaMultivectorPGA3D(&self.0 + rhs)
// 	}
// }

// impl<T:SeaCoefficient> SubAssign<T> for SeaMultivectorPGA3D<T> {
//     fn sub_assign(&mut self, rhs: T){
// 		self.0 -= rhs
// 	}
// }
// impl<T:SeaCoefficient> Sub<T> for SeaMultivectorPGA3D<T> {
//     type Output = SeaMultivectorPGA3D<T>;

//     fn sub(self, rhs: T) -> Self::Output {
// 		&self - rhs
// 	}
// }
// impl<T:SeaCoefficient> Sub<T> for &SeaMultivectorPGA3D<T> {
//     type Output = SeaMultivectorPGA3D<T>;

//     fn sub(self, rhs: T) -> Self::Output {
// 		SeaMultivectorPGA3D(&self.0 - rhs)
// 	}
// }

// impl<T:SeaCoefficient> MulAssign<T> for SeaMultivectorPGA3D<T> {
//     fn mul_assign(&mut self, rhs: T){
// 		self.0 *= rhs
// 	}
// }
// impl<T:SeaCoefficient> Mul<T> for SeaMultivectorPGA3D<T> {
//     type Output = SeaMultivectorPGA3D<T>;

//     fn mul(self, rhs: T) -> Self::Output {
// 		&self * rhs
// 	}
// }
// impl<T:SeaCoefficient> Mul<T> for &SeaMultivectorPGA3D<T> {
//     type Output = SeaMultivectorPGA3D<T>;

//     fn mul(self, rhs: T) -> Self::Output {
// 		SeaMultivectorPGA3D(&self.0 * rhs)
// 	}
// }

// impl<T:SeaCoefficient> DivAssign<T> for SeaMultivectorPGA3D<T> {
//     fn div_assign(&mut self, rhs: T){
// 		self.0 /= rhs
// 	}
// }
// impl<T:SeaCoefficient> Div<T> for SeaMultivectorPGA3D<T> {
//     type Output = SeaMultivectorPGA3D<T>;

//     fn div(self, rhs: T) -> Self::Output {
// 		&self / rhs
// 	}
// }
// impl<T:SeaCoefficient> Div<T> for &SeaMultivectorPGA3D<T> {
//     type Output = SeaMultivectorPGA3D<T>;

//     fn div(self, rhs: T) -> Self::Output {
// 		SeaMultivectorPGA3D(&self.0 / rhs)
// 	}
// }


#[cfg(test)]
mod tests {
	use std::f64::consts::TAU;

	use super::*;

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
			let scene_point_trans = camera_trans.into() * scene_point * camera_trans_inv.into();
	
			let fp_to_scene = focal_point&scene_point_trans;
	
			let screen_point = fp_to_scene^image_plane;
	
			assert_eq!(
				screen_point.point_to_euclidean(),
				[0.0, 1.0, 0.0]
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
			let scene_point_trans = camera_trans.into() * scene_point * camera_trans_inv.into();
	
			let fp_to_scene = focal_point&scene_point_trans;
	
			let screen_point = fp_to_scene^image_plane;
	
			assert_eq!(
				screen_point.point_to_euclidean(),
				[-1.0, 0.0, 0.0]
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
			let scene_point_trans = camera_trans.into() * scene_point * camera_trans_inv.into();
	
			let fp_to_scene = &focal_point&scene_point_trans;
	
			let mut screen_point = fp_to_scene^&image_plane;
			
			screen_point.simplify();

			assert_eq!(
				screen_point.point_to_euclidean(),
				[-1.0, 0.0, 0.0]
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
			let scene_point_trans = camera_trans.into() * scene_point * camera_trans_inv.into();
	
			let fp_to_scene = &focal_point&scene_point_trans;
	
			let screen_point = fp_to_scene^&image_plane;
	
			assert_eq!(
				screen_point.point_to_euclidean(),
				[0.0, -1.0, 0.0]
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
			let scene_point_trans = camera_trans.into() * scene_point * camera_trans_inv.into();
	
			let fp_to_scene = &focal_point&scene_point_trans;
	
			let mut screen_point = fp_to_scene^&image_plane;
			
			screen_point.simplify();

			assert_eq!(
				screen_point.point_to_euclidean(),
				[0.0, -1.0, 0.0]
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
			let scene_point_trans = camera_trans.into() * scene_point * camera_trans_inv.into();

			let mut screen_point = focal_point&scene_point_trans^image_plane;

			screen_point.simplify();

			assert_eq!(
				screen_point.point_to_euclidean(),
				[0.0, -1.0, 0.0]
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
			let scene_point_trans = camera_trans.into() * scene_point * camera_trans_inv.into();

			let _screen_point2 = &focal_point & &scene_point_trans;

			let mut screen_point = focal_point&scene_point_trans^image_plane;

			screen_point.simplify();

			assert_eq!(
				screen_point.point_to_euclidean(),
				[-0.08682,0.4924, 0.0]
			);
		}
	}
}