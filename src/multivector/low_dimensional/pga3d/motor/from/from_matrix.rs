use super::*;

impl<T:SeaCoefficient + 'static> From<Matrix3<T>> for SeaMotorPGA3D<T> {
	fn from(rot: Matrix3<T>) -> Self {
		Self::from(&rot)
	}		
}
impl<T:SeaCoefficient + 'static> From<&Matrix3<T>> for SeaMotorPGA3D<T> {
	fn from(rot: &Matrix3<T>) -> Self {

		let two = T::two();
		let four = T::four();

		let trace = rot.trace();
		let mut quat = [T::zero(); 4];
		if trace > T::zero() {
			let s = T::one() / (two*(trace + T::one()).sqrt());
			quat[0] = ( rot[(2, 1)] - rot[(1, 2)] ) * s;
			quat[1] = ( rot[(0, 2)] - rot[(2, 0)] ) * s;
			quat[2] = ( rot[(1, 0)] - rot[(0, 1)] ) * s; 
			quat[3] = T::one() / (four*s);
		}
		else if rot[(0, 0)] > rot[(1, 1)] && rot[(0, 0)] > rot[(2,2)] {
			let s = two * ( T::one() + rot[(0,0)] - rot[(1,1)] - rot[(2,2)]).sqrt();
			quat[0] = s / four;
			quat[1] = ( rot[(0, 1)] + rot[(1, 0)] ) / s; 
			quat[2] = ( rot[(0, 2)] + rot[(2, 0)] ) / s;
			quat[3] = ( rot[(2, 1)] - rot[(1,2)] ) / s;
			
		}
		else if rot[(1, 1)] > rot[(2, 2)] {
			let s = two * ( T::one() + rot[(1,1)] - rot[(0,0)] - rot[(2,2)]).sqrt();
			quat[0] = ( rot[(0, 1)] + rot[(1, 0)] ) / s;
			quat[1] = s / four; 
			quat[2] = ( rot[(1, 2)] + rot[(2,1)] ) / s;
			quat[3] = ( rot[(0, 2)] - rot[(2, 0)] ) / s;
		}
		else {
			let s = two * ( T::one() + rot[(2,2)] - rot[(0,0)] - rot[(1,1)]).sqrt();
			quat[0] = ( rot[(0, 2)] + rot[(2, 0)] ) / s;
			quat[1] = ( rot[(1, 2)] + rot[(2, 1)] ) / s;
			quat[2] = s / four;
			quat[3] = ( rot[(1, 0)] - rot[(0, 1)] ) / s; 
		}

		let mut q = Self::from_quaternion(-quat[0], quat[1], -quat[2], quat[3]);
		q.quaternion_angle_positivize();
		q
	}
}
impl<T:SeaCoefficient + 'static> From<Matrix4<T>> for SeaMotorPGA3D<T> {
	fn from(value: Matrix4<T>) -> Self {
		Self::from(&value)
	}
}
impl<T:SeaCoefficient + 'static> From<&Matrix4<T>> for SeaMotorPGA3D<T> {
	fn from(value: &Matrix4<T>) -> Self {
		let rot_mat = Matrix3x4::identity() * (value * Matrix4x3::identity());
		let rot_mot = Self::from(rot_mat);

		let trans_mot = Self::new_translation_motor(value[(0,3)], value[(1,3)], value[(2,3)]);

		trans_mot * rot_mot
	}
}

#[cfg(test)]
mod tests {
	use std::f64::consts::TAU;
	
	use super::super::tests::assert_approx_eq_vec4;
	use super::super::tests::assert_approx_eq_mat4;
	use super::multivector::tests::assert_approx_eq_mv;

	use super::*;
	
	#[test]
	fn from_matrix3() {
		let rot_mat = Matrix4::from_euler_angles(TAU/4.0, 0.0, 0.0);
		let rot_quat = SeaMotorPGA3D::<f64>::from(rot_mat);

		let p = SeaMultivectorPGA3D::new_point(1.0, 2.0, 3.0);

		let p_vec = p.point_to_homogeneous();

		let p_mat = rot_mat * p_vec;
		let p_mv = rot_quat / p;

		assert_approx_eq_vec4(
			p_mv.point_to_homogeneous(),
			p_mat
		);

		let rot_mat = Matrix4::from_euler_angles(0.0, TAU/4.0, 0.0);
		let rot_quat = SeaMotorPGA3D::<f64>::from(rot_mat);

		let p = SeaMultivectorPGA3D::new_point(1.0, 2.0, 3.0);

		let p_vec = p.point_to_homogeneous();

		let p_mat = rot_mat * p_vec;
		let p_mv = rot_quat / p;

		assert_approx_eq_vec4(
			p_mv.point_to_homogeneous(),
			p_mat
		);

		let rot_mat = Matrix4::from_euler_angles(0.0, 0.0, TAU/4.0);
		let rot_quat = SeaMotorPGA3D::<f64>::from(rot_mat);

		let p = SeaMultivectorPGA3D::new_point(1.0, 2.0, 3.0);

		let p_vec = p.point_to_homogeneous();

		let p_mat = rot_mat * p_vec;
		let p_mv = rot_quat / p;

		assert_approx_eq_vec4(
			p_mv.point_to_homogeneous(),
			p_mat
		);

		let rot_mat = Matrix4::from_euler_angles(-TAU/3.1, TAU/4.9,-TAU/2.8);
		let rot_quat = SeaMotorPGA3D::<f64>::from(rot_mat);

		let p = SeaMultivectorPGA3D::new_point(1.0, 2.0, 3.0);

		let p_vec = p.point_to_homogeneous();

		let p_mat = rot_mat * p_vec;
		let p_mv = rot_quat / p;

		assert_approx_eq_vec4(
			p_mv.point_to_homogeneous(),
			p_mat
		);
	}

	#[test]
	fn from_matrix4() {
			let rot_mat = Matrix4::from_euler_angles(TAU/3.1, TAU/4.9,TAU/2.8);
			let trans_mat = Matrix4::new_translation(&Vector3::new(1.2, 2.3, 3.4));
			let full_mat = trans_mat * rot_mat;
			let trans_mot = SeaMotorPGA3D::<f64>::from(full_mat);

			let p = SeaMultivectorPGA3D::new_point(1.0, 2.0, 3.0);

			let p_vec = p.point_to_homogeneous();

			let p_mat = full_mat * p_vec;
			let p_mv = trans_mot / p;

			assert_approx_eq_vec4(
				p_mv.point_to_homogeneous(),
				p_mat
			);
		}

}
