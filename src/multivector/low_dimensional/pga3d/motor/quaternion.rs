use super::*;

// quaternion
impl<T:SeaCoefficient + 'static> SeaMotorPGA3D<T> {
	pub fn from_quaternion(qx: T,qy: T,qz: T,qw: T) -> Self{
		let mut mv = SeaCoreMultivector4D::<T>::new(
			ROTATION_MOTOR_BASIS_BLADES,
			vec![qw, qz, qy, qx]
		);
		mv.remove_zeros();
		Self(mv.into())
	}

	pub fn from_dual_quaternion(dx: T,dy: T,dz: T,dw: T) -> Self{
		let mut mv = SeaCoreMultivector4D::<T>::new(
			TRANSLATION_MOTOR_BASIS_BLADES,
			vec![dx, dy, dz, dw]
		);
		mv.remove_zeros();
		Self(mv.into())
	}

	pub fn from_combined_quaternion(dx: T,dy: T,dz: T,dw: T, qx: T,qy: T,qz: T,qw: T) -> Self{
		let mut mv = SeaCoreMultivector4D::<T>::new(
			COMBINED_MOTOR_BASIS_BLADES,
			vec![qw, qx, qy, qz, dx, dy, dz, dw]
		);
		mv.remove_zeros();
		Self(mv.into())
	}

	pub fn quaternion_angle_positivize(&mut self) {
		if self.scalar().is_sign_negative() {
			self.negate();
		}
	}

	pub fn quaternion_angle_positivized(&self) -> Self{
		let mut out = self.clone();
		out.quaternion_angle_positivize();
		out
	}

	fn to_dual_quaternion(&self) -> DualQuaternion<T> {
		let rw = self.basis_blade_coefficient(&SeaBasisBladeIndexLow::SCALAR);
		let rz = -self.basis_blade_coefficient(&SeaBasisBladeIndexLow::E12);
		let ry = self.basis_blade_coefficient(&SeaBasisBladeIndexLow::E13);
		let rx = -self.basis_blade_coefficient(&SeaBasisBladeIndexLow::E23);
		let tx = self.basis_blade_coefficient(&SeaBasisBladeIndexLow::E14);
		let ty = self.basis_blade_coefficient(&SeaBasisBladeIndexLow::E24);
		let tz = self.basis_blade_coefficient(&SeaBasisBladeIndexLow::E34);
		let tw = self.basis_blade_coefficient(&SeaBasisBladeIndexLow::E1234);

		let real = Quaternion::new(rx, ry, rz, rw);
		let dual = Quaternion::new(tx, ty, tz, tw);

		DualQuaternion::from_real_and_dual(real, dual)
	}
}

#[cfg(test)]
mod tests {
	use std::f64::consts::TAU;
	
	pub use super::super::tests::assert_approx_eq_vec4;
	pub use super::super::tests::assert_approx_eq_mat4;
	use super::multivector::tests::assert_approx_eq_mv;

	use super::*;
		
	#[test]
	fn roll() {
		let rot_mat = Matrix4::from_euler_angles(TAU/3.1,0.0,0.0);
		let rot_mat_quat = SeaMotorPGA3D::<f64>::from(rot_mat);
		let rot_quat = {
			SeaMotorPGA3D::new_rotation_motor(TAU/3.1,0.0, 0.0) // pos neg
		}; 

		assert_approx_eq_mv(rot_mat_quat.into(), rot_quat.into());

		let rot_mat = Matrix4::from_euler_angles(-TAU/3.1,0.0,0.0);
		let rot_mat_quat = SeaMotorPGA3D::<f64>::from(rot_mat);
		let rot_quat = {
			SeaMotorPGA3D::new_rotation_motor(-TAU/3.1,0.0, 0.0) // pos pos
		}; 

		assert_approx_eq_mv(rot_mat_quat.into(), rot_quat.into());
	}

	#[test]
	fn pitch() {
		let rot_mat = Matrix4::from_euler_angles(0.0,TAU/4.9,0.0);
		let rot_mat_quat = SeaMotorPGA3D::<f64>::from(rot_mat);
		let rot_quat = {
			SeaMotorPGA3D::new_rotation_motor(0.0,TAU/4.9, 0.0) // pos pos
		}; 

		assert_approx_eq_mv(rot_mat_quat.into(), rot_quat.into());

		let rot_mat = Matrix4::from_euler_angles(0.0,-TAU/4.9,0.0);
		let rot_mat_quat = SeaMotorPGA3D::<f64>::from(rot_mat);
		let rot_quat = {
			SeaMotorPGA3D::new_rotation_motor(0.0,-TAU/4.9, 0.0) //pos neg
		}; 

		assert_approx_eq_mv(rot_mat_quat.into(), rot_quat.into());
	}

	#[test]
	fn yaw() {
		let rot_mat = Matrix4::from_euler_angles(0.0,0.0 ,TAU/2.8);
		let rot_mat_quat = SeaMotorPGA3D::<f64>::from(rot_mat);
		let rot_quat = {
			SeaMotorPGA3D::new_rotation_motor(0.0,0.0, TAU/2.8) // pos neg
		}; 

		assert_approx_eq_mv(rot_mat_quat.into(), rot_quat.into());

		let rot_mat = Matrix4::from_euler_angles(0.0,0.0 ,TAU/2.8);
		// dbg!(rot_mat);
		let rot_mat_quat = SeaMotorPGA3D::<f64>::from(rot_mat); // neg neg
		let rot_quat = {
			SeaMotorPGA3D::new_rotation_motor(0.0,0.0, TAU/2.8) // pos pos
		}; 

		assert_approx_eq_mv(rot_mat_quat.into(), rot_quat.into());
	}

	#[test]
	fn roll_pitch() {
		let roll_mat = Matrix4::from_euler_angles(TAU/3.1,TAU/4.9,0.0);

		let rot_mat_quat = SeaMotorPGA3D::<f64>::from(roll_mat);
		let rot_quat = {
			SeaMotorPGA3D::new_rotation_motor(0.0,TAU/4.9, 0.0)
			*
			SeaMotorPGA3D::new_rotation_motor(TAU/3.1,0.0, 0.0)
		}; 

		assert_approx_eq_mv(rot_mat_quat.into(), rot_quat.into());

		let roll_mat = Matrix4::from_euler_angles(-TAU/3.1,TAU/4.9,0.0);

		let rot_mat_quat = SeaMotorPGA3D::<f64>::from(roll_mat);
		let rot_quat = {
			SeaMotorPGA3D::new_rotation_motor(0.0,TAU/4.9, 0.0)
			*
			SeaMotorPGA3D::new_rotation_motor(-TAU/3.1,0.0, 0.0)
		}; 

		assert_approx_eq_mv(rot_mat_quat.into(), rot_quat.into());

		let roll_mat = Matrix4::from_euler_angles(TAU/3.1,-TAU/4.9,0.0);

		let rot_mat_quat = SeaMotorPGA3D::<f64>::from(roll_mat);
		let rot_quat = {
			SeaMotorPGA3D::new_rotation_motor(0.0,-TAU/4.9, 0.0)
			*
			SeaMotorPGA3D::new_rotation_motor(TAU/3.1,0.0, 0.0)
		}; 

		assert_approx_eq_mv(rot_mat_quat.into(), rot_quat.into());

		let roll_mat = Matrix4::from_euler_angles(-TAU/3.1,-TAU/4.9,0.0);

		let rot_mat_quat = SeaMotorPGA3D::<f64>::from(roll_mat);
		let rot_quat = {
			SeaMotorPGA3D::new_rotation_motor(0.0,-TAU/4.9, 0.0)
			*
			SeaMotorPGA3D::new_rotation_motor(-TAU/3.1,0.0, 0.0)
		}; 

		assert_approx_eq_mv(rot_mat_quat.into(), rot_quat.into());
	}

	#[test]
	fn roll_yaw() {
		let roll_mat = Matrix4::from_euler_angles(TAU/3.1,0.0,TAU/2.8);

		let rot_mat_quat = SeaMotorPGA3D::<f64>::from(roll_mat);
		let rot_quat = {
			SeaMotorPGA3D::new_rotation_motor(0.0,0.0, TAU/2.8)
			*
			SeaMotorPGA3D::new_rotation_motor(TAU/3.1,0.0, 0.0)
		}; 

		assert_approx_eq_mv(rot_mat_quat.into(), rot_quat.into());

		let roll_mat = Matrix4::from_euler_angles(-TAU/3.1,0.0,TAU/2.8);

		let rot_mat_quat = SeaMotorPGA3D::<f64>::from(roll_mat);
		let rot_quat = {
			SeaMotorPGA3D::new_rotation_motor(0.0,0.0, TAU/2.8)
			*
			SeaMotorPGA3D::new_rotation_motor(-TAU/3.1,0.0, 0.0)
		}; 

		assert_approx_eq_mv(rot_mat_quat.into(), rot_quat.into());

		let roll_mat = Matrix4::from_euler_angles(TAU/3.1,0.0,-TAU/2.8);

		let rot_mat_quat = SeaMotorPGA3D::<f64>::from(roll_mat);
		let rot_quat = {
			SeaMotorPGA3D::new_rotation_motor(0.0,0.0, -TAU/2.8)
			*
			SeaMotorPGA3D::new_rotation_motor(TAU/3.1,0.0, 0.0)
		}; 

		assert_approx_eq_mv(rot_mat_quat.into(), rot_quat.into());

		let roll_mat = Matrix4::from_euler_angles(-TAU/3.1,0.0,-TAU/2.8);

		let rot_mat_quat = SeaMotorPGA3D::<f64>::from(roll_mat);
		let rot_quat = {
			SeaMotorPGA3D::new_rotation_motor(0.0,0.0, -TAU/2.8)
			*
			SeaMotorPGA3D::new_rotation_motor(-TAU/3.1,0.0, 0.0)
		}; 

		assert_approx_eq_mv(rot_mat_quat.into(), rot_quat.into());
	}

	#[test]
	fn pitch_yaw() {
		let roll_mat = Matrix4::from_euler_angles(TAU/3.1,TAU/4.9,0.0);

		let rot_mat_quat = SeaMotorPGA3D::<f64>::from(roll_mat);
		let rot_quat = {
			SeaMotorPGA3D::new_rotation_motor(0.0,TAU/4.9, 0.0)
			*
			SeaMotorPGA3D::new_rotation_motor(TAU/3.1,0.0, 0.0)
		}; 

		assert_approx_eq_mv(rot_mat_quat.into(), rot_quat.into());

		let roll_mat = Matrix4::from_euler_angles(-TAU/3.1,TAU/4.9,0.0);

		let rot_mat_quat = SeaMotorPGA3D::<f64>::from(roll_mat);
		let rot_quat = {
			SeaMotorPGA3D::new_rotation_motor(0.0,TAU/4.9, 0.0)
			*
			SeaMotorPGA3D::new_rotation_motor(-TAU/3.1,0.0, 0.0)
		}; 

		assert_approx_eq_mv(rot_mat_quat.into(), rot_quat.into());

		let roll_mat = Matrix4::from_euler_angles(TAU/3.1,-TAU/4.9,0.0);

		let rot_mat_quat = SeaMotorPGA3D::<f64>::from(roll_mat);
		let rot_quat = {
			SeaMotorPGA3D::new_rotation_motor(0.0,-TAU/4.9, 0.0)
			*
			SeaMotorPGA3D::new_rotation_motor(TAU/3.1,0.0, 0.0)
		}; 

		assert_approx_eq_mv(rot_mat_quat.into(), rot_quat.into());

		let roll_mat = Matrix4::from_euler_angles(-TAU/3.1,-TAU/4.9,0.0);

		let rot_mat_quat = SeaMotorPGA3D::<f64>::from(roll_mat);
		let rot_quat = {
			SeaMotorPGA3D::new_rotation_motor(0.0,-TAU/4.9, 0.0)
			*
			SeaMotorPGA3D::new_rotation_motor(-TAU/3.1,0.0, 0.0)
		}; 

		assert_approx_eq_mv(rot_mat_quat.into(), rot_quat.into());
	}

	#[test]
	fn roll_pitch_yaw() {
		let roll_mat = Matrix4::from_euler_angles(TAU/4.0,TAU/4.0, TAU/4.0);

		let rot_mat_quat = SeaMotorPGA3D::<f64>::from(roll_mat);

		let mut rot_quat = {
			SeaMotorPGA3D::new_rotation_motor(0.0,0.0, TAU/4.0)
			*
			SeaMotorPGA3D::new_rotation_motor(0.0,TAU/4.0, 0.0)
			*
			SeaMotorPGA3D::new_rotation_motor(TAU/4.0,0.0, 0.0)
		}; 

		rot_quat.quaternion_angle_positivize();

		assert_approx_eq_mv(rot_mat_quat.into(), rot_quat.into());

		let roll_mat = Matrix4::from_euler_angles(-TAU/3.1,TAU/4.9,TAU/2.8);

		let rot_mat_quat = SeaMotorPGA3D::<f64>::from(roll_mat);
		let mut rot_quat = {
			SeaMotorPGA3D::new_rotation_motor(0.0,0.0, TAU/2.8)
			*
			SeaMotorPGA3D::new_rotation_motor(0.0,TAU/4.9, 0.0)
			*
			SeaMotorPGA3D::new_rotation_motor(-TAU/3.1,0.0, 0.0)
		}; 

		rot_quat.quaternion_angle_positivize();

		assert_approx_eq_mv(rot_mat_quat.into(), rot_quat.into());

		let roll_mat = Matrix4::from_euler_angles(TAU/3.1,-TAU/4.9,TAU/2.8);

		let rot_mat_quat = SeaMotorPGA3D::<f64>::from(roll_mat);
		let mut rot_quat = {
			SeaMotorPGA3D::new_rotation_motor(0.0,0.0, TAU/2.8)
			*
			SeaMotorPGA3D::new_rotation_motor(0.0,-TAU/4.9, 0.0)
			*
			SeaMotorPGA3D::new_rotation_motor(TAU/3.1,0.0, 0.0)
		};

		rot_quat.quaternion_angle_positivize();

		assert_approx_eq_mv(rot_mat_quat.into(), rot_quat.into());

		let roll_mat = Matrix4::from_euler_angles(-TAU/3.1,-TAU/4.9,TAU/2.8);

		let rot_mat_quat = SeaMotorPGA3D::<f64>::from(roll_mat);
		let mut rot_quat = {
			SeaMotorPGA3D::new_rotation_motor(0.0,0.0, TAU/2.8)
			*
			SeaMotorPGA3D::new_rotation_motor(0.0,-TAU/4.9, 0.0)
			*
			SeaMotorPGA3D::new_rotation_motor(-TAU/3.1,0.0, 0.0)
		}; 
		
		rot_quat.quaternion_angle_positivize();

		assert_approx_eq_mv(rot_mat_quat.into(), rot_quat.into());

		let roll_mat = Matrix4::from_euler_angles(TAU/3.1,TAU/4.9,-TAU/2.8);

		let rot_mat_quat = SeaMotorPGA3D::<f64>::from(roll_mat);
		let mut rot_quat = {
			SeaMotorPGA3D::new_rotation_motor(0.0,0.0, -TAU/2.8)
			*
			SeaMotorPGA3D::new_rotation_motor(0.0,TAU/4.9, 0.0)
			*
			SeaMotorPGA3D::new_rotation_motor(TAU/3.1,0.0, 0.0)
		}; 
		
		rot_quat.quaternion_angle_positivize();

		assert_approx_eq_mv(rot_mat_quat.into(), rot_quat.into());

		let roll_mat = Matrix4::from_euler_angles(-TAU/3.1,TAU/4.9,-TAU/2.8);

		let rot_mat_quat = SeaMotorPGA3D::<f64>::from(roll_mat);
		let mut rot_quat = {
			SeaMotorPGA3D::new_rotation_motor(0.0,0.0, -TAU/2.8)
			*
			SeaMotorPGA3D::new_rotation_motor(0.0,TAU/4.9, 0.0)
			*
			SeaMotorPGA3D::new_rotation_motor(-TAU/3.1,0.0, 0.0)
		}; 
		
		rot_quat.quaternion_angle_positivize();

		assert_approx_eq_mv(rot_mat_quat.into(), rot_quat.into());

		let roll_mat = Matrix4::from_euler_angles(TAU/3.1,-TAU/4.9,-TAU/2.8);

		let rot_mat_quat = SeaMotorPGA3D::<f64>::from(roll_mat);
		let mut rot_quat = {
			SeaMotorPGA3D::new_rotation_motor(0.0,0.0, -TAU/2.8)
			*
			SeaMotorPGA3D::new_rotation_motor(0.0,-TAU/4.9, 0.0)
			*
			SeaMotorPGA3D::new_rotation_motor(TAU/3.1,0.0, 0.0)
		}; 
		
		rot_quat.quaternion_angle_positivize();

		assert_approx_eq_mv(rot_mat_quat.into(), rot_quat.into());

		let roll_mat = Matrix4::from_euler_angles(-TAU/3.1,-TAU/4.9,-TAU/2.8);

		let rot_mat_quat = SeaMotorPGA3D::<f64>::from(roll_mat);
		let mut rot_quat = {
			SeaMotorPGA3D::new_rotation_motor(0.0,0.0, -TAU/2.8)
			*
			SeaMotorPGA3D::new_rotation_motor(0.0,-TAU/4.9, 0.0)
			*
			SeaMotorPGA3D::new_rotation_motor(-TAU/3.1,0.0, 0.0)
		}; 
		
		rot_quat.quaternion_angle_positivize();

		assert_approx_eq_mv(rot_mat_quat.into(), rot_quat.into());
	}
}
