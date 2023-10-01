use std::{f64::consts::TAU, arch::x86_64::_MM_FLUSH_ZERO_MASK};

use super::*;

mod scalar_operations;
mod linear_operations;
mod non_linear_operations;
mod non_geometric_linear_operations;


#[derive(Debug, Clone, PartialEq)]
pub struct SeaMotorPGA3D<T:SeaCoefficient>(pub SeaMultivectorPGA3D<T>);
impl<T:SeaCoefficient> SeaMotorPGA3D<T> {
	pub fn identity() -> Self {
		Self::from(T::one())
	}

	pub fn from_quaternion(qx: T,qy: T,qz: T,qw: T) -> Self{
		// println!("{} {} {} {}", qx, qy, qz, qw);
		// println!("{} {} {} {}", -qx, qy, -qz, qw);
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

	pub fn quaternion_mirror_y_axis(&mut self){
		
	}

	pub fn quaternion_mirrored_y_axis(&self) -> Self{
		let mut out = self.clone();
		out.quaternion_mirror_y_axis();
		out
	}
}

impl<T:SeaCoefficient> SeaMotorPGA3D<T> {
	pub fn log(&self) -> SeaMultivectorPGA3D<T> {
		let b = &self.grade(&SeaGrade(2));
		let s = (-b | b).scalar().sqrt();

		// pure translation
		if s.is_zero(){
			return (&**self)/self.scalar() - T::one();
		}

		// turn
		if self.scalar().is_zero() {
			let tau_div_4 = T::from(TAU/4.0).unwrap();
			return -self.pseudoscalar()/s + tau_div_4;
		}

		let two = T::from(2.0).unwrap();
		let p = &( (-b ^ b) / (two*s) );
		
		let lhs = (p / self.scalar()) + (s / self.scalar()).atan();
		let rhs = b*((-p+s)/(s*s));

		lhs * rhs
	}
}
impl<T:SeaCoefficient + 'static> SeaMotorPGA3D<T> {

	pub fn to_homogeneous_transformation_matrices(&self) -> (Matrix4<T>, Matrix4<T>) {
		let rw = self.basis_blade_coefficient(&SeaBasisBladeIndexLow::SCALAR);
		let rz = -self.basis_blade_coefficient(&SeaBasisBladeIndexLow::E12);
		let ry = self.basis_blade_coefficient(&SeaBasisBladeIndexLow::E13);
		let rx = -self.basis_blade_coefficient(&SeaBasisBladeIndexLow::E23);
		let tx = self.basis_blade_coefficient(&SeaBasisBladeIndexLow::E14);
		let ty = self.basis_blade_coefficient(&SeaBasisBladeIndexLow::E24);
		let tz = self.basis_blade_coefficient(&SeaBasisBladeIndexLow::E34);
		let tw = self.basis_blade_coefficient(&SeaBasisBladeIndexLow::E1234);

		let zero = T::zero();
		let one = T::one();
		let two = T::from(2).unwrap();

		let a: Matrix4<T> = Matrix4::new(
			one - two*(ry*ry + rz*rz), two*rx*ry, two*rx*rz, two*(ry*tz - rz*ty),
			two*rx*ry, one - two*(rx*rx + rz*rz), two*ry*rz, two*(rz*tx - rx*tz),
			two*rx*rz, two*ry*rz, one - two*(rx*rx + ry*ry), two*(rx*ty - ry*tx),
			zero, zero, zero, one
		);

		let b: Matrix4<T> = Matrix4::new(
			zero, -two*rz*rw, two*ry*rw, two*(rw*tx - rx*tw),
			two*rz*rw, zero, -two*rx*rw, two*(rw*ty - ry*tw),
			-two*ry*rw, two*rx*rw, zero, two*(rw*tz - rz*tw),
			zero, zero, zero, zero
		);

		let m = a+b;
		let m_inv = a-b;

		(m, m_inv)
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

impl<T: SeaCoefficient> Deref for SeaMotorPGA3D<T> {
	type Target = SeaMultivectorPGA3D<T>;
	fn deref(&self) -> &Self::Target {
		&self.0
	}
}
impl<T: SeaCoefficient> DerefMut for SeaMotorPGA3D<T> {
	fn deref_mut(&mut self) -> &mut Self::Target {
		&mut self.0
	}
}

impl<T:SeaCoefficient> Default for SeaMotorPGA3D<T> {
	fn default() -> Self {
		Self(SeaMultivectorPGA3D::<T>::default())
	}
}

impl<T:SeaCoefficient> Motor<T> for SeaMotorPGA3D<T> {}

impl<T:SeaCoefficient> Versor<T> for SeaMotorPGA3D<T> {
	fn inverse(&self) -> Self {
		Self((**self).reversion_conjugated())
	}
}

impl<T:SeaCoefficient> Multivector<T> for SeaMotorPGA3D<T> {
	type MultivectorIndexType = SeaBasisBladeIndexLow;
	type GeneralMultivectorOutput = SeaMultivectorPGA3D<T>;

	fn basis_blade_coefficient(&self, index: &Self::MultivectorIndexType) -> T {
		(**self).basis_blade_coefficient(index)
	}

	fn scalar(&self) -> T {
		(**self).scalar()
	}

	fn pseudoscalar(&self) -> Self::GeneralMultivectorOutput {
		(**self).pseudoscalar()
	}

	fn grade(&self, grade: &SeaGrade) -> Self::GeneralMultivectorOutput {
		(**self).grade(grade)
	}

	fn exp(&self) -> Self {
		Self((**self).exp())
	}
}
impl<T:SeaCoefficient> MultivectorMut<T> for SeaMotorPGA3D<T> {
	fn negate(&mut self) {
		self.0.negate();
	}
}

impl<T: SeaCoefficient> MultivectorConstMetric<T> for SeaMotorPGA3D<T> {
	const METRIC: SeaMetric = SeaMultivectorPGA3D::<T>::METRIC;
}

impl<T:SeaCoefficient> MultivectorIter for SeaMotorPGA3D<T> {
	fn iter(&self) -> <&Self as IntoIterator>::IntoIter {
		self.into_iter()
	}
}

impl<T:SeaCoefficient> MultivectorIterMut for SeaMotorPGA3D<T> {
	fn iter_mut(&mut self) -> <&mut Self as IntoIterator>::IntoIter{
		self.into_iter()
	}
}

impl<'a, T:SeaCoefficient> IntoIterator for &'a SeaMotorPGA3D<T> {
	type IntoIter =  SeaCoreSparseMultivector4DIter<'a, T>;
	type Item = SeaBasisBladeLowRef<'a, T>;

	fn into_iter(self) -> Self::IntoIter {
		self.0.iter()
	}

}
impl<'a, T:SeaCoefficient> IntoIterator for &'a mut SeaMotorPGA3D<T> {
	type IntoIter =  SeaCoreSparseMultivector4DIterMut<'a, T>;
	type Item = SeaBasisBladeLowRefMut<'a, T>;

	fn into_iter(self) -> Self::IntoIter {
		self.0.iter_mut()
	}

}

impl<'a, T:SeaCoefficient> MultivectorGradeIter for SeaMotorPGA3D<T> {
	fn grade_iter(&self, grade: &SeaGrade) -> <&Self as MultivectorGradeIntoIterator>::MultivectorGradeIntoIter {
		self.grade_into_iter(grade)
	}
}

impl<'a, T:SeaCoefficient> MultivectorGradeIterMut for SeaMotorPGA3D<T> {
	fn grade_iter_mut(&mut self, grade: &SeaGrade) -> <&mut Self as MultivectorGradeIntoIterator>::MultivectorGradeIntoIter {
		self.grade_into_iter(grade)
	}
}

impl<'a, T:SeaCoefficient> MultivectorGradeIntoIterator for &'a SeaMotorPGA3D<T> {
	type MultivectorGradeIntoIter = SeaCoreSparseMultivector4DIter<'a, T>;

	fn grade_into_iter(self, grade: &SeaGrade) -> Self::MultivectorGradeIntoIter {
		self.0.grade_iter(grade)
	}
}
impl<'a, T:SeaCoefficient> MultivectorGradeIntoIterator for &'a mut SeaMotorPGA3D<T> {
	type MultivectorGradeIntoIter = SeaCoreSparseMultivector4DIterMut<'a, T>;

	fn grade_into_iter(self, grade: &SeaGrade) -> Self::MultivectorGradeIntoIter {
		self.0.grade_iter_mut(grade)
	}
}

impl<T:SeaCoefficient> From<T> for SeaMotorPGA3D<T> {
	fn from(value: T) -> Self {
		SeaMotorPGA3D(SeaMultivectorPGA3D::<T>::new_lifted_scalar(value))	
	}
}

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

		let trans_mot = SeaMultivectorPGA3D::new_translation_motor(value[(0,3)], value[(1,3)], value[(2,3)]);

		trans_mot * rot_mot
	}
}

impl<T:SeaCoefficient> Neg for SeaMotorPGA3D<T> {
	type Output = Self;

	fn neg(self) -> Self::Output {
		Self(-&*self)
	}
}
impl<T:SeaCoefficient> Neg for &SeaMotorPGA3D<T> {
	type Output = SeaMotorPGA3D<T>;

	fn neg(self) -> Self::Output {
		SeaMotorPGA3D::<T>(-&**self)
	}
}

impl<T:SeaCoefficient> Not for SeaMotorPGA3D<T> {
	type Output = SeaMotorPGA3D<T>;

	fn not(self) -> Self::Output {
		self.dual()
	}
}
impl<T:SeaCoefficient> Not for &SeaMotorPGA3D<T> {
	type Output = SeaMotorPGA3D<T>;

	fn not(self) -> Self::Output {
		self.dual()
	}
}

#[cfg(test)]
mod tests {
	use std::f64::consts::TAU;
	
	use super::super::tests::assert_approx_eq_vec4;
	use super::super::tests::assert_approx_eq_mat4;
	use super::multivector::tests::assert_approx_eq_mv;

	use super::*;

	mod matrix{
		use super::*;

		#[test]
		fn translation() {
			let p = SeaMultivectorPGA3D::new_point(1.0, 0.0, 0.0);

			let l = SeaMultivectorPGA3D::new_ideal_line(0.0,1.0,0.0);
			let r = SeaMultivectorPGA3D::new_motor(&l, TAU/4.0);

			let p_trans = &r * &p * r.inverse();

			let(m, m_inv) = r.to_homogeneous_transformation_matrices();

			let p_homo_trans = m*p.point_to_homogeneous();

			assert_approx_eq_vec4(
				p_homo_trans,
				p_trans.point_to_homogeneous()
			);
			assert_eq!(
				m * m_inv,
				Matrix4::identity()
			);
			assert_eq!(
				m * m_inv,
				m_inv * m
			);


			let p = SeaMultivectorPGA3D::new_point(5.0, -3.0, 7.7);

			let l = SeaMultivectorPGA3D::new_ideal_line(1.0,3.0,2.0);
			let r = SeaMultivectorPGA3D::new_motor(&l, TAU/3.2);

			let p_trans = &r * &p * r.inverse();

			let(m, m_inv) = r.to_homogeneous_transformation_matrices();

			let p_homo_trans = m*p.point_to_homogeneous();

			assert_approx_eq_vec4(
				p_homo_trans,
				p_trans.point_to_homogeneous()
			);
			assert_eq!(
				m * m_inv,
				Matrix4::identity()
			);
			assert_eq!(
				m * m_inv,
				m_inv * m
			);

			let p = SeaMultivectorPGA3D::new_point(5.0, -3.0, 7.7);

			let l = SeaMultivectorPGA3D::new_ideal_line(1.0,3.0,2.0);
			let r1 = SeaMultivectorPGA3D::new_motor(&l, TAU/3.2);
			let r2 = SeaMultivectorPGA3D::new_motor(&l, -TAU/3.2);

			let p_trans1 = &r1 * &p * r1.inverse();
			let p_trans2 = &r2 * &p * r2.inverse();

			let(m1, m1_inv) = r1.to_homogeneous_transformation_matrices();

			let p_homo_trans1 = m1*p.point_to_homogeneous();
			let p_homo_trans2 = m1_inv*p.point_to_homogeneous();

			assert_approx_eq_vec4(
				p_homo_trans1,
				p_trans1.point_to_homogeneous()
			);

			assert_approx_eq_vec4(
				p_homo_trans2,
				p_trans2.point_to_homogeneous()
			);
		}

		#[test]
		fn rotation() {
			let p = SeaMultivectorPGA3D::new_point(1.0, 0.0, 0.0);

			let l = SeaMultivectorPGA3D::new_line(0.0,1.0,0.0);
			let r = SeaMultivectorPGA3D::new_motor(&l, TAU/4.0);

			let p_trans = &r * &p * r.inverse();

			let(m, m_inv) = r.to_homogeneous_transformation_matrices();

			let p_homo_trans = m*p.point_to_homogeneous();

			assert_approx_eq_vec4(
				p_homo_trans,
				p_trans.point_to_homogeneous()
			);
			assert_approx_eq_mat4(
				m * m_inv,
				Matrix4::identity()
			);
			assert_approx_eq_mat4(
				m * m_inv,
				m_inv * m
			);


			let p = SeaMultivectorPGA3D::new_point(5.0, -3.0, 7.7);

			let l = SeaMultivectorPGA3D::new_line(1.0,3.0,2.0);
			let r = SeaMultivectorPGA3D::new_motor(&l, TAU/3.2);

			let p_trans = &r * &p * r.inverse();

			let(m, m_inv) = r.to_homogeneous_transformation_matrices();

			let p_homo_trans = m*p.point_to_homogeneous();

			assert_approx_eq_vec4(
				p_homo_trans,
				p_trans.point_to_homogeneous()
			);
			assert_approx_eq_mat4(
				m * m_inv,
				Matrix4::identity()
			);
			assert_approx_eq_mat4(
				m * m_inv,
				m_inv * m
			);

			let p = SeaMultivectorPGA3D::new_point(5.0, -3.0, 7.7);

			let l = SeaMultivectorPGA3D::new_line(1.0,3.0,2.0);
			let r1 = SeaMultivectorPGA3D::new_motor(&l, TAU/3.2);
			let r2 = SeaMultivectorPGA3D::new_motor(&l, -TAU/3.2);

			let p_trans1 = &r1 * &p * r1.inverse();
			let p_trans2 = &r2 * &p * r2.inverse();

			let(m1, m1_inv) = r1.to_homogeneous_transformation_matrices();

			let p_homo_trans1 = m1*p.point_to_homogeneous();
			let p_homo_trans2 = m1_inv*p.point_to_homogeneous();

			assert_approx_eq_vec4(
				p_homo_trans1,
				p_trans1.point_to_homogeneous()
			);

			assert_approx_eq_vec4(
				p_homo_trans2,
				p_trans2.point_to_homogeneous()
			);
		}

		
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

	mod quaternion{
		use super::*;
		
		#[test]
		fn roll() {
			let rot_mat = Matrix4::from_euler_angles(TAU/3.1,0.0,0.0);
			let rot_mat_quat = SeaMotorPGA3D::<f64>::from(rot_mat);
			let rot_quat = {
				SeaMultivectorPGA3D::new_rotation_motor(TAU/3.1,0.0, 0.0) // pos neg
			}; 

			assert_approx_eq_mv(rot_mat_quat.into(), rot_quat.into());

			let rot_mat = Matrix4::from_euler_angles(-TAU/3.1,0.0,0.0);
			let rot_mat_quat = SeaMotorPGA3D::<f64>::from(rot_mat);
			let rot_quat = {
				SeaMultivectorPGA3D::new_rotation_motor(-TAU/3.1,0.0, 0.0) // pos pos
			}; 

			assert_approx_eq_mv(rot_mat_quat.into(), rot_quat.into());
		}

		#[test]
		fn pitch() {
			let rot_mat = Matrix4::from_euler_angles(0.0,TAU/4.9,0.0);
			let rot_mat_quat = SeaMotorPGA3D::<f64>::from(rot_mat);
			let rot_quat = {
				SeaMultivectorPGA3D::new_rotation_motor(0.0,TAU/4.9, 0.0) // pos pos
			}; 

			assert_approx_eq_mv(rot_mat_quat.into(), rot_quat.into());

			let rot_mat = Matrix4::from_euler_angles(0.0,-TAU/4.9,0.0);
			let rot_mat_quat = SeaMotorPGA3D::<f64>::from(rot_mat);
			let rot_quat = {
				SeaMultivectorPGA3D::new_rotation_motor(0.0,-TAU/4.9, 0.0) //pos neg
			}; 

			assert_approx_eq_mv(rot_mat_quat.into(), rot_quat.into());
		}

		#[test]
		fn yaw() {
			let rot_mat = Matrix4::from_euler_angles(0.0,0.0 ,TAU/2.8);
			let rot_mat_quat = SeaMotorPGA3D::<f64>::from(rot_mat);
			let rot_quat = {
				SeaMultivectorPGA3D::new_rotation_motor(0.0,0.0, TAU/2.8) // pos neg
			}; 

			assert_approx_eq_mv(rot_mat_quat.into(), rot_quat.into());

			let rot_mat = Matrix4::from_euler_angles(0.0,0.0 ,TAU/2.8);
			// dbg!(rot_mat);
			let rot_mat_quat = SeaMotorPGA3D::<f64>::from(rot_mat); // neg neg
			let rot_quat = {
				SeaMultivectorPGA3D::new_rotation_motor(0.0,0.0, TAU/2.8) // pos pos
			}; 

			assert_approx_eq_mv(rot_mat_quat.into(), rot_quat.into());
		}

		#[test]
		fn roll_pitch() {
			let roll_mat = Matrix4::from_euler_angles(TAU/3.1,TAU/4.9,0.0);

			let rot_mat_quat = SeaMotorPGA3D::<f64>::from(roll_mat);
			let rot_quat = {
				SeaMultivectorPGA3D::new_rotation_motor(0.0,TAU/4.9, 0.0)
				*
				SeaMultivectorPGA3D::new_rotation_motor(TAU/3.1,0.0, 0.0)
			}; 

			assert_approx_eq_mv(rot_mat_quat.into(), rot_quat.into());

			let roll_mat = Matrix4::from_euler_angles(-TAU/3.1,TAU/4.9,0.0);

			let rot_mat_quat = SeaMotorPGA3D::<f64>::from(roll_mat);
			let rot_quat = {
				SeaMultivectorPGA3D::new_rotation_motor(0.0,TAU/4.9, 0.0)
				*
				SeaMultivectorPGA3D::new_rotation_motor(-TAU/3.1,0.0, 0.0)
			}; 

			assert_approx_eq_mv(rot_mat_quat.into(), rot_quat.into());

			let roll_mat = Matrix4::from_euler_angles(TAU/3.1,-TAU/4.9,0.0);

			let rot_mat_quat = SeaMotorPGA3D::<f64>::from(roll_mat);
			let rot_quat = {
				SeaMultivectorPGA3D::new_rotation_motor(0.0,-TAU/4.9, 0.0)
				*
				SeaMultivectorPGA3D::new_rotation_motor(TAU/3.1,0.0, 0.0)
			}; 

			assert_approx_eq_mv(rot_mat_quat.into(), rot_quat.into());

			let roll_mat = Matrix4::from_euler_angles(-TAU/3.1,-TAU/4.9,0.0);

			let rot_mat_quat = SeaMotorPGA3D::<f64>::from(roll_mat);
			let rot_quat = {
				SeaMultivectorPGA3D::new_rotation_motor(0.0,-TAU/4.9, 0.0)
				*
				SeaMultivectorPGA3D::new_rotation_motor(-TAU/3.1,0.0, 0.0)
			}; 

			assert_approx_eq_mv(rot_mat_quat.into(), rot_quat.into());
		}

		#[test]
		fn roll_yaw() {
			let roll_mat = Matrix4::from_euler_angles(TAU/3.1,0.0,TAU/2.8);

			let rot_mat_quat = SeaMotorPGA3D::<f64>::from(roll_mat);
			let rot_quat = {
				SeaMultivectorPGA3D::new_rotation_motor(0.0,0.0, TAU/2.8)
				*
				SeaMultivectorPGA3D::new_rotation_motor(TAU/3.1,0.0, 0.0)
			}; 

			assert_approx_eq_mv(rot_mat_quat.into(), rot_quat.into());

			let roll_mat = Matrix4::from_euler_angles(-TAU/3.1,0.0,TAU/2.8);

			let rot_mat_quat = SeaMotorPGA3D::<f64>::from(roll_mat);
			let rot_quat = {
				SeaMultivectorPGA3D::new_rotation_motor(0.0,0.0, TAU/2.8)
				*
				SeaMultivectorPGA3D::new_rotation_motor(-TAU/3.1,0.0, 0.0)
			}; 

			assert_approx_eq_mv(rot_mat_quat.into(), rot_quat.into());

			let roll_mat = Matrix4::from_euler_angles(TAU/3.1,0.0,-TAU/2.8);

			let rot_mat_quat = SeaMotorPGA3D::<f64>::from(roll_mat);
			let rot_quat = {
				SeaMultivectorPGA3D::new_rotation_motor(0.0,0.0, -TAU/2.8)
				*
				SeaMultivectorPGA3D::new_rotation_motor(TAU/3.1,0.0, 0.0)
			}; 

			assert_approx_eq_mv(rot_mat_quat.into(), rot_quat.into());

			let roll_mat = Matrix4::from_euler_angles(-TAU/3.1,0.0,-TAU/2.8);

			let rot_mat_quat = SeaMotorPGA3D::<f64>::from(roll_mat);
			let rot_quat = {
				SeaMultivectorPGA3D::new_rotation_motor(0.0,0.0, -TAU/2.8)
				*
				SeaMultivectorPGA3D::new_rotation_motor(-TAU/3.1,0.0, 0.0)
			}; 

			assert_approx_eq_mv(rot_mat_quat.into(), rot_quat.into());
		}

		#[test]
		fn pitch_yaw() {
			let roll_mat = Matrix4::from_euler_angles(TAU/3.1,TAU/4.9,0.0);

			let rot_mat_quat = SeaMotorPGA3D::<f64>::from(roll_mat);
			let rot_quat = {
				SeaMultivectorPGA3D::new_rotation_motor(0.0,TAU/4.9, 0.0)
				*
				SeaMultivectorPGA3D::new_rotation_motor(TAU/3.1,0.0, 0.0)
			}; 

			assert_approx_eq_mv(rot_mat_quat.into(), rot_quat.into());

			let roll_mat = Matrix4::from_euler_angles(-TAU/3.1,TAU/4.9,0.0);

			let rot_mat_quat = SeaMotorPGA3D::<f64>::from(roll_mat);
			let rot_quat = {
				SeaMultivectorPGA3D::new_rotation_motor(0.0,TAU/4.9, 0.0)
				*
				SeaMultivectorPGA3D::new_rotation_motor(-TAU/3.1,0.0, 0.0)
			}; 

			assert_approx_eq_mv(rot_mat_quat.into(), rot_quat.into());

			let roll_mat = Matrix4::from_euler_angles(TAU/3.1,-TAU/4.9,0.0);

			let rot_mat_quat = SeaMotorPGA3D::<f64>::from(roll_mat);
			let rot_quat = {
				SeaMultivectorPGA3D::new_rotation_motor(0.0,-TAU/4.9, 0.0)
				*
				SeaMultivectorPGA3D::new_rotation_motor(TAU/3.1,0.0, 0.0)
			}; 

			assert_approx_eq_mv(rot_mat_quat.into(), rot_quat.into());

			let roll_mat = Matrix4::from_euler_angles(-TAU/3.1,-TAU/4.9,0.0);

			let rot_mat_quat = SeaMotorPGA3D::<f64>::from(roll_mat);
			let rot_quat = {
				SeaMultivectorPGA3D::new_rotation_motor(0.0,-TAU/4.9, 0.0)
				*
				SeaMultivectorPGA3D::new_rotation_motor(-TAU/3.1,0.0, 0.0)
			}; 

			assert_approx_eq_mv(rot_mat_quat.into(), rot_quat.into());
		}

		#[test]
		fn roll_pitch_yaw() {
			let roll_mat = Matrix4::from_euler_angles(TAU/4.0,TAU/4.0, TAU/4.0);

			let rot_mat_quat = SeaMotorPGA3D::<f64>::from(roll_mat);

			let mut rot_quat = {
				SeaMultivectorPGA3D::new_rotation_motor(0.0,0.0, TAU/4.0)
				*
				SeaMultivectorPGA3D::new_rotation_motor(0.0,TAU/4.0, 0.0)
				*
				SeaMultivectorPGA3D::new_rotation_motor(TAU/4.0,0.0, 0.0)
			}; 

			rot_quat.quaternion_angle_positivize();

			assert_approx_eq_mv(rot_mat_quat.into(), rot_quat.into());

			let roll_mat = Matrix4::from_euler_angles(-TAU/3.1,TAU/4.9,TAU/2.8);

			let rot_mat_quat = SeaMotorPGA3D::<f64>::from(roll_mat);
			let mut rot_quat = {
				SeaMultivectorPGA3D::new_rotation_motor(0.0,0.0, TAU/2.8)
				*
				SeaMultivectorPGA3D::new_rotation_motor(0.0,TAU/4.9, 0.0)
				*
				SeaMultivectorPGA3D::new_rotation_motor(-TAU/3.1,0.0, 0.0)
			}; 

			rot_quat.quaternion_angle_positivize();

			assert_approx_eq_mv(rot_mat_quat.into(), rot_quat.into());

			let roll_mat = Matrix4::from_euler_angles(TAU/3.1,-TAU/4.9,TAU/2.8);

			let rot_mat_quat = SeaMotorPGA3D::<f64>::from(roll_mat);
			let mut rot_quat = {
				SeaMultivectorPGA3D::new_rotation_motor(0.0,0.0, TAU/2.8)
				*
				SeaMultivectorPGA3D::new_rotation_motor(0.0,-TAU/4.9, 0.0)
				*
				SeaMultivectorPGA3D::new_rotation_motor(TAU/3.1,0.0, 0.0)
			};

			rot_quat.quaternion_angle_positivize();

			assert_approx_eq_mv(rot_mat_quat.into(), rot_quat.into());

			let roll_mat = Matrix4::from_euler_angles(-TAU/3.1,-TAU/4.9,TAU/2.8);

			let rot_mat_quat = SeaMotorPGA3D::<f64>::from(roll_mat);
			let mut rot_quat = {
				SeaMultivectorPGA3D::new_rotation_motor(0.0,0.0, TAU/2.8)
				*
				SeaMultivectorPGA3D::new_rotation_motor(0.0,-TAU/4.9, 0.0)
				*
				SeaMultivectorPGA3D::new_rotation_motor(-TAU/3.1,0.0, 0.0)
			}; 
			
			rot_quat.quaternion_angle_positivize();

			assert_approx_eq_mv(rot_mat_quat.into(), rot_quat.into());

			let roll_mat = Matrix4::from_euler_angles(TAU/3.1,TAU/4.9,-TAU/2.8);

			let rot_mat_quat = SeaMotorPGA3D::<f64>::from(roll_mat);
			let mut rot_quat = {
				SeaMultivectorPGA3D::new_rotation_motor(0.0,0.0, -TAU/2.8)
				*
				SeaMultivectorPGA3D::new_rotation_motor(0.0,TAU/4.9, 0.0)
				*
				SeaMultivectorPGA3D::new_rotation_motor(TAU/3.1,0.0, 0.0)
			}; 
			
			rot_quat.quaternion_angle_positivize();

			assert_approx_eq_mv(rot_mat_quat.into(), rot_quat.into());

			let roll_mat = Matrix4::from_euler_angles(-TAU/3.1,TAU/4.9,-TAU/2.8);

			let rot_mat_quat = SeaMotorPGA3D::<f64>::from(roll_mat);
			let mut rot_quat = {
				SeaMultivectorPGA3D::new_rotation_motor(0.0,0.0, -TAU/2.8)
				*
				SeaMultivectorPGA3D::new_rotation_motor(0.0,TAU/4.9, 0.0)
				*
				SeaMultivectorPGA3D::new_rotation_motor(-TAU/3.1,0.0, 0.0)
			}; 
			
			rot_quat.quaternion_angle_positivize();

			assert_approx_eq_mv(rot_mat_quat.into(), rot_quat.into());

			let roll_mat = Matrix4::from_euler_angles(TAU/3.1,-TAU/4.9,-TAU/2.8);

			let rot_mat_quat = SeaMotorPGA3D::<f64>::from(roll_mat);
			let mut rot_quat = {
				SeaMultivectorPGA3D::new_rotation_motor(0.0,0.0, -TAU/2.8)
				*
				SeaMultivectorPGA3D::new_rotation_motor(0.0,-TAU/4.9, 0.0)
				*
				SeaMultivectorPGA3D::new_rotation_motor(TAU/3.1,0.0, 0.0)
			}; 
			
			rot_quat.quaternion_angle_positivize();

			assert_approx_eq_mv(rot_mat_quat.into(), rot_quat.into());

			let roll_mat = Matrix4::from_euler_angles(-TAU/3.1,-TAU/4.9,-TAU/2.8);

			let rot_mat_quat = SeaMotorPGA3D::<f64>::from(roll_mat);
			let mut rot_quat = {
				SeaMultivectorPGA3D::new_rotation_motor(0.0,0.0, -TAU/2.8)
				*
				SeaMultivectorPGA3D::new_rotation_motor(0.0,-TAU/4.9, 0.0)
				*
				SeaMultivectorPGA3D::new_rotation_motor(-TAU/3.1,0.0, 0.0)
			}; 
			
			rot_quat.quaternion_angle_positivize();

			assert_approx_eq_mv(rot_mat_quat.into(), rot_quat.into());
		}

	}
	
	#[test]
	fn log() {
		let l = SeaMultivectorPGA3D::new_ideal_line(1.0,TAU/4.0,10.0);
		let r = SeaMultivectorPGA3D::new_motor(&l, 2.0);

		let mut lout = r.log();

		lout.simplify();

		assert_approx_eq_mv(l, lout);

		let l = SeaMultivectorPGA3D::new_line(0.0,TAU/4.0,0.0);
		let r = SeaMultivectorPGA3D::new_motor(&l, 2.0);

		let mut lout = r.log();

		lout.simplify();

		assert_approx_eq_mv(l, lout);

		let mut l = SeaMultivectorPGA3D::new_line(TAU/8.0,TAU/4.0,TAU/3.2);
		l.scalar_normalize();
		let r = SeaMultivectorPGA3D::new_motor(&l, 2.0);

		let mut lout = r.log();

		lout.simplify();

		assert_approx_eq_mv(l, lout);
	}
}