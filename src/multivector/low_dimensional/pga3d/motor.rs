use std::f64::consts::TAU;

use super::*;

mod scalar_operations;
mod linear_operations;
mod non_linear_operations;
mod non_geometric_linear_operations;

#[derive(Debug, PartialEq)]
pub struct SeaMotorPGA3D<T:SeaCoefficient>(pub SeaMultivectorPGA3D<T>);
impl<T:SeaCoefficient> SeaMotorPGA3D<T> {
	fn log(&self) -> SeaMultivectorPGA3D<T> {
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

	fn to_homogeneous_transformation_matrices(&self) -> (Matrix4<T>, Matrix4<T>) {
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
impl<T:SeaCoefficient> MultivectorMut<T> for SeaMotorPGA3D<T> {}

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

	mod to_matrix{
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