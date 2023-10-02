use super::*;

//matrix
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
}


#[cfg(test)]
mod tests {
	use std::f64::consts::TAU;
	
	use super::super::tests::assert_approx_eq_vec4;
	use super::super::tests::assert_approx_eq_mat4;
	use super::multivector::tests::assert_approx_eq_mv;

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