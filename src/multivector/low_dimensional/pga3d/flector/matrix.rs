use super::*;

//matrix
impl<T:SeaCoefficient + 'static> SeaFlectorPGA3D<T> {
	pub fn to_homogeneous_transformation_matrices(&self) -> (Matrix4<T>, Matrix4<T>) {
		let hx = self.basis_blade_coefficient(&SeaBasisBladeIndexLow::E1);
		let hy = self.basis_blade_coefficient(&SeaBasisBladeIndexLow::E2);
		let hz = self.basis_blade_coefficient(&SeaBasisBladeIndexLow::E3);
		let hw = self.basis_blade_coefficient(&SeaBasisBladeIndexLow::E4);
		let sw = self.basis_blade_coefficient(&SeaBasisBladeIndexLow::E123);
		let sz = -self.basis_blade_coefficient(&SeaBasisBladeIndexLow::E124);
		let sy = self.basis_blade_coefficient(&SeaBasisBladeIndexLow::E134);
		let sx = -self.basis_blade_coefficient(&SeaBasisBladeIndexLow::E234);

		let zero = T::zero();
		let one = T::one();
		let two = T::two();

		let a: Matrix4<T> = Matrix4::new(
			two*(hy*hy + hz*hz) - one, -two*hx*hy, -two*hx*hz, two*(sx*sw - hx*hw),
			-two*hx*hy, two*(hx*hx + hz*hz) - one, -two*hy*hz, two*(sy*sw - hy*hw),
			-two*hx*hz, -two*hy*hz, two*(hx*hx + hy*hy) - one, two*(sz*sw - hz*hw),
			zero, zero, zero, one
		);
		
		let b: Matrix4<T> = Matrix4::new(
			zero, two*hz*sw, -two*hy*sw, two*(hy*sz - hz*sy),
			-two*hz*sw, zero, two*hx*sw, two*(hz*sx - hx*sz),
			two*hy*sw, -two*hx*sw, zero, two*(hx*sy - hy*sx),
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
	fn basic_plane() {
		let p = SeaMultivectorPGA3D::new_point(1.0, 2.0, 3.0);

		let f = SeaMultivectorPGA3D::new_plane(0.0,1.0,0.0,0.0);

		let(m, m_inv) = f.to_homogeneous_transformation_matrices();

		let p_mat_homo = m*p.point_to_homogeneous();

		let p_trans = &f / &p;
		let p_trans_homo = p_trans.point_to_homogeneous();

		assert_approx_eq_vec4(
			p_mat_homo,
			p_trans_homo
		);
		assert_eq!(
			m,
			Matrix4::new_nonuniform_scaling(&Vector3::new(1.0,-1.0,1.0))
		);
		assert_eq!(
			m * m_inv,
			Matrix4::identity()
		);
		assert_eq!(
			m * m_inv,
			m_inv * m
		);
	}

	#[test]
	fn different_plane() {
		let p = SeaMultivectorPGA3D::new_point(1.0, 2.0, 3.0);

		let mut f = SeaMultivectorPGA3D::new_plane(1.2,3.2342,5.534, 2.775);
		f.scalar_normalize();

		let(m, m_inv) = f.to_homogeneous_transformation_matrices();

		let p_mat_homo = m*p.point_to_homogeneous();

		let p_trans = &f / &p;
		let p_trans_homo = p_trans.point_to_homogeneous();

		assert_approx_eq_vec4(
			p_mat_homo,
			p_trans_homo
		);
		assert_approx_eq_mat4(
			m * m_inv,
			Matrix4::identity()
		);
		assert_eq!(
			m * m_inv,
			m_inv * m
		);
	}

	#[test]
	fn basic_point() {
		let p = SeaMultivectorPGA3D::new_point(1.0, 2.0, 3.0);

		let f = SeaMultivectorPGA3D::new_point(0.0, 0.0, 0.0);

		let(m, m_inv) = f.to_homogeneous_transformation_matrices();

		let p_mat_homo = m*p.point_to_homogeneous();

		let p_trans = &f / &p;
		let p_trans_homo = p_trans.point_to_homogeneous();

		assert_approx_eq_vec4(
			p_mat_homo,
			Vector4::new(-1.0, -2.0, -3.0, 1.0)
		);

		assert_approx_eq_vec4(
			p_mat_homo,
			p_trans_homo
		);
		assert_approx_eq_mat4(
			m * m_inv,
			Matrix4::identity()
		);
		assert_eq!(
			m * m_inv,
			m_inv * m
		);
	}

	#[test]
	fn different_points() {
		let p = SeaMultivectorPGA3D::new_point(1.0, 2.0, 3.0);

		let f = SeaMultivectorPGA3D::new_point(1.2, 3.4, 7.5);

		let(m, m_inv) = f.to_homogeneous_transformation_matrices();

		let p_mat_homo = m*p.point_to_homogeneous();

		let p_trans = &f / &p;
		let p_trans_homo = p_trans.point_to_homogeneous();

		assert_approx_eq_vec4(
			p_mat_homo,
			Vector4::new(1.4, 4.8, 12.0, 1.0)
		);
		assert_approx_eq_vec4(
			p_mat_homo,
			p_trans_homo
		);
		assert_approx_eq_mat4(
			m * m_inv,
			Matrix4::identity()
		);
		assert_eq!(
			m * m_inv,
			m_inv * m
		);

		let p = SeaMultivectorPGA3D::new_point(1.0, 2.0, 3.0);

		let f = SeaMultivectorPGA3D::new_point(0.4, 1.1, 3.1);

		let(m, m_inv) = f.to_homogeneous_transformation_matrices();

		let p_mat_homo = m*p.point_to_homogeneous();

		let p_trans = &f / &p;
		let p_trans_homo = p_trans.point_to_homogeneous();

		assert_approx_eq_vec4(
			p_mat_homo,
			Vector4::new(-0.2, 0.2, 3.2, 1.0)
		);
		assert_approx_eq_vec4(
			p_mat_homo,
			p_trans_homo
		);
		assert_approx_eq_mat4(
			m * m_inv,
			Matrix4::identity()
		);
		assert_eq!(
			m * m_inv,
			m_inv * m
		);
	}
}