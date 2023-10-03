use super::*;

impl<T:SeaCoefficient + 'static> From<Matrix3<T>> for SeaFlectorPGA3D<T> {
	fn from(plane_mat: Matrix3<T>) -> Self {
		Self::from(&plane_mat)
	}		
}
impl<T:SeaCoefficient + 'static> From<&Matrix3<T>> for SeaFlectorPGA3D<T> {
	fn from(plane_mat: &Matrix3<T>) -> Self {
		let zero = T::zero();
		let one = T::one();
		let two = T::two();
		let four = T::four();

		let mut line = [T::zero(); 3];

		let sum_sqares = ( ( plane_mat[(0, 0)] +  plane_mat[(1, 1)] + plane_mat[(2, 2)] - one)/four) + one;
			
		let sum_two_x_squares = ( plane_mat[(1, 1)] +  plane_mat[(2,2)] + two )/two;
		let sum_two_y_squares = ( plane_mat[(0, 0)] +  plane_mat[(2,2)] + two )/two;
		let two_z_squares = ( plane_mat[(0, 0)] +  plane_mat[(1,1)] + two )/two;

		let x_squared = sum_two_x_squares - sum_sqares;
		let y_squared = sum_two_y_squares - sum_sqares;
		let z_squared = two_z_squares - sum_sqares;

		line[0] = x_squared.sqrt(); 
		line[1] = y_squared.sqrt(); 
		line[2] = z_squared.sqrt();

		let x_y = -(plane_mat[(1,0)] + plane_mat[(0,1)])/four;
		let x_z = -(plane_mat[(2,0)] + plane_mat[(0,2)])/four;

		if x_y < zero {
			line[1] = -line[1];
		}; 
		if x_z < zero {
			line[2] = -line[2];
		}; 

		Self::from_origin_plane(-line[0], -line[1], -line[2])
	}
}
impl<T:SeaCoefficient + 'static> From<Matrix4<T>> for SeaFlectorPGA3D<T> {
	fn from(value: Matrix4<T>) -> Self {
		Self::from(&value)
	}
}
impl<T:SeaCoefficient + 'static> From<&Matrix4<T>> for SeaFlectorPGA3D<T> {
	fn from(flec_mat: &Matrix4<T>) -> Self {
		let two = T::two();
		let four = T::four();

		let plane_mat = Matrix3x4::identity() * flec_mat * Matrix4x3::identity();
		let plane_flec = Self::from(&plane_mat);

		let hx = plane_flec.basis_blade_coefficient(&SeaBasisBladeIndexLow::E1);
		let hy = plane_flec.basis_blade_coefficient(&SeaBasisBladeIndexLow::E2);
		let hz = plane_flec.basis_blade_coefficient(&SeaBasisBladeIndexLow::E3);
		let hw = plane_flec.basis_blade_coefficient(&SeaBasisBladeIndexLow::E4);

		let sum = hx + hy + hz;
		
		let sw = if sum.is_zero() {
			T::one()
		}
		else {
			((plane_mat[(0,1)] - plane_mat[(1,0)]) +
			(plane_mat[(2,0)] - plane_mat[(0,2)]) +
			(plane_mat[(1,2)] - plane_mat[(2,1)]))/(four*sum)
		};

		if sw.is_zero() {
			let hw = -(flec_mat[(0, 3)] + flec_mat[(1, 3)] + flec_mat[(2, 3)])/(two*(hx + hy + hz));
			
			if hw.is_sign_negative() {
				return Self::from_plane(-hx, -hy, -hz, -hw);
			}
			return Self::from_plane(hx, hy, hz, hw);
		}

		let d1 = flec_mat[(0, 3)]/two + hx*hw;
		let d2 = flec_mat[(1,3)]/two + hy*hw;
		let d3 = flec_mat[(2,3)]/two + hz*hw;

		let p1 = Self::from_plane(sw, hz, hy,  -d1);
		let p2 = Self::from_plane(hz, sw, hx,  -d2);
		let p3 = Self::from_plane(hy, hx, sw,  -d3);

		let point_flec = p1 ^ p2 ^ p3;

		plane_flec + point_flec
	}
}

#[cfg(test)]
mod tests {
	use std::f64::consts::TAU;
	
	use funty::Floating;
use multivector::tests::assert_approx_eq_vec3;

	use super::super::tests::assert_approx_eq_vec4;
	use super::super::tests::assert_approx_eq_mat4;
	use super::multivector::tests::assert_approx_eq_mv;

	use super::*;

	#[test]
	fn from_matrix3() {
		let f = SeaMultivectorPGA3D::new_plane(-1.0,2.0,3.0,0.0);

		let(m, _) = f.to_homogeneous_transformation_matrices();

		let pm =  Matrix3x4::identity() * &m * Matrix4x3::identity();
		let fout = SeaFlectorPGA3D::from(pm);

		assert_approx_eq_mv(
			f.into(), 
			fout.into()
		);
	}

	#[test]
	fn from_matrix4(){
		let p = SeaMultivectorPGA3D::new_point(-7.0,5.0,2.0);
		
		let(m, _) = p.to_homogeneous_transformation_matrices();

		let pfout = SeaFlectorPGA3D::from(m);
		
		assert_approx_eq_mv(
			p.into(), 
			pfout.into()
		);


		let p = SeaMultivectorPGA3D::new_plane(-1.0,2.0,3.0,1e-15);

		let(m, _) = p.to_homogeneous_transformation_matrices();

		let pfout = SeaFlectorPGA3D::from(m);
		
		assert_approx_eq_mv(
			p.into(), 
			pfout.into()
		);

		let p = SeaMultivectorPGA3D::new_plane(-3.5,12.3,-2.5,1e-15 );
		
		let(m, _) = p.to_homogeneous_transformation_matrices();
		
		let pfout = SeaFlectorPGA3D::from(m);
		
		assert_approx_eq_mv(
			p.into(), 
			pfout.into()
		);
	}

	#[test]
	fn test_matrix_transform() {
		let w = SeaMultivectorPGA3D::new_point(1.0,2.0,3.0);
		let mut p = SeaMultivectorPGA3D::new_plane(-1.0,2.0,3.0,0.0);
		p.scalar_normalize();

		let wm = w.point_to_homogeneous();
		let(m, _) = p.to_homogeneous_transformation_matrices();
		
		let out_mv = p / w;
		let out_mat = m * wm;
		
		assert_approx_eq_vec4(
			out_mat, 
			out_mv.point_to_homogeneous()
		);

		let w = SeaMultivectorPGA3D::new_point(1.0,2.0,3.0);
		let p = SeaMultivectorPGA3D::new_point(-7.0,5.0,2.0);
		let f = SeaMultivectorPGA3D::new_plane(-3.5,12.3,-2.5,1e-15 );
		let mut pw = f + p;
		pw.scalar_normalize();

		let wm = w.point_to_homogeneous();
		let(m, _) = pw.to_homogeneous_transformation_matrices();
		
		let out_mv = pw / w;
		let out_mat = m * wm;
		
		assert_approx_eq_vec4(
			out_mat,
			out_mv.point_to_homogeneous()
		);
	}
}