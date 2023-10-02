use super::*;

// vector
impl<T:SeaCoefficient> SeaMultivectorPGA3D<T>{
	pub fn point_to_homogeneous(&self) -> Vector4<T> {
		let x = -self.0.basis_blade_coefficient(&SeaBasisBladeIndexLow::E234);
		let y =  self.0.basis_blade_coefficient(&SeaBasisBladeIndexLow::E134);
		let z = -self.0.basis_blade_coefficient(&SeaBasisBladeIndexLow::E124);
		let w =  self.0.basis_blade_coefficient(&SeaBasisBladeIndexLow::E123);

		Vector4::new(x,y,z,w)
	}

	pub fn point_to_euclidean(&self) -> Vector3<T> {
		let x = -self.0.basis_blade_coefficient(&SeaBasisBladeIndexLow::E234);
		let y =  self.0.basis_blade_coefficient(&SeaBasisBladeIndexLow::E134);
		let z = -self.0.basis_blade_coefficient(&SeaBasisBladeIndexLow::E124);
		let w =  self.0.basis_blade_coefficient(&SeaBasisBladeIndexLow::E123);

		Vector3::new(x/w, y/w, z/w)
	}
}

#[cfg(test)]
pub mod tests {
	use std::f64::consts::TAU;

	use super::super::tests::assert_approx_eq_vec3;

	use super::*;

	#[test]
	fn point_to_vec() {
		let p1 = SeaMultivectorPGA3D::new_point(1.0,2.0,3.0);
		let v = p1.point_to_euclidean();
		assert_approx_eq_vec3(
			v,
			Vector3::new(1.0, 2.0, 3.0)
		);
	}
}