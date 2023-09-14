use super::*;

impl<T:SeaCoefficient> SandwichProduct<Self> for SeaMotorPGA3D<T> {
	type Output = Self;

	fn sandwich_product(&self, rhs: Self) -> Self::Output {
		let a_mv: &SeaCoreMultivector4D<T> = (&**self).into();
		let b_mv: SeaCoreMultivector4D<T> = rhs.0.into();
		Self(
			SeaMultivectorPGA3D::from(
				a_mv.sandwich_product(&b_mv, &METRIC)
			)
		)
	}
}
impl<T:SeaCoefficient> SandwichProduct<&Self> for SeaMotorPGA3D<T> {
	type Output = Self;

	fn sandwich_product(&self, rhs: &Self) -> Self::Output {
		let a_mv: &SeaCoreMultivector4D<T> = (&**self).into();
		let b_mv: &SeaCoreMultivector4D<T> = (&**rhs).into();
		Self(
			SeaMultivectorPGA3D::from(
				a_mv.sandwich_product(&b_mv, &METRIC)
			)
		)
	}
}
impl<T:SeaCoefficient> Div<Self> for SeaMotorPGA3D<T> {
	type Output = SeaMotorPGA3D<T>;
	fn div(self, rhs: Self) -> Self::Output {
		self.sandwich_product(rhs)
	}
}
impl<T:SeaCoefficient> Div<&Self> for SeaMotorPGA3D<T> {
	type Output = SeaMotorPGA3D<T>;
	fn div(self, rhs: &Self) -> Self::Output {
		self.sandwich_product(rhs)
	}
}
impl<T:SeaCoefficient> Div<SeaMotorPGA3D<T>> for &SeaMotorPGA3D<T> {
	type Output = SeaMotorPGA3D<T>;
	fn div(self, rhs: SeaMotorPGA3D<T>) -> Self::Output {
		self.sandwich_product(rhs)
	}
}
impl<T:SeaCoefficient> Div<Self> for &SeaMotorPGA3D<T> {
	type Output = SeaMotorPGA3D<T>;
	fn div(self, rhs: Self) -> Self::Output {
		self.sandwich_product(rhs)
	}
}

// Multivector
impl<T:SeaCoefficient> SandwichProduct<SeaMultivectorPGA3D<T>> for SeaMotorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;

	fn sandwich_product(&self, rhs: SeaMultivectorPGA3D<T>) -> Self::Output {
		let a_mv: &SeaCoreMultivector4D<T> = (&**self).into();
		let b_mv: SeaCoreMultivector4D<T> = rhs.into();
		let out_mv = a_mv.sandwich_product(&b_mv, &METRIC);
		SeaMultivectorPGA3D::from(out_mv)
	}
}
impl<T:SeaCoefficient> SandwichProduct<&SeaMultivectorPGA3D<T>> for SeaMotorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;

	fn sandwich_product(&self, rhs: &SeaMultivectorPGA3D<T>) -> Self::Output {
		let a_mv: &SeaCoreMultivector4D<T> = (&**self).into();
		let b_mv: &SeaCoreMultivector4D<T> = (&*rhs).into();
		let out_mv = a_mv.sandwich_product(&b_mv, &METRIC);
		SeaMultivectorPGA3D::from(out_mv)
	}
}
impl<T:SeaCoefficient> Div<SeaMultivectorPGA3D<T>> for SeaMotorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn div(self, rhs: SeaMultivectorPGA3D<T>) -> Self::Output {
		self.sandwich_product(rhs)
	}
}
impl<T:SeaCoefficient> Div<&SeaMultivectorPGA3D<T>> for SeaMotorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn div(self, rhs: &SeaMultivectorPGA3D<T>) -> Self::Output {
		self.sandwich_product(rhs)
	}
}
impl<T:SeaCoefficient> Div<SeaMultivectorPGA3D<T>> for &SeaMotorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn div(self, rhs: SeaMultivectorPGA3D<T>) -> Self::Output {
		self.sandwich_product(rhs)
	}
}
impl<T:SeaCoefficient> Div<&SeaMultivectorPGA3D<T>> for &SeaMotorPGA3D<T> {
	type Output = SeaMultivectorPGA3D<T>;
	fn div(self, rhs: &SeaMultivectorPGA3D<T>) -> Self::Output {
		self.sandwich_product(rhs)
	}
}

#[cfg(test)]
mod tests {
	use std::f64::consts::TAU;

	use super::super::super::super::tests::assert_approx_eq_vec3;

	use super::*;

	#[test]
	fn rotation() {
		let translation_dir = SeaMultivectorPGA3D::new_ideal_line(0.0,0.0,0.0);
		let camera_pos = SeaMultivectorPGA3D::new_motor(&translation_dir, -3.0);

		let rotation_axis = SeaMultivectorPGA3D::new_line(0.0,0.0,1.0);
		let camera_rot = SeaMultivectorPGA3D::new_motor(&rotation_axis, -TAU/4.0);

		let camera_trans = camera_rot * camera_pos;
		
		let scene_point = SeaMultivectorPGA3D::new_point(0.0,3.0,2.0);
		let scene_point_trans = camera_trans / scene_point;

		assert_approx_eq_vec3(
			scene_point_trans.point_to_euclidean(),
			Vector3::new(-3.0, 0.0, 2.0)
		);
	}

	#[test]
	fn translation() {
		let translation_dir = SeaMultivectorPGA3D::new_ideal_line(0.0,1.0,1.0);
		let camera_pos = SeaMultivectorPGA3D::new_motor(&translation_dir, -3.0);

		let rotation_axis = SeaMultivectorPGA3D::new_line(0.0,0.0,0.0);
		let camera_rot = SeaMultivectorPGA3D::new_motor(&rotation_axis, -TAU/4.0);

		let camera_trans = camera_rot * camera_pos;
		
		let scene_point = SeaMultivectorPGA3D::new_point(0.0,3.0,2.0);
		let scene_point_trans = camera_trans / scene_point;

		assert_approx_eq_vec3(
			scene_point_trans.point_to_euclidean(),
			Vector3::new(0.0, 0.0, -1.0) 
		);
	}

	#[test]
	fn combined() {
		let translation_dir = SeaMultivectorPGA3D::new_ideal_line(0.0,1.0,1.0);
		let camera_pos = SeaMultivectorPGA3D::new_motor(&translation_dir, -3.0);

		let rotation_axis = SeaMultivectorPGA3D::new_line(0.0,1.0,0.0);
		let camera_rot = SeaMultivectorPGA3D::new_motor(&rotation_axis, -TAU/4.0);

		let camera_trans = camera_rot * camera_pos;
		
		let scene_point = SeaMultivectorPGA3D::new_point(0.0,3.0,2.0);
		let scene_point_trans = camera_trans / scene_point;

		assert_approx_eq_vec3(
			scene_point_trans.point_to_euclidean(),
			Vector3::new(-1.0, 0.0, 0.0) 
		);
	}
}