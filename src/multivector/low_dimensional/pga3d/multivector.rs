use super::*;

mod scalar_operations;
mod linear_operations;
mod non_linear_operations;
mod non_geometric_linear_operations;

mod matrix;

mod iterators;
mod from;
mod into;

#[derive(Debug, PartialEq, Clone)]
pub struct SeaMultivectorPGA3D<T:SeaCoefficient>(SeaCoreMultivector4D<T>);
impl<T:SeaCoefficient> SeaMultivectorPGA3D<T>{
	pub fn new_lifted_scalar(s:T) -> Self {
		let mv = SeaCoreMultivector4D::<T>::new(
			SeaBasisBladeSet::SCALAR,
			vec![s]
		);
		Self(mv)
	}

	pub fn new_point(x:T, y:T, z:T) -> Self {
		let mut mv = SeaCoreMultivector4D::<T>::new(
			SeaBasisBladeSet::grade(&SeaGrade(3)),
			vec![T::one(),-z,y,-x]
		);
		mv.remove_zeros();
		Self(mv)
	}

	pub fn new_plane(x:T, y:T, z:T, w:T) -> Self {
		let mut mv = SeaCoreMultivector4D::<T>::new(
			SeaBasisBladeSet::grade(&SeaGrade(1)),
			vec![x,y,z,w]
		);
		mv.remove_zeros();
		Self(mv)
	}

	pub fn new_line(x:T, y:T, z:T) -> Self {
		let mut mv = SeaCoreMultivector4D::<T>::new(
			LINE_BASIS_BLADES,
			vec![-z,y,-x]
		);
		mv.remove_zeros();
		Self(mv)
	}

	pub fn new_ideal_line(x:T, y:T, z:T) -> Self {
		let mut mv = SeaCoreMultivector4D::<T>::new(
			IDEAL_LINE_BASIS_BLADES,
			vec![x,y,z]
		);
		mv.remove_zeros();
		Self(mv)
	}

	pub fn new_motor(line:&Self, value:T) -> SeaMotorPGA3D<T> {
		let two = T::two();
		let inner = (line*value/two).exp();
		SeaMotorPGA3D(inner)
	}

	pub fn identity() -> Self {
		Self::from(T::one())
	}

	pub fn is_zero(&self) -> bool{
		self.0.is_zero()
	}

	pub fn simplify(&mut self) {
		self.0.remove_approx_zeros(T::from(1e-15).unwrap());
		self.0.round_to_decimals(14);
	}

	pub fn simplify_epsilon(&mut self, epsilon: T) {
		self.0.remove_approx_zeros(epsilon);
		self.0.round_to_decimals(14);
	}

	pub fn raw(&self) -> (u16, [T;16]) {
		let mut out_mat = [T::zero(); 16];
		for (i,o) in out_mat.iter_mut().enumerate() {
			*o = self.0.basis_blade_coefficient(&SeaBasisBladeIndexLow(i as u8))
		}
		let out_bbb = self.0.basis_blade_set;
		(out_bbb, out_mat)
	}

}

impl<T:SeaCoefficient> Multivector<T> for SeaMultivectorPGA3D<T> {
	type MultivectorIndexType = SeaBasisBladeIndexLow;
	type GeneralMultivectorOutput = Self;

	fn basis_blade_coefficient(&self, index: &Self::MultivectorIndexType) -> T {
		self.0.basis_blade_coefficient(index)
	}

	fn scalar(&self) -> T {
		self.basis_blade_coefficient(&SeaBasisBladeIndexLow::SCALAR)
	}

	fn pseudoscalar(&self) -> Self {
		self.grade(&SeaGrade(Self::METRIC.dimensions))
	}

	fn grade(&self, grade: &SeaGrade) -> Self {
		Self(self.0.grade(grade))
	}

	fn exp(&self) -> Self {
		Self(self.0.exp(&Self::METRIC))
	}
}
impl<T:SeaCoefficient> MultivectorMut<T> for SeaMultivectorPGA3D<T> {}

impl<T: SeaCoefficient> MultivectorConstMetric<T> for SeaMultivectorPGA3D<T> {
	const METRIC: SeaMetric = METRIC_PGA3D;
}

impl<T:SeaCoefficient> Default for SeaMultivectorPGA3D<T> {
	fn default() -> Self {
		Self(SeaCoreMultivector4D::<T>::default())
	}
}


#[cfg(test)]
pub mod tests {
	use std::f64::consts::TAU;

	pub use super::super::tests::assert_approx_eq_vec3;

	use super::*;

	pub fn assert_approx_eq_mv<T: SeaCoefficient>(a: SeaMultivectorPGA3D<T>, b: SeaMultivectorPGA3D<T>) {
		assert_eq!(a.0.basis_blade_set, b.0.basis_blade_set, "{:?} != {:?}", a, b);
		for (ca, cb) in a.0.coefficients.iter().zip(&b.0.coefficients) {
			assert!((*ca - *cb).abs() < T::from(1e-14).unwrap(), "{:?} != {:?}", a, b)
		}
	}

    #[test]
    fn create_plane_from_points() {
		let p1 = SeaMultivectorPGA3D::new_point(1.0,0.0,0.0);
		let p2 = SeaMultivectorPGA3D::new_point(0.0,1.0,0.0);
		let p3 = SeaMultivectorPGA3D::new_point(0.0,0.0,0.0);
		let c = p1&p2&p3;
		let out = SeaMultivectorPGA3D::new_plane(0.0,0.0,1.0, 0.0);

		assert_eq!(
			c,
			out
		);

		let p1 = SeaMultivectorPGA3D::new_point(1.0,2.0,3.0);
		let p2 = SeaMultivectorPGA3D::new_point(2.0,1.0,4.0);
		let p3 = SeaMultivectorPGA3D::new_point(2.0,5.0,7.0);
		let c = p1&p2&p3;
		let out = SeaMultivectorPGA3D::new_plane(-7.0,-3.0,4.0, 1.0);

		assert_eq!(
			c,
			out
		);

		let p1 = SeaMultivectorPGA3D::new_point(1.0,2.0,3.0);
		let p2 = SeaMultivectorPGA3D::new_point(2.0,1.0,4.0);
		let p3 = SeaMultivectorPGA3D::new_point(2.0,5.0,7.0);
		let c = p1&p2&p3;
		let out = SeaMultivectorPGA3D::new_plane(-7.0,-3.0,4.0, 1.0);

		assert_eq!(
			c,
			out
		);
	}

	#[test]
	fn project_point_onto_plane() {
		let plane = SeaMultivectorPGA3D::new_plane(1.0,0.0,0.0, 0.0);
		let p = SeaMultivectorPGA3D::new_point(1.0,1.0,1.0);
		let proj = (&plane|p)*plane;
		
		let out = Vector3::new(0.0, 1.0, 1.0);
		assert_approx_eq_vec3(
			proj.point_to_euclidean(),
			out
		);

		let plane = SeaMultivectorPGA3D::new_plane(0.0,1.0,0.0, 0.0);
		let p = SeaMultivectorPGA3D::new_point(1.0,1.0,1.0);
		let proj = (&plane|p)*plane;
		
		let out = Vector3::new(1.0, 0.0, 1.0);
		assert_approx_eq_vec3(
			proj.point_to_euclidean(),
			out
		);

		let plane = SeaMultivectorPGA3D::new_plane(0.0,0.0,1.0, 0.0);
		let p = SeaMultivectorPGA3D::new_point(1.0,1.0,1.0);
		let proj = (&plane|p)*plane;
		
		let out = Vector3::new(1.0, 1.0, 0.0);
		assert_approx_eq_vec3(
			proj.point_to_euclidean(),
			out
		);

		let plane = SeaMultivectorPGA3D::new_plane(-7.0,-3.0,4.0, 1.0);
		let p4 = SeaMultivectorPGA3D::new_point(0.0,0.0,0.0);
		let proj = (&plane|p4)*plane;

		let out = Vector3::new(0.0945945945945946, 0.04054054054054054, -0.05405405405405406);
		assert_approx_eq_vec3(
			proj.point_to_euclidean(),
			out
		);

		let plane = SeaMultivectorPGA3D::new_plane(1.0,0.0,0.0, -1.0);
		let p4 = SeaMultivectorPGA3D::new_point(100.0,0.0,0.0);
		let proj = (&plane|p4)*plane;

		let out = Vector3::new(1.0, 0.0, 0.0);
		assert_approx_eq_vec3(
			proj.point_to_euclidean(),
			out
		);
	}

	#[test]
	fn rotor() {
		let l = SeaMultivectorPGA3D::new_line(5.0,3.0,2.0);
		let r = SeaMultivectorPGA3D::new_motor(&l, TAU/1.6);
		let r_inv = r.inverse();
		let mut approx_one = r * r_inv;
		approx_one.simplify();
		assert_eq!(
			SeaMultivectorPGA3D::from(approx_one),
			SeaMultivectorPGA3D::new_lifted_scalar(1.0),
		);


		let l = SeaMultivectorPGA3D::new_line(0.0,0.0,1.0);
		let r = SeaMultivectorPGA3D::new_motor(&l, TAU/2.0);
		let r_inv = r.inverse();
		let p = SeaMultivectorPGA3D::new_point(1.0, 0.0, 0.0);
		let p_rot = r * p * r_inv;

		let out = Vector3::new(-1.0, 0.0, 0.0);

		assert_approx_eq_vec3(
			p_rot.point_to_euclidean(),
			out
		)



	}

	#[test]
	fn translator() {
		let l = SeaMultivectorPGA3D::new_ideal_line(2.0,3.0,5.0);
		let r = SeaMultivectorPGA3D::new_motor(&l, TAU/3.2);
		let r_inv = r.inverse();
		let mut approx_one = r * r_inv;
		approx_one.simplify();
		assert_eq!(
			SeaMultivectorPGA3D::from(approx_one),
			SeaMultivectorPGA3D::new_lifted_scalar(1.0),
		);


		let l = SeaMultivectorPGA3D::new_ideal_line(1.0,0.0,0.0);
		let r = SeaMultivectorPGA3D::new_motor(&l, 4.0);
		let r_inv = r.inverse();
		let p = SeaMultivectorPGA3D::new_point(0.0, 0.0, 0.0);
		let p_trans = r * p * r_inv;

		let out = Vector3::new(4.0, 0.0, 0.0);

		assert_approx_eq_vec3(
			p_trans.point_to_euclidean(),
			out
		)



	}
}