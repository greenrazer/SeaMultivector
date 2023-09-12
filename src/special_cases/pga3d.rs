use std::ops::{Deref, DerefMut, Shl, Shr, BitOr, BitXor, Mul, Not, BitAnd, Index, AddAssign, Add, Neg, Sub, Div, SubAssign, MulAssign, DivAssign};

use crate::coefficient::SeaCoefficient;
use crate::multivector_product_traits::ScalarProduct;
use crate::multivector_traits::Multivector;
use crate::grade::SeaGrade;
use crate::low_dimensional_multivector::Multivector4D;

use crate::low_dimensional_multivector::SeaBasisBladeSet;
use crate::low_dimensional_multivector::SeaBasisBladeIndexLow;
use crate::multivector_traits::MultivectorMut;

#[derive(Debug, PartialEq)]
pub struct SeaMotorPGA3D<T:SeaCoefficient>(pub MultivectorPGA3D<T>);
impl<T:SeaCoefficient> SeaMotorPGA3D<T> {
	pub fn inversed(&self) -> Self {
		Self(MultivectorPGA3D((**self).0.reversion_conjugated()))
	}

	pub fn into(self) -> MultivectorPGA3D<T> {
		self.0
	}
}
impl<T: SeaCoefficient> Deref for SeaMotorPGA3D<T> {
	type Target = MultivectorPGA3D<T>;
	fn deref(&self) -> &Self::Target {
		&self.0
	}
}
impl<T: SeaCoefficient> DerefMut for SeaMotorPGA3D<T> {
	fn deref_mut(&mut self) -> &mut Self::Target {
		&mut self.0
	}
}

#[derive(Debug, PartialEq, Clone)]
pub struct MultivectorPGA3D<T:SeaCoefficient>(Multivector4D<T>);

const LINE_BASIS_BLADES: u16 = 0x0068;
const IDEAL_LINE_BASIS_BLADES: u16 = 0x1600;

impl<T:SeaCoefficient> MultivectorPGA3D<T>{
	pub fn new() -> Self{
		Self(Multivector4D::new(3, 0, 1))
	}

	pub fn new_lifted_scalar(s:T) -> Self {
		let mut mv = Self::new();
		mv.0.multivector.basis_blade_bits = SeaBasisBladeSet::SCALAR;
		mv.0.multivector.coefficients = vec![s];
		mv
	}

	pub fn new_point(x:T, y:T, z:T) -> Self {
		let mut mv = Self::new();
		mv.0.multivector.basis_blade_bits = SeaBasisBladeSet::grade(&SeaGrade(3));
		mv.0.multivector.coefficients = vec![T::one(),-z,y,-x];
		mv.0.multivector.remove_zeros();
		mv
	}

	pub fn new_plane(x:T, y:T, z:T, w:T) -> Self {
		let mut mv = Self::new();
		mv.0.multivector.basis_blade_bits = SeaBasisBladeSet::grade(&SeaGrade(1));
		mv.0.multivector.coefficients = vec![x,y,z,w];
		mv.0.multivector.remove_zeros();
		mv
	}

	pub fn new_line(x:T, y:T, z:T) -> Self {
		let mut mv = Self::new();
		mv.0.multivector.basis_blade_bits = LINE_BASIS_BLADES;
		mv.0.multivector.coefficients = vec![z,-y,x];
		mv.0.multivector.remove_zeros();
		mv
	}

	pub fn new_ideal_line(x:T, y:T, z:T) -> Self {
		let mut mv = Self::new();
		mv.0.multivector.basis_blade_bits = IDEAL_LINE_BASIS_BLADES;
		mv.0.multivector.coefficients = vec![x,y,z];
		mv.0.multivector.remove_zeros();
		mv
	}

	pub fn new_combined_line(x:T, y:T, z:T, xy:T, xz:T, yz:T) -> Self {
		let mut mv = Self::new();
		mv.0.multivector.basis_blade_bits = SeaBasisBladeSet::grade(&SeaGrade(2));
		mv.0.multivector.coefficients = vec![xy, xz, yz, x,y,z];
		mv.0.multivector.remove_zeros();
		mv
	}

	pub fn new_motor(line:&Self, value:T) -> SeaMotorPGA3D<T> {
		SeaMotorPGA3D(Self((line*value/T::from(2.0).unwrap()).0.exp()))
	}

	pub fn is_zero(&self) -> bool{
		self.0.multivector.basis_blade_bits.is_empty()
	}

	pub fn point_to_vec4(&self) -> [T;4] {
		let mut out = [
			-self.0.basis_blade_coefficient(&SeaBasisBladeIndexLow(14)),
			self.0.basis_blade_coefficient(&SeaBasisBladeIndexLow(13)),
			-self.0.basis_blade_coefficient(&SeaBasisBladeIndexLow(11)),
			self.0.basis_blade_coefficient(&SeaBasisBladeIndexLow(7))
		];

		for o in &mut out {
			if o.abs() < T::from(1e-15).unwrap() {
				*o = T::zero()
			}
		}
		out
	}

	pub fn point_to_vec3(&self) -> [T;3] {
		let sf = self.0.basis_blade_coefficient(&SeaBasisBladeIndexLow(7));
		let mut out = [
			-self.0.basis_blade_coefficient(&SeaBasisBladeIndexLow(14))/sf,
			self.0.basis_blade_coefficient(&SeaBasisBladeIndexLow(13))/sf,
			-self.0.basis_blade_coefficient(&SeaBasisBladeIndexLow(11))/sf,
		];

		for o in &mut out {
			if o.abs() < T::from(1e-15).unwrap() {
				*o = T::zero()
			}
		}
		out
	}

	fn simplify(&mut self) {
		self.0.remove_approx_zeros(T::from(1e-15).unwrap());
		self.0.round_to_decimals(14);
	}

	pub fn raw(&self) -> (u16, [T;16]) {
		let mut out_mat = [T::zero(); 16];
		for i in 0..16 {
			out_mat[i] = self.0.basis_blade_coefficient(&SeaBasisBladeIndexLow(i as u8))
		}
		let out_bbb = self.0.multivector.basis_blade_bits;
		(out_bbb, out_mat)
	}

	pub fn scalar_norm(&self) -> T{
		self.0.scalar_norm()
	}

	pub fn scalar_normalized(&self) -> Self{
		Self(self.0.scalar_normalized())
	}

	pub fn scalar_normalize(&mut self) {
		self.0.scalar_normalize();
	}

	pub fn reversion_conjugated(&self) -> Self{
		Self(self.0.reversion_conjugated())
	}

	pub fn exp(&self) -> Self {
		Self(self.0.exp())
	}

	pub fn scalar_product(&self, rhs: &Self) -> T {
		self.0.scalar_product(&rhs.0)
	}
}

impl<T:SeaCoefficient> From<SeaMotorPGA3D<T>> for MultivectorPGA3D<T> {
	fn from(value: SeaMotorPGA3D<T>) -> Self {
		value.0
	}
}
impl<T:SeaCoefficient> From<&SeaMotorPGA3D<T>> for MultivectorPGA3D<T> {
	fn from(value: &SeaMotorPGA3D<T>) -> Self {
		value.0.clone()
	}
}

impl<T:SeaCoefficient> Not for MultivectorPGA3D<T> {
	type Output = MultivectorPGA3D<T>;

	fn not(self) -> Self::Output {
		Self(self.0.dual())
	}
}
impl<T:SeaCoefficient> Not for &MultivectorPGA3D<T> {
	type Output = MultivectorPGA3D<T>;

	fn not(self) -> Self::Output {
		MultivectorPGA3D(self.0.dual())
	}
}

impl<T:SeaCoefficient> Neg for &MultivectorPGA3D<T> {
	type Output = MultivectorPGA3D<T>;

	fn neg(self) -> Self::Output {
		MultivectorPGA3D(-&(self.0))
	}
}
impl<T:SeaCoefficient> Neg for MultivectorPGA3D<T> {
	type Output = MultivectorPGA3D<T>;

	fn neg(self) -> Self::Output {
		-&self
	}
}

impl<T:SeaCoefficient> Shl<Self> for MultivectorPGA3D<T> {
	type Output = MultivectorPGA3D<T>;
	fn shl(self, rhs: Self) -> Self::Output {
		MultivectorPGA3D(self.0 << rhs.0)
	}
}
impl<T:SeaCoefficient> Shl<&Self> for MultivectorPGA3D<T> {
	type Output = MultivectorPGA3D<T>;
	fn shl(self, rhs: &Self) -> Self::Output {
		MultivectorPGA3D(self.0 << &rhs.0)
	}
}
impl<T:SeaCoefficient> Shl<MultivectorPGA3D<T>> for &MultivectorPGA3D<T> {
	type Output = MultivectorPGA3D<T>;
	fn shl(self, rhs: MultivectorPGA3D<T>) -> Self::Output {
		MultivectorPGA3D(&self.0 << rhs.0)
	}
}
impl<T:SeaCoefficient> Shl<Self> for &MultivectorPGA3D<T> {
	type Output = MultivectorPGA3D<T>;
	fn shl(self, rhs: Self) -> Self::Output {
		MultivectorPGA3D(&self.0 << &rhs.0)
	}
}

impl<T:SeaCoefficient> Shr<Self> for MultivectorPGA3D<T> {
	type Output = MultivectorPGA3D<T>;
	fn shr(self, rhs: Self) -> Self::Output {
		MultivectorPGA3D(self.0 >> rhs.0)
	}
}
impl<T:SeaCoefficient> Shr<&Self> for MultivectorPGA3D<T> {
	type Output = MultivectorPGA3D<T>;
	fn shr(self, rhs: &Self) -> Self::Output {
		MultivectorPGA3D(self.0 >> &rhs.0)
	}
}
impl<T:SeaCoefficient> Shr<MultivectorPGA3D<T>> for &MultivectorPGA3D<T> {
	type Output = MultivectorPGA3D<T>;
	fn shr(self, rhs: MultivectorPGA3D<T>) -> Self::Output {
		MultivectorPGA3D(&self.0 >> rhs.0)
	}
}
impl<T:SeaCoefficient> Shr<Self> for &MultivectorPGA3D<T> {
	type Output = MultivectorPGA3D<T>;
	fn shr(self, rhs: Self) -> Self::Output {
		MultivectorPGA3D(&self.0 >> &rhs.0)
	}
}

impl<T:SeaCoefficient> BitOr<Self> for MultivectorPGA3D<T> {
	type Output = MultivectorPGA3D<T>;
	fn bitor(self, rhs: Self) -> Self::Output {
		MultivectorPGA3D(self.0 | rhs.0)
	}
}
impl<T:SeaCoefficient> BitOr<&Self> for MultivectorPGA3D<T> {
	type Output = MultivectorPGA3D<T>;
	fn bitor(self, rhs: &Self) -> Self::Output {
		MultivectorPGA3D(self.0 | &rhs.0)
	}
}
impl<T:SeaCoefficient> BitOr<MultivectorPGA3D<T>> for &MultivectorPGA3D<T> {
	type Output = MultivectorPGA3D<T>;
	fn bitor(self, rhs: MultivectorPGA3D<T>) -> Self::Output {
		MultivectorPGA3D(&self.0 | rhs.0)
	}
}
impl<T:SeaCoefficient> BitOr<Self> for &MultivectorPGA3D<T> {
	type Output = MultivectorPGA3D<T>;
	fn bitor(self, rhs: Self) -> Self::Output {
		MultivectorPGA3D(&self.0 | &rhs.0)
	}
}

impl<T:SeaCoefficient> BitXor<Self> for MultivectorPGA3D<T> {
	type Output = MultivectorPGA3D<T>;
	fn bitxor(self, rhs: Self) -> Self::Output {
		MultivectorPGA3D(self.0 ^ rhs.0)
	}
}
impl<T:SeaCoefficient> BitXor<&Self> for MultivectorPGA3D<T> {
	type Output = MultivectorPGA3D<T>;
	fn bitxor(self, rhs: &Self) -> Self::Output {
		MultivectorPGA3D(self.0 ^ &rhs.0)
	}
}
impl<T:SeaCoefficient> BitXor<MultivectorPGA3D<T>> for &MultivectorPGA3D<T> {
	type Output = MultivectorPGA3D<T>;
	fn bitxor(self, rhs: MultivectorPGA3D<T>) -> Self::Output {
		MultivectorPGA3D(&self.0 ^ rhs.0)
	}
}
impl<T:SeaCoefficient> BitXor<Self> for &MultivectorPGA3D<T> {
	type Output = MultivectorPGA3D<T>;
	fn bitxor(self, rhs: Self) -> Self::Output {
		MultivectorPGA3D(&self.0 ^ &rhs.0)
	}
}

impl<T:SeaCoefficient> BitAnd<Self> for MultivectorPGA3D<T> {
	type Output = MultivectorPGA3D<T>;
	
	fn bitand(self, rhs: Self) -> Self::Output {
		MultivectorPGA3D(self.0 & rhs.0)
	}
}
impl<T:SeaCoefficient> BitAnd<&Self> for MultivectorPGA3D<T> {
	type Output = MultivectorPGA3D<T>;
	
	fn bitand(self, rhs: &Self) -> Self::Output {
		MultivectorPGA3D(self.0 & &rhs.0)
	}
}
impl<T:SeaCoefficient> BitAnd<MultivectorPGA3D<T>> for &MultivectorPGA3D<T> {
	type Output = MultivectorPGA3D<T>;
	
	fn bitand(self, rhs: MultivectorPGA3D<T>) -> Self::Output {
		MultivectorPGA3D(&self.0 & rhs.0)
	}
}
impl<T:SeaCoefficient> BitAnd<Self> for &MultivectorPGA3D<T> {
	type Output = MultivectorPGA3D<T>;
	
	fn bitand(self, rhs: Self) -> Self::Output {
		MultivectorPGA3D(&self.0 & &rhs.0)
	}
}

impl<T:SeaCoefficient> Mul<Self> for MultivectorPGA3D<T> {
	type Output = MultivectorPGA3D<T>;
	fn mul(self, rhs: Self) -> Self::Output {
		MultivectorPGA3D(self.0 * rhs.0)
	}
}
impl<T:SeaCoefficient> Mul<&Self> for MultivectorPGA3D<T> {
	type Output = MultivectorPGA3D<T>;
	fn mul(self, rhs: &Self) -> Self::Output {
		MultivectorPGA3D(self.0 * &rhs.0)
	}
}
impl<T:SeaCoefficient> Mul<MultivectorPGA3D<T>> for &MultivectorPGA3D<T> {
	type Output = MultivectorPGA3D<T>;
	fn mul(self, rhs: MultivectorPGA3D<T>) -> Self::Output {
		MultivectorPGA3D(&self.0 * rhs.0)
	}
}
impl<T:SeaCoefficient> Mul<Self> for &MultivectorPGA3D<T> {
	type Output = MultivectorPGA3D<T>;
	fn mul(self, rhs: Self) -> Self::Output {
		MultivectorPGA3D(&self.0 * &rhs.0)
	}
}

impl<T:SeaCoefficient> AddAssign<T> for MultivectorPGA3D<T> {
    fn add_assign(&mut self, rhs: T){
		self.0 += rhs
	}
}
impl<T:SeaCoefficient> Add<T> for MultivectorPGA3D<T> {
    type Output = MultivectorPGA3D<T>;

    fn add(self, rhs: T) -> Self::Output {
		&self + rhs
	}
}
impl<T:SeaCoefficient> Add<T> for &MultivectorPGA3D<T> {
    type Output = MultivectorPGA3D<T>;

    fn add(self, rhs: T) -> Self::Output {
		MultivectorPGA3D(&self.0 + rhs)
	}
}

impl<T:SeaCoefficient> SubAssign<T> for MultivectorPGA3D<T> {
    fn sub_assign(&mut self, rhs: T){
		self.0 -= rhs
	}
}
impl<T:SeaCoefficient> Sub<T> for MultivectorPGA3D<T> {
    type Output = MultivectorPGA3D<T>;

    fn sub(self, rhs: T) -> Self::Output {
		&self - rhs
	}
}
impl<T:SeaCoefficient> Sub<T> for &MultivectorPGA3D<T> {
    type Output = MultivectorPGA3D<T>;

    fn sub(self, rhs: T) -> Self::Output {
		MultivectorPGA3D(&self.0 - rhs)
	}
}

impl<T:SeaCoefficient> MulAssign<T> for MultivectorPGA3D<T> {
    fn mul_assign(&mut self, rhs: T){
		self.0 *= rhs
	}
}
impl<T:SeaCoefficient> Mul<T> for MultivectorPGA3D<T> {
    type Output = MultivectorPGA3D<T>;

    fn mul(self, rhs: T) -> Self::Output {
		&self * rhs
	}
}
impl<T:SeaCoefficient> Mul<T> for &MultivectorPGA3D<T> {
    type Output = MultivectorPGA3D<T>;

    fn mul(self, rhs: T) -> Self::Output {
		MultivectorPGA3D(&self.0 * rhs)
	}
}

impl<T:SeaCoefficient> DivAssign<T> for MultivectorPGA3D<T> {
    fn div_assign(&mut self, rhs: T){
		self.0 /= rhs
	}
}
impl<T:SeaCoefficient> Div<T> for MultivectorPGA3D<T> {
    type Output = MultivectorPGA3D<T>;

    fn div(self, rhs: T) -> Self::Output {
		&self / rhs
	}
}
impl<T:SeaCoefficient> Div<T> for &MultivectorPGA3D<T> {
    type Output = MultivectorPGA3D<T>;

    fn div(self, rhs: T) -> Self::Output {
		MultivectorPGA3D(&self.0 / rhs)
	}
}


#[cfg(test)]
mod tests {
	use std::f64::consts::TAU;

	use super::*;

    #[test]
    fn create_plane_from_points() {
		let p1 = MultivectorPGA3D::new_point(1.0,0.0,0.0);
		let p2 = MultivectorPGA3D::new_point(0.0,1.0,0.0);
		let p3 = MultivectorPGA3D::new_point(0.0,0.0,0.0);
		let c = p1&p2&p3;
		let out = MultivectorPGA3D::new_plane(0.0,0.0,1.0, 0.0);

		assert_eq!(
			c,
			out
		);

		let p1 = MultivectorPGA3D::new_point(1.0,2.0,3.0);
		let p2 = MultivectorPGA3D::new_point(2.0,1.0,4.0);
		let p3 = MultivectorPGA3D::new_point(2.0,5.0,7.0);
		let c = p1&p2&p3;
		let out = MultivectorPGA3D::new_plane(-7.0,-3.0,4.0, 1.0);

		assert_eq!(
			c,
			out
		);

		let p1 = MultivectorPGA3D::new_point(1.0,2.0,3.0);
		let p2 = MultivectorPGA3D::new_point(2.0,1.0,4.0);
		let p3 = MultivectorPGA3D::new_point(2.0,5.0,7.0);
		let c = p1&p2&p3;
		let out = MultivectorPGA3D::new_plane(-7.0,-3.0,4.0, 1.0);

		assert_eq!(
			c,
			out
		);
	}

	#[test]
	fn point_to_vec() {
		let p1 = MultivectorPGA3D::new_point(1.0,2.0,3.0);
		let v = p1.point_to_vec3();
		assert_eq!(
			v,
			[1.0, 2.0, 3.0]
		);
	}

	#[test]
	fn project_point_onto_plane() {
		let plane = MultivectorPGA3D::new_plane(1.0,0.0,0.0, 0.0);
		let p = MultivectorPGA3D::new_point(1.0,1.0,1.0);
		let proj = (&plane|p)*plane;
		
		let out = [0.0, 1.0, 1.0];
		assert_eq!(
			proj.point_to_vec3(),
			out
		);

		let plane = MultivectorPGA3D::new_plane(0.0,1.0,0.0, 0.0);
		let p = MultivectorPGA3D::new_point(1.0,1.0,1.0);
		let proj = (&plane|p)*plane;
		
		let out = [1.0, 0.0, 1.0];
		assert_eq!(
			proj.point_to_vec3(),
			out
		);

		let plane = MultivectorPGA3D::new_plane(0.0,0.0,1.0, 0.0);
		let p = MultivectorPGA3D::new_point(1.0,1.0,1.0);
		let proj = (&plane|p)*plane;
		
		let out = [1.0, 1.0, 0.0];
		assert_eq!(
			proj.point_to_vec3(),
			out
		);

		let plane = MultivectorPGA3D::new_plane(-7.0,-3.0,4.0, 1.0);
		let p4 = MultivectorPGA3D::new_point(0.0,0.0,0.0);
		let proj = (&plane|p4)*plane;

		let out = [0.0945945945945946, 0.04054054054054054, -0.05405405405405406];
		assert_eq!(
			proj.point_to_vec3(),
			out
		);

		let plane = MultivectorPGA3D::new_plane(1.0,0.0,0.0, -1.0);
		let p4 = MultivectorPGA3D::new_point(100.0,0.0,0.0);
		let proj = (&plane|p4)*plane;

		let out = [1.0, 0.0, 0.0];
		assert_eq!(
			proj.point_to_vec3(),
			out
		);
	}

	#[test]
	fn rotor() {
		let l = MultivectorPGA3D::new_line(5.0,3.0,2.0);
		let r = MultivectorPGA3D::new_motor(&l, TAU/1.6);
		let r_inv = r.inversed();
		let mut approx_one = r.into() * r_inv.into();
		approx_one.simplify();
		assert_eq!(
			approx_one,
			MultivectorPGA3D::new_lifted_scalar(1.0),
		);


		let l = MultivectorPGA3D::new_line(0.0,0.0,1.0);
		let r = MultivectorPGA3D::new_motor(&l, TAU/2.0);
		let r_inv = r.inversed();
		let p = MultivectorPGA3D::new_point(1.0, 0.0, 0.0);
		let p_rot = r.into() * p * r_inv.into();

		let out = [-1.0, 0.0, 0.0];

		assert_eq!(
			p_rot.point_to_vec3(),
			out
		)



	}

	#[test]
	fn translator() {
		let l = MultivectorPGA3D::new_ideal_line(2.0,3.0,5.0);
		let r = MultivectorPGA3D::new_motor(&l, TAU/3.2);
		let r_inv = r.inversed();
		let mut approx_one = r.into() * r_inv.into();
		approx_one.simplify();
		assert_eq!(
			approx_one,
			MultivectorPGA3D::new_lifted_scalar(1.0),
		);


		let l = MultivectorPGA3D::new_ideal_line(1.0,0.0,0.0);
		let r = MultivectorPGA3D::new_motor(&l, 4.0);
		let r_inv = r.inversed();
		let p = MultivectorPGA3D::new_point(0.0, 0.0, 0.0);
		let p_trans = r.into() * p * r_inv.into();

		let out = [4.0, 0.0, 0.0];

		assert_eq!(
			p_trans.point_to_vec3(),
			out
		)



	}

	mod camera {
		use super::*;

		#[test]
		fn simple_projection() {
			let image_plane = MultivectorPGA3D::new_plane(0.0,0.0,1.0, 0.0);
			let focal_point = MultivectorPGA3D::new_point(0.0, 0.0, -1.0);
	
			let translation_dir = MultivectorPGA3D::new_ideal_line(1.0,1.0,1.0);
			let camera_pos = MultivectorPGA3D::new_motor(&translation_dir, 0.0);
	
			let rotation_axis = MultivectorPGA3D::new_line(1.0,0.0,0.0);
			let camera_rot = MultivectorPGA3D::new_motor(&rotation_axis, 0.0);
	
			let camera_trans = SeaMotorPGA3D(camera_rot.into() * camera_pos.into());
			let camera_trans_inv = camera_trans.inversed();
			
			let scene_point = MultivectorPGA3D::new_point(0.0,3.0,2.0);
			let scene_point_trans = camera_trans.into() * scene_point * camera_trans_inv.into();
	
			let fp_to_scene = focal_point&scene_point_trans;
	
			let screen_point = fp_to_scene^image_plane;
	
			assert_eq!(
				screen_point.point_to_vec3(),
				[0.0, 1.0, 0.0]
			);
		}

		#[test]
		fn translation() {
			let image_plane = MultivectorPGA3D::new_plane(0.0,0.0,1.0, 0.0);
			let focal_point = MultivectorPGA3D::new_point(0.0, 0.0, -1.0);
	
			let translation_dir = MultivectorPGA3D::new_ideal_line(1.0,1.0,0.0);
			let camera_pos = MultivectorPGA3D::new_motor(&translation_dir, -3.0);
	
			let rotation_axis = MultivectorPGA3D::new_line(1.0,0.0,0.0);
			let camera_rot = MultivectorPGA3D::new_motor(&rotation_axis, 0.0);
	
			let camera_trans = SeaMotorPGA3D(camera_rot.into() * camera_pos.into());
			let camera_trans_inv = camera_trans.inversed();
			
			let scene_point = MultivectorPGA3D::new_point(0.0,3.0,2.0);
			let scene_point_trans = camera_trans.into() * scene_point * camera_trans_inv.into();
	
			let fp_to_scene = focal_point&scene_point_trans;
	
			let screen_point = fp_to_scene^image_plane;
	
			assert_eq!(
				screen_point.point_to_vec3(),
				[-1.0, 0.0, 0.0]
			);
		}

		
		#[test]
		fn rotation() {
			let image_plane = MultivectorPGA3D::new_plane(0.0,0.0,1.0, 0.0);
			let focal_point = MultivectorPGA3D::new_point(0.0, 0.0, -1.0);
	
			let translation_dir = MultivectorPGA3D::new_ideal_line(0.0,0.0,0.0);
			let camera_pos = MultivectorPGA3D::new_motor(&translation_dir, -3.0);
	
			let rotation_axis = MultivectorPGA3D::new_line(0.0,0.0,1.0);
			let camera_rot = MultivectorPGA3D::new_motor(&rotation_axis, -TAU/4.0);
	
			let camera_trans = SeaMotorPGA3D(camera_rot.into() * camera_pos.into());
			let camera_trans_inv = camera_trans.inversed();
			
			let scene_point = MultivectorPGA3D::new_point(0.0,3.0,2.0);
			let scene_point_trans = camera_trans.into() * scene_point * camera_trans_inv.into();
	
			let fp_to_scene = &focal_point&scene_point_trans;
	
			let mut screen_point = fp_to_scene^&image_plane;
			
			screen_point.simplify();

			assert_eq!(
				screen_point.point_to_vec3(),
				[-1.0, 0.0, 0.0]
			);

			let image_plane = MultivectorPGA3D::new_plane(0.0,0.0,1.0, 0.0);
			let focal_point = MultivectorPGA3D::new_point(0.0, 0.0, -1.0);
	
			let translation_dir = MultivectorPGA3D::new_ideal_line(0.0,0.0,0.0);
			let camera_pos = MultivectorPGA3D::new_motor(&translation_dir, 0.0);
	
			let rotation_axis = MultivectorPGA3D::new_line(0.0,0.0,1.0);
			let camera_rot = MultivectorPGA3D::new_motor(&rotation_axis, -TAU/2.0);
	
			let camera_trans = SeaMotorPGA3D(camera_rot.into() * camera_pos.into());
			let camera_trans_inv = camera_trans.inversed();
			
			let scene_point = MultivectorPGA3D::new_point(0.0,3.0,2.0);
			let scene_point_trans = camera_trans.into() * scene_point * camera_trans_inv.into();
	
			let fp_to_scene = &focal_point&scene_point_trans;
	
			let screen_point = fp_to_scene^&image_plane;
	
			assert_eq!(
				screen_point.point_to_vec3(),
				[0.0, -1.0, 0.0]
			);
		}

		
		#[test]
		fn rotate_and_trans() {
			let image_plane = MultivectorPGA3D::new_plane(0.0,0.0,1.0, 0.0);
			let focal_point = MultivectorPGA3D::new_point(0.0, 0.0, -1.0);
	
			let translation_dir = MultivectorPGA3D::new_ideal_line(1.0,1.0,0.0);
			let camera_pos = MultivectorPGA3D::new_motor(&translation_dir, -3.0);
	
			let rotation_axis = MultivectorPGA3D::new_line(0.0,0.0,1.0);
			let camera_rot = MultivectorPGA3D::new_motor(&rotation_axis, -TAU/4.0);
			
			let camera_trans = SeaMotorPGA3D(camera_rot.into() * camera_pos.into());
			let camera_trans_inv = camera_trans.inversed();
			
			let scene_point = MultivectorPGA3D::new_point(0.0,3.0,2.0);
			let scene_point_trans = camera_trans.into() * scene_point * camera_trans_inv.into();
	
			let fp_to_scene = &focal_point&scene_point_trans;
	
			let mut screen_point = fp_to_scene^&image_plane;
			
			screen_point.simplify();

			assert_eq!(
				screen_point.point_to_vec3(),
				[0.0, -1.0, 0.0]
			);
		}
	
		#[test]
		fn precomputed() {
			let image_plane = MultivectorPGA3D::new_plane(0.0,0.0,1.0, 0.0);
			let focal_point = MultivectorPGA3D::new_point(0.0, 0.0, -1.0);
	
			let translation_dir = MultivectorPGA3D::new_ideal_line(1.0,1.0,0.0);
			let camera_pos = MultivectorPGA3D::new_motor(&translation_dir, -3.0);
	
			let rotation_axis = MultivectorPGA3D::new_line(0.0,0.0,1.0);
			let camera_rot = MultivectorPGA3D::new_motor(&rotation_axis, -TAU/4.0);
			
			let camera_trans = SeaMotorPGA3D(camera_rot.into() * camera_pos.into());
			let camera_trans_inv = camera_trans.inversed();
			
			let scene_point = MultivectorPGA3D::new_point(0.0,3.0,2.0);
			let scene_point_trans = camera_trans.into() * scene_point * camera_trans_inv.into();

			let mut screen_point = focal_point&scene_point_trans^image_plane;

			screen_point.simplify();

			assert_eq!(
				screen_point.point_to_vec3(),
				[0.0, -1.0, 0.0]
			);
		}

		#[test]
		fn sanity_check() {
			let image_plane = MultivectorPGA3D::new_plane(0.0,0.0,1.0, 0.0);
			let mut focal_point = MultivectorPGA3D::new_point(0.0, 0.0, -1.0);
			focal_point.simplify();

			let translation_dir = MultivectorPGA3D::new_ideal_line(0.0,0.0,0.0);
			let mut camera_pos = MultivectorPGA3D::new_motor(&translation_dir, 0.0);
	
			let rotation_axis = MultivectorPGA3D::new_line(0.0,0.0,0.0);
			let mut camera_rot = MultivectorPGA3D::new_motor(&rotation_axis, 0.0);
			
			let mut camera_trans = SeaMotorPGA3D(camera_rot.into() * camera_pos.into());
			let camera_trans_inv = camera_trans.inversed();


			
			let scene_point = MultivectorPGA3D::new_point(-0.08682,0.4924,0.00);
			let scene_point_trans = camera_trans.into() * scene_point * camera_trans_inv.into();

			let mut screen_point2 = &focal_point & &scene_point_trans;

			let mut screen_point = focal_point&scene_point_trans^image_plane;

			screen_point.simplify();

			assert_eq!(
				screen_point.point_to_vec3(),
				[-0.08682,0.4924, 0.0]
			);
		}
	}
}