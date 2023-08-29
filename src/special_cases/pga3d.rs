use std::ops::{Deref, DerefMut, Shl, Shr, BitOr, BitXor, Mul, Not, BitAnd, Index, AddAssign, Add, Neg, Sub, Div, SubAssign, MulAssign, DivAssign};

use crate::coefficient::SeaCoefficient;
use crate::multivector_traits::Multivector;
use crate::grade::SeaGrade;
use crate::low_dimensional_multivector::Multivector4D;

use crate::low_dimensional_multivector::SeaBasisBladeBits;
use crate::low_dimensional_multivector::SeaBasisBladeIndexLow;

#[derive(Debug, PartialEq)]
struct MultivectorPGA3D<T:SeaCoefficient>(Multivector4D<T>);

const POINT_BASIS_BLADES: SeaBasisBladeBits<u16> = SeaBasisBladeBits::<u16>::GRADES[3];
const PLANE_BASIS_BLADES: SeaBasisBladeBits<u16> = SeaBasisBladeBits::<u16>::GRADES[1];
const LINE_BASIS_BLADES: SeaBasisBladeBits<u16> = SeaBasisBladeBits::<u16>(0x0068);
const IDEAL_LINE_BASIS_BLADES: SeaBasisBladeBits<u16> = SeaBasisBladeBits::<u16>(0x1600);

impl<T:SeaCoefficient> MultivectorPGA3D<T>{
	fn new() -> Self{
		Self(Multivector4D::new(3, 0, 1))
	}

	fn new_point(x:T, y:T, z:T) -> Self {
		let mut mv = Self::new();
		mv.0.multivector.basis_blade_bits = POINT_BASIS_BLADES;
		mv.0.multivector.coefficients = vec![T::one(),-z,y,-x];
		mv.0.multivector.simplify();
		mv
	}

	fn new_plane(x:T, y:T, z:T, w:T) -> Self {
		let mut mv = Self::new();
		mv.0.multivector.basis_blade_bits = PLANE_BASIS_BLADES;
		mv.0.multivector.coefficients = vec![x,y,z,w];
		mv.0.multivector.simplify();
		mv
	}

	fn new_line(a:T, b:T, c:T) -> Self {
		let mut mv = Self::new();
		mv.0.multivector.basis_blade_bits = LINE_BASIS_BLADES;
		mv.0.multivector.coefficients = vec![a,-b,c];
		mv.0.multivector.simplify();
		mv
	}

	fn new_ideal_line(d:T, e:T, f:T) -> Self {
		let mut mv = Self::new();
		mv.0.multivector.basis_blade_bits = IDEAL_LINE_BASIS_BLADES;
		mv.0.multivector.coefficients = vec![d,e,f];
		mv.0.multivector.simplify();
		mv
	}

	fn point_to_vec3(&self) -> [T;3] {
		let sf = self.0.basis_blade_coefficient(&SeaBasisBladeIndexLow(7));
		[
			-self.0.basis_blade_coefficient(&SeaBasisBladeIndexLow(14))/sf,
			self.0.basis_blade_coefficient(&SeaBasisBladeIndexLow(13))/sf,
			-self.0.basis_blade_coefficient(&SeaBasisBladeIndexLow(11))/sf,
		]
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


#[cfg(test)]
mod tests {
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
	fn point_to_vec() {
		let p1 = MultivectorPGA3D::new_point(1.0,2.0,3.0);
		let v = p1.point_to_vec3();
		assert_eq!(
			v,
			[
				1.0,
				2.0,
				3.0
			]
		);
	}
}