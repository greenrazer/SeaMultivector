use super::*;

//log
impl<T:SeaCoefficient> SeaMotorPGA3D<T> {
	pub fn log(&self) -> SeaMotorPGA3D<T> {
		let b = &self.grade(&SeaGrade(2));
		let s = (-b | b).scalar().sqrt();

		// pure translation
		if s.is_zero(){
			return Self((&**self)/self.scalar() - T::one());
		}

		// turn
		if self.scalar().is_zero() {
			let tau_div_4 = T::from(TAU/4.0).unwrap();
			return  Self(-self.pseudoscalar()/s + tau_div_4);
		}

		let two = T::from(2.0).unwrap();
		let p = &( (-b ^ b) / (two*s) );
		
		let lhs = (p / self.scalar()) + (s / self.scalar()).atan();
		let rhs = b*((-p+s)/(s*s));

		Self(lhs * rhs)
	}
}

#[cfg(test)]
mod tests {
	use std::f64::consts::TAU;
	
	pub use super::super::tests::assert_approx_eq_vec4;
	pub use super::super::tests::assert_approx_eq_mat4;
	use super::multivector::tests::assert_approx_eq_mv;

	use super::*;

	#[test]
	fn log() {
		let l = SeaMultivectorPGA3D::new_ideal_line(1.0,TAU/4.0,10.0);
		let r = SeaMultivectorPGA3D::new_motor(&l, 2.0);

		let mut lout = r.log();

		lout.simplify();

		assert_approx_eq_mv(l.into(), lout.into());

		let l = SeaMultivectorPGA3D::new_line(0.0,TAU/4.0,0.0);
		let r = SeaMultivectorPGA3D::new_motor(&l, 2.0);

		let mut lout = r.log();

		lout.simplify();

		assert_approx_eq_mv(l.into(), lout.into());

		let mut l = SeaMultivectorPGA3D::new_line(TAU/8.0,TAU/4.0,TAU/3.2);
		l.scalar_normalize();
		let r = SeaMultivectorPGA3D::new_motor(&l, 2.0);

		let mut lout = r.log();

		lout.simplify();

		assert_approx_eq_mv(l.into(), lout.into());
	}
}