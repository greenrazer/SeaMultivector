use std::{f64::consts::TAU, arch::x86_64::_MM_FLUSH_ZERO_MASK};

use super::*;

mod scalar_operations;
mod linear_operations;
mod non_linear_operations;
mod non_geometric_linear_operations;

mod quaternion;
mod matrix;

mod iterators;
mod from;


#[derive(Debug, Clone, PartialEq)]
pub struct SeaMotorPGA3D<T:SeaCoefficient>(pub SeaMultivectorPGA3D<T>);
impl<T:SeaCoefficient> SeaMotorPGA3D<T> {
	pub fn identity() -> Self {
		Self::from(T::one())
	}

	pub fn new_translation_motor(x:T, y:T, z:T) -> SeaMotorPGA3D<T> {
		let line = SeaMultivectorPGA3D::new_ideal_line(x, y, z);
		SeaMultivectorPGA3D::new_motor(&line, T::one())
	}

	pub fn new_rotation_motor(rx:T, ry:T, rz:T) -> SeaMotorPGA3D<T> {
		let line = SeaMultivectorPGA3D::new_line(rx, ry, rz);
		SeaMultivectorPGA3D::new_motor(&line, T::one())
	}

	pub fn new_combined_motor(x:T, y:T, z:T, rx:T, ry:T, rz:T, order: MotorCreationOrder) -> SeaMotorPGA3D<T> {
		let motor_trans = Self::new_translation_motor(x, y, z);
		let motor_rot = Self::new_rotation_motor(rx, ry, rz);
		match order {
			MotorCreationOrder::TranslationFirst => motor_rot * motor_trans,
			MotorCreationOrder::RotationFirst => motor_trans * motor_rot
		}
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
}
impl<T:SeaCoefficient> MultivectorMut<T> for SeaMotorPGA3D<T> {}

impl<T: SeaCoefficient> MultivectorConstMetric<T> for SeaMotorPGA3D<T> {
	const METRIC: SeaMetric = SeaMultivectorPGA3D::<T>::METRIC;
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

impl<T:SeaCoefficient> Default for SeaMotorPGA3D<T> {
	fn default() -> Self {
		Self(SeaMultivectorPGA3D::<T>::default())
	}
}
