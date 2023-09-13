use super::*;

mod scalar_operations;
mod linear_operations;
mod non_geometric_linear_operations;

#[derive(Debug, PartialEq)]
pub struct SeaMotorPGA3D<T:SeaCoefficient>(pub SeaMultivectorPGA3D<T>);
impl<T:SeaCoefficient> SeaMotorPGA3D<T> {
	pub fn into(self) -> SeaMultivectorPGA3D<T> {
		self.0
	}
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

impl<T:SeaCoefficient> Motor<T> for SeaMotorPGA3D<T> {}

impl<T:SeaCoefficient> Versor<T> for SeaMotorPGA3D<T> {
	fn inverse(&self) -> Self {
		Self((**self).reversion_conjugated())
	}
}

impl<T:SeaCoefficient> Multivector<T> for SeaMotorPGA3D<T> {
	type MultivectorIndexType = SeaBasisBladeIndexLow;

	fn basis_blade_coefficient(&self, index: &Self::MultivectorIndexType) -> T {
		(**self).basis_blade_coefficient(index)
	}


	fn reversion_conjugated(&self) -> Self {
		Self((**self).reversion_conjugated())
	}

	fn involution_conjugated(&self) -> Self {
		Self((**self).involution_conjugated())
	}

	fn clifford_conjugated(&self) -> Self {
		Self((**self).clifford_conjugated())
	}


	fn poincare_dual(&self) -> Self {
		Self((**self).poincare_dual())
	}

	fn right_complimented(&self) -> Self {
		Self((**self).right_complimented())
	}

	fn left_complimented(&self) -> Self {
		Self((**self).left_complimented())
	}

	fn dual(&self) -> Self {
		Self((**self).dual())
	}

	fn undual(&self) -> Self {
		Self((**self).undual())
	}

	fn scalar_norm(&self) -> T {
		(**self).scalar_norm()
	}

	fn scalar_normalized(&self) -> Self {
		Self((**self).scalar_normalized())
	}


	fn exp(&self) -> Self {
		Self((**self).exp())
	}
}
impl<T:SeaCoefficient> MultivectorMut<T> for SeaMotorPGA3D<T> {
	
	fn reversion_conjugation(&mut self) {
		(**self).reversion_conjugation()
	}

	fn involution_conjugation(&mut self) {
		(**self).involution_conjugation()
	}

	fn clifford_conjugation(&mut self) {
		(**self).clifford_conjugation()
	}
	

	fn poincare_dualization(&mut self) {
		(**self).poincare_dualization()
	}

	fn right_compliment(&mut self) {
		(**self).right_compliment()
	}

	fn left_compliment(&mut self) {
		(**self).left_compliment()
	}

	fn dualization(&mut self) {
		(**self).dualization()
	}

	fn undualization(&mut self) {
		(**self).undualization()
	}
	
	fn scalar_normalize(&mut self) {
		(**self).scalar_normalize()
	}

}

impl<T:SeaCoefficient> Neg for SeaMotorPGA3D<T> {
	type Output = Self;

	fn neg(self) -> Self::Output {
		Self(-&*self)
	}
}
impl<T:SeaCoefficient> Neg for &SeaMotorPGA3D<T> {
	type Output = SeaMotorPGA3D<T>;

	fn neg(self) -> Self::Output {
		SeaMotorPGA3D::<T>(-&**self)
	}
}

impl<T:SeaCoefficient> Not for SeaMotorPGA3D<T> {
	type Output = SeaMotorPGA3D<T>;

	fn not(self) -> Self::Output {
		self.dual()
	}
}
impl<T:SeaCoefficient> Not for &SeaMotorPGA3D<T> {
	type Output = SeaMotorPGA3D<T>;

	fn not(self) -> Self::Output {
		self.dual()
	}
}
