use super::*;

// Regressive Product
impl<T:SeaCoefficient, I:SeaBasisBladeBitsInner> SeaMultivectorLowSparseCoefficents<T, I> {
	pub fn regressive_product(&self, rhs: &Self, metric: &SeaMetric) -> Self {
		let (a_dual, b_dual) = (self.dual(metric), rhs.dual(metric));

		let mut out_dual = &a_dual ^ &b_dual;
		
		out_dual.undualization(metric);
		
		out_dual
	}
}
