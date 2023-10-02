use super::*;

impl<'a, T:SeaCoefficient, I:SeaBasisBladeSet> MultivectorGradeIter for SeaLowDimensionalMultivectorCustomMetric<T, I> {
	fn grade_iter(&self, grade: &SeaGrade) -> <&Self as MultivectorGradeIntoIterator>::MultivectorGradeIntoIter {
		self.grade_into_iter(grade)
	}
}

impl<'a, T:SeaCoefficient, I:SeaBasisBladeSet> MultivectorGradeIterMut for SeaLowDimensionalMultivectorCustomMetric<T, I> {
	fn grade_iter_mut(&mut self, grade: &SeaGrade) -> <&mut Self as MultivectorGradeIntoIterator>::MultivectorGradeIntoIter {
		self.grade_into_iter(grade)
	}
}

impl<'a, T:SeaCoefficient, I:SeaBasisBladeSet> MultivectorGradeIntoIterator for &'a SeaLowDimensionalMultivectorCustomMetric<T, I> {
	type MultivectorGradeIntoIter = SparseMultivectorIter<'a, T, I>;

	fn grade_into_iter(self, grade: &SeaGrade) -> Self::MultivectorGradeIntoIter {
		self.multivector.grade_iter(grade)
	}
}
impl<'a, T:SeaCoefficient, I:SeaBasisBladeSet> MultivectorGradeIntoIterator for &'a mut SeaLowDimensionalMultivectorCustomMetric<T, I> {
	type MultivectorGradeIntoIter = SparseMultivectorIterMut<'a, T, I>;

	fn grade_into_iter(self, grade: &SeaGrade) -> Self::MultivectorGradeIntoIter {
		self.multivector.grade_iter_mut(grade)
	}
}