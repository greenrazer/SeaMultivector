use super::*;

impl<'a, T:SeaCoefficient> MultivectorGradeIter for SeaFlectorPGA3D<T> {
	fn grade_iter(&self, grade: &SeaGrade) -> <&Self as MultivectorGradeIntoIterator>::MultivectorGradeIntoIter {
		self.grade_into_iter(grade)
	}
}

impl<'a, T:SeaCoefficient> MultivectorGradeIterMut for SeaFlectorPGA3D<T> {
	fn grade_iter_mut(&mut self, grade: &SeaGrade) -> <&mut Self as MultivectorGradeIntoIterator>::MultivectorGradeIntoIter {
		self.grade_into_iter(grade)
	}
}

impl<'a, T:SeaCoefficient> MultivectorGradeIntoIterator for &'a SeaFlectorPGA3D<T> {
	type MultivectorGradeIntoIter = SeaCoreSparseMultivector4DIter<'a, T>;

	fn grade_into_iter(self, grade: &SeaGrade) -> Self::MultivectorGradeIntoIter {
		self.0.grade_iter(grade)
	}
}
impl<'a, T:SeaCoefficient> MultivectorGradeIntoIterator for &'a mut SeaFlectorPGA3D<T> {
	type MultivectorGradeIntoIter = SeaCoreSparseMultivector4DIterMut<'a, T>;

	fn grade_into_iter(self, grade: &SeaGrade) -> Self::MultivectorGradeIntoIter {
		self.0.grade_iter_mut(grade)
	}
}