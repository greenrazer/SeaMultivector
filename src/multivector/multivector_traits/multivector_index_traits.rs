use super::SeaGrade;

pub trait MultivectorIter:
	where
		for<'a> &'a Self: IntoIterator
{
	fn iter(&self) -> <&Self as IntoIterator>::IntoIter;
}

pub trait MultivectorIterMut:
	where
		for<'a> &'a mut Self: IntoIterator
{
	fn iter_mut(&mut self) -> <&mut Self as IntoIterator>::IntoIter;
}

pub trait MultivectorGradeIntoIterator:
{
	type MultivectorGradeIntoIter;
	fn grade_into_iter(self, grade: &SeaGrade) -> Self::MultivectorGradeIntoIter;
}

pub trait MultivectorGradeIter:
	where
			for<'a> &'a Self: MultivectorGradeIntoIterator
{
	fn grade_iter(&self, grade: &SeaGrade) -> <&Self as MultivectorGradeIntoIterator>::MultivectorGradeIntoIter;
}

pub trait MultivectorGradeIterMut:
	where
			for<'a> &'a mut Self: MultivectorGradeIntoIterator
{
	fn grade_iter_mut(&mut self, grade: &SeaGrade) -> <&mut Self as MultivectorGradeIntoIterator>::MultivectorGradeIntoIter;
}