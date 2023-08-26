#[derive(Debug, PartialEq)]
pub enum ParseMultivectorError{
	InvalidBasisBladeIndex,
	InvalidBasisBladeCoefficient,
	InvalidParentheses,
	InvalidBasisVector,
	InvalidBasisVectorZero,
	InvalidBasisVectorOrdering,
	InvalidBasisVectorOutOfRange,
	DuplicateBasisVector,
	InvalidBasisBladeForAlgebra,
	DuplicateBasisBlade,
	MissingMetric,
	InvalidMetricDimensions,
	InvalidMetric(ParseMetricError),
}

#[derive(Debug, PartialEq)]
pub enum ParseMetricError {
	MissingNumPositive,
	InvalidNumPositive,
	InvalidNumNegitive,
	InvalidNumZero,
	InvalidMetricLength
}

#[derive(Debug, PartialEq)]
pub enum IntoVersorError {
	NonVersorMultivector
}