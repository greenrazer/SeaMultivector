use super::*;

use std::str::FromStr;

use crate::error::ParseError;

impl<T: Coefficient, M: Metric> FromStr for Multivector<T, M> {
    type Err = ParseError;

    /// Reads `[p,q,r>] term (+|- term)*`, where a term is a number, `e` followed by ascending
    /// basis vector numbers, or a number directly followed by such a blade: `2`, `-e13`,
    /// `+ 0.5e124`. An empty right-hand side is zero.
    ///
    /// The `p,q,r>` prefix is required for [`Dynamic`](crate::Dynamic) metrics and optional for
    /// [`Signature`](crate::Signature) metrics. Numbers have no exponent notation, because `1e5`
    /// means 1 times the blade e5.
    ///
    /// # Errors
    ///
    /// A [`ParseError`] naming the first problem: a bad prefix, a malformed blade, a blade outside
    /// the algebra, a repeated blade, or an unexpected character.
    fn from_str(s: &str) -> Result<Self, ParseError> {
        let (prefix, body, offset) = match s.split_once('>') {
            Some((prefix, body)) => (Some(prefix), body, prefix.len() + 1),
            None => (None, s, 0),
        };
        let metric = M::from_prefix(prefix)?;

        let mut dense = [T::zero(); 64];
        let mut stored = M::Bits::EMPTY;
        parse_terms(body, offset, |index, value| {
            if u32::from(index.bits()) >= M::Bits::BLADES {
                return Err(ParseError::InvalidBasisBladeForAlgebra);
            }
            if stored.has(index) {
                return Err(ParseError::DuplicateBasisBlade);
            }
            stored |= M::Bits::from_index(index);
            dense[index.bits() as usize] = value;
            Ok(())
        })?;
        let coefficients = stored
            .iter()
            .map(|index| dense[index.bits() as usize])
            .collect();
        Ok(Self::from_raw(metric, stored, coefficients))
    }
}

/// Scans the terms of `body` and hands each blade and signed coefficient to `term`. `offset` is
/// the position of `body` in the full input, for error positions.
fn parse_terms<T: Coefficient>(
    body: &str,
    offset: usize,
    mut term: impl FnMut(BladeIndex, T) -> Result<(), ParseError>,
) -> Result<(), ParseError> {
    let bytes = body.as_bytes();
    let unexpected = |i: usize| ParseError::UnexpectedToken {
        position: offset + i,
    };
    let skip_space = |mut i: usize| {
        while i < bytes.len() && bytes[i].is_ascii_whitespace() {
            i += 1;
        }
        i
    };

    let mut i = skip_space(0);
    let mut first = true;
    while i < bytes.len() {
        // Every term but the first needs a sign.
        let negative = match bytes[i] {
            sign @ (b'+' | b'-') => {
                i = skip_space(i + 1);
                sign == b'-'
            }
            _ if first => false,
            _ => return Err(unexpected(i)),
        };
        first = false;

        let number_start = i;
        while i < bytes.len() && (bytes[i].is_ascii_digit() || bytes[i] == b'.') {
            i += 1;
        }
        let number = &body[number_start..i];

        let index = if i < bytes.len() && bytes[i] == b'e' {
            let digits_start = i + 1;
            i = digits_start;
            while i < bytes.len() && bytes[i].is_ascii_alphanumeric() {
                i += 1;
            }
            BladeIndex::from_digits(&body[digits_start..i])?
        } else if number.is_empty() {
            return Err(unexpected(i));
        } else {
            BladeIndex::SCALAR
        };

        let magnitude = if number.is_empty() {
            T::one()
        } else {
            number
                .parse::<T>()
                .map_err(|_| ParseError::InvalidBasisBladeCoefficient)?
        };
        term(index, if negative { -magnitude } else { magnitude })?;
        i = skip_space(i);
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::metric::{Dyn3, Pga3};

    fn terms(body: &str) -> Result<Vec<(u8, f64)>, ParseError> {
        let mut out = Vec::new();
        parse_terms(body, 0, |index, value| {
            out.push((index.bits(), value));
            Ok(())
        })?;
        Ok(out)
    }

    #[test]
    fn reads_signs_numbers_and_blades() {
        assert_eq!(terms(""), Ok(vec![]));
        assert_eq!(terms("  3 "), Ok(vec![(0, 3.0)]));
        assert_eq!(terms("-e1"), Ok(vec![(1, -1.0)]));
        assert_eq!(
            terms("1 - 2.5e12 +e3"),
            Ok(vec![(0, 1.0), (3, -2.5), (4, 1.0)])
        );
        assert_eq!(terms("+.5e123"), Ok(vec![(7, 0.5)]));
        // `1e5` is 1 times e5, not scientific notation.
        assert_eq!(terms("1e5"), Ok(vec![(16, 1.0)]));
    }

    #[test]
    fn rejects_input_the_grammar_does_not_cover() {
        let cases = [
            ("e1 + * e2", ParseError::UnexpectedToken { position: 5 }),
            ("e1 e2", ParseError::UnexpectedToken { position: 3 }),
            ("e1 + ", ParseError::UnexpectedToken { position: 5 }),
            ("+ - e1", ParseError::UnexpectedToken { position: 2 }),
            ("1.2.3", ParseError::InvalidBasisBladeCoefficient),
            ("2e", ParseError::InvalidBasisBladeIndex),
            ("ex", ParseError::InvalidBasisVector),
        ];
        for (body, err) in cases {
            assert_eq!(terms(body), Err(err), "{body:?}");
        }
    }

    #[test]
    fn checks_the_algebra_and_the_prefix() {
        assert_eq!(
            "3>e4".parse::<Multivector<f64, Dyn3>>(),
            Err(ParseError::InvalidBasisBladeForAlgebra)
        );
        assert_eq!(
            "3>e1 + e1".parse::<Multivector<f64, Dyn3>>(),
            Err(ParseError::DuplicateBasisBlade)
        );
        assert_eq!(
            "e1".parse::<Multivector<f64, Dyn3>>(),
            Err(ParseError::MissingMetric)
        );
        assert!("e14".parse::<Multivector<f64, Pga3>>().is_ok());
        assert_eq!(
            "3>e1".parse::<Multivector<f64, Pga3>>(),
            Err(ParseError::MismatchedMetric)
        );
    }

    #[test]
    fn reports_positions_in_the_full_input() {
        assert_eq!(
            "3,0,1>e1 ? e2".parse::<Multivector<f64, Pga3>>(),
            Err(ParseError::UnexpectedToken { position: 9 })
        );
    }
}
