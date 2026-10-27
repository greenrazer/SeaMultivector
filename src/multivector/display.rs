use super::*;

use std::fmt;

impl<T: Coefficient, M: Metric> fmt::Display for Multivector<T, M> {
    /// Writes the syntax [`FromStr`](std::str::FromStr) reads, with the metric prefix:
    /// `3,0,1>1 - 2e12 + e34`. Every stored blade is written, including zeros, so parsing the
    /// output gives back an equal value; zero is written as just the prefix. A precision such as
    /// `{:.3}` applies to every coefficient.
    ///
    /// # Examples
    ///
    /// ```
    /// use seagrass_math::prelude::*;
    ///
    /// let a = Multivector::<f64, Pga3>::from_str("1 - 2e12 + e34").unwrap();
    /// assert_eq!(a.to_string(), "3,0,1>1 - 2e12 + e34");
    /// assert_eq!(a.to_string().parse::<Multivector<f64, Pga3>>().unwrap(), a);
    /// assert_eq!(format!("{:.1}", a), "3,0,1>1.0 - 2.0e12 + 1.0e34");
    /// ```
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}>", self.metric)?;
        for (n, (index, c)) in self.iter().enumerate() {
            let sign = match (n, c.is_sign_negative()) {
                (0, false) => "",
                (0, true) => "-",
                (_, false) => " + ",
                (_, true) => " - ",
            };
            f.write_str(sign)?;
            let magnitude = c.abs();
            let implicit_one =
                index != BladeIndex::SCALAR && magnitude == T::one() && f.precision().is_none();
            if !implicit_one {
                match f.precision() {
                    Some(precision) => write!(f, "{magnitude:.precision$}")?,
                    None => write!(f, "{magnitude}")?,
                }
            }
            if index != BladeIndex::SCALAR {
                write!(f, "{index}")?;
            }
        }
        Ok(())
    }
}
