#![allow(dead_code)]
#![allow(unused_imports)]

mod errors;
mod metric;
mod coefficient;
mod low_dimensional_multivector;

pub mod prelude {
    pub use std::str::FromStr;
    pub use super::low_dimensional_multivector::{
        Multivector3D, 
        Multivector4D, 
        Multivector5D, 
        Multivector6D
    };
}

#[cfg(test)]
mod tests {
    use super::prelude::*;

    #[test]
    fn multivector_type_wrapper_works() {
		let mv = Multivector4D::<f32>::from_str("3,1,0>e1 + e12");
		assert!(mv.is_ok());

		let mv = Multivector4D::<f32>::from_str("3,1,0>e1 + e12");
		assert!(mv.is_ok());
	}
}
