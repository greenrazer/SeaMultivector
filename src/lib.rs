#![allow(dead_code)]
#![allow(unused_imports)]

mod errors;
mod metric;
mod coefficient;
mod low_dimensional_multivector;
mod multivector_traits;
mod multivector_product_traits;


pub mod prelude {
    pub use std::str::FromStr;
    pub use super::low_dimensional_multivector::{
        Multivector3D, 
        Multivector4D, 
        Multivector5D, 
        Multivector6D
    };
}


// Multivector core
// conjugations
// geometric product
// BasisBlades LeftContration, RightContraction, Inner, Outer
// LeftContration, RightContraction
// Inner

// not done:

//   Outer ^

// Mul for all products
// all products for mv with metric
//  just to see how fast implmentation goes
// test our implementation in lib.rs

// Dual !

// Add basisblade for multivector

// Regressive Product
// Norm


// grade iterator
// collect into graded multivector (grade, array<T>) dense

//scalar is T
// vector is always just array<T>
//pseudovector needs a wrapper<array<T>>
//pseudoscalar needs a wrapper<T>
//Versor/ Motor/Glide Wrappers
//invere sandwich exponential logarithm


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
