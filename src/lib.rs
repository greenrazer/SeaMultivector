#![allow(dead_code)]
#![allow(unused_imports)]

mod errors;
mod metric;
mod coefficient;
mod low_dimensional_multivector;
mod multivector_traits;
mod multivector_product_traits;
mod grade;

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
// Outer ^
// Mul for all products
// Some Documentation
// Dual !
// grade iterator
// Add basisblade for multivector
// Regressive Product
// Norm

// not done:
// exp

//Versor/ Motor/Glide Wrappers
//invere sandwich exponential logarithm
//operations of bundles of multivectors with the same metric and basis blade bits


#[cfg(test)]
mod tests {
    use super::prelude::*;

    pub fn catch_unwind_silent<F: FnOnce() -> R + std::panic::UnwindSafe, R>(f: F) -> std::thread::Result<R> {
		let prev_hook = std::panic::take_hook();
		std::panic::set_hook(Box::new(|_| {}));
		let result = std::panic::catch_unwind(f);
		std::panic::set_hook(prev_hook);
		result
	}

    #[test]
    fn multivector_type_wrapper_works() {
		let mv = Multivector4D::<f32>::from_str("3,1,0>e1 + e12");
		assert!(mv.is_ok());
	}

    #[test]
    fn multivector_products_works() {
		let mv_a = Multivector4D::<f32>::from_str("3,1,0>e1 + e12").unwrap();
		let mv_b = Multivector4D::<f32>::from_str("3,1,0>e3 + e123 + e4").unwrap();

        let out = Multivector4D::<f32>::from_str("3,1,0>-e3 +e13 +e14 +e23 +e123 +e124").unwrap();
		assert_eq!(
            &mv_a * &mv_b,
            out
        );

        let mv_a = Multivector4D::<f32>::from_str("3,1,0>e1 + e12").unwrap();
		let mv_b = Multivector4D::<f32>::from_str("3,1,0>e3 + e123 + e4").unwrap();

        let out = Multivector4D::<f32>::from_str("3,1,0>-e3 +e23").unwrap();
		assert_eq!(
            &mv_a | &mv_b,
            out
        );

        let mv_a = Multivector4D::<f32>::from_str("3,1,0>e1 + e12").unwrap();
		let mv_b = Multivector4D::<f32>::from_str("3,1,0>e3 + e123 + e4").unwrap();

        let out = Multivector4D::<f32>::from_str("3,1,0>e13 +e14 +e123 +e124").unwrap();
		assert_eq!(
            &mv_a ^ &mv_b,
            out
        );

        let mv_a = Multivector4D::<f32>::from_str("3,1,0>e1 + e12").unwrap();
		let mv_b = Multivector4D::<f32>::from_str("3,1,0>e3 + e123 + e4").unwrap();

        let out = Multivector4D::<f32>::from_str("3,1,0>-e3 +e23").unwrap();
		assert_eq!(
            &mv_a << &mv_b,
            out
        );

        let mv_a = Multivector4D::<f32>::from_str("3,1,0>e1 + e12").unwrap();
		let mv_b = Multivector4D::<f32>::from_str("3,1,0>e3 + e123 + e4").unwrap();

        let out = Multivector4D::<f32>::from_str("3,1,0>").unwrap();
		assert_eq!(
            &mv_a >> &mv_b,
            out
        );
	}
}
