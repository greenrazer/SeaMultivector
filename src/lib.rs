#![allow(dead_code)]
#![allow(unused_imports)]

mod errors;
mod coefficient;
mod multivector;

pub mod prelude {
    pub use std::str::FromStr;
    
    pub use super::multivector::multivector_traits::*;
    pub use super::multivector::multivector_product_traits::*;

    pub use super::multivector::low_dimensional::{
        SeaMultivector3D, 
        SeaMultivector4D, 
        SeaMultivector5D, 
        SeaMultivector6D
    };

    pub use super::multivector::low_dimensional::pga3d::SeaMultivectorPGA3D;
    pub use super::multivector::low_dimensional::pga3d::SeaMotorPGA3D;
    pub use super::multivector::low_dimensional::basis_blade_set::SeaBasisBladeSet;
}

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
		let mv = SeaMultivector4D::<f32>::from_str("3,1,0>e1 + e12");
		assert!(mv.is_ok());
	}

    #[test]
    fn multivector_products_works() {
		let mv_a = SeaMultivector4D::<f32>::from_str("3,1,0>e1 + e12").unwrap();
		let mv_b = SeaMultivector4D::<f32>::from_str("3,1,0>e3 + e123 + e4").unwrap();

        let out = SeaMultivector4D::<f32>::from_str("3,1,0>-e3 +e13 +e14 +e23 +e123 +e124").unwrap();
		assert_eq!(
            &mv_a * &mv_b,
            out
        );

        let mv_a = SeaMultivector4D::<f32>::from_str("3,1,0>e1 + e12").unwrap();
		let mv_b = SeaMultivector4D::<f32>::from_str("3,1,0>e3 + e123 + e4").unwrap();

        let out = SeaMultivector4D::<f32>::from_str("3,1,0>-e3 +e23").unwrap();
		assert_eq!(
            &mv_a | &mv_b,
            out
        );

        let mv_a = SeaMultivector4D::<f32>::from_str("3,1,0>e1 + e12").unwrap();
		let mv_b = SeaMultivector4D::<f32>::from_str("3,1,0>e3 + e123 + e4").unwrap();

        let out = SeaMultivector4D::<f32>::from_str("3,1,0>e13 +e14 +e123 +e124").unwrap();
		assert_eq!(
            &mv_a ^ &mv_b,
            out
        );

        let mv_a = SeaMultivector4D::<f32>::from_str("3,1,0>e1 + e12").unwrap();
		let mv_b = SeaMultivector4D::<f32>::from_str("3,1,0>e3 + e123 + e4").unwrap();

        let out = SeaMultivector4D::<f32>::from_str("3,1,0>-e3 +e23").unwrap();
		assert_eq!(
            &mv_a << &mv_b,
            out
        );

        let mv_a = SeaMultivector4D::<f32>::from_str("3,1,0>e1 + e12").unwrap();
		let mv_b = SeaMultivector4D::<f32>::from_str("3,1,0>e3 + e123 + e4").unwrap();

        let out = SeaMultivector4D::<f32>::from_str("3,1,0>").unwrap();
		assert_eq!(
            &mv_a >> &mv_b,
            out
        );
	}
}
