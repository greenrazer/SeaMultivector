//! Code written against 0.1 keeps compiling through the deprecated aliases and method names.
#![allow(deprecated)] // The point of this file is to use the deprecated names.

use seagrass_math::prelude::*;

#[test]
fn runtime_metric_aliases_and_renamed_methods_still_work() {
    let a = SeaMultivector4D::<f32>::from_str("3,1,0>e1 + e12").unwrap();
    let b = SeaMultivector4D::<f32>::from_str("3,1,0>e3 + e123 + e4").unwrap();
    assert_eq!(
        &a * &b,
        SeaMultivector4D::<f32>::from_str("3,1,0>-e3 + e13 + e14 + e23 + e123 + e124").unwrap()
    );
    assert_eq!(a.reversion_conjugated(), a.reversed());
    assert_eq!(a.clifford_conjugated(), a.conjugated());
    assert_eq!(a.right_complimented(), a.right_complemented());
    assert_eq!(a.scalar_norm(), a.norm());
    assert_eq!(a.basis_blade_coefficient(&SeaBasisBladeIndexLow::E12), 1.0);

    let mut c = a.clone();
    c.reversion_conjugation();
    assert_eq!(c, a.reversed());
}

#[test]
fn pga3d_aliases_still_work() {
    let motor = SeaMotorPGA3D::new_translation_motor(1.0, 2.0, 3.0);
    let point: SeaFlectorPGA3D<f64> = SeaMultivectorPGA3D::new_point(0.0, 0.0, 0.0);
    let moved = &motor / &point;
    assert_eq!(moved.point_to_euclidean(), Vector3::new(1.0, 2.0, 3.0));
    assert_eq!(motor.inverse().reversion_conjugated(), motor);
}

#[test]
fn error_aliases_name_the_new_types() {
    let err: ParseMultivectorError = SeaMultivector3D::<f64>::from_str("e1").unwrap_err();
    assert_eq!(err, ParseError::MissingMetric);
    let err: ParseMetricError = "a".parse::<Dyn3>().unwrap_err();
    assert_eq!(err, MetricError::InvalidNumPositive);
}
