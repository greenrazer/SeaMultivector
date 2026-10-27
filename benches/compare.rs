//! The same operations in Seagrass and in other Rust crates, on one machine with one harness.
//!
//! Results are recorded in `BENCHMARKS.md`. Each competitor runs in its documented fastest
//! single-value form: `clifford` and `geometric_algebra` through their generated PGA types, nalgebra
//! through matrices, isometries and dual quaternions.

use std::hint::black_box;

use criterion::{Criterion, criterion_group, criterion_main};
use nalgebra::{Isometry3, Point3, Translation3, UnitDualQuaternion, UnitQuaternion, Vector4};
use seagrass_math::prelude::*;

fn motor_times_motor(c: &mut Criterion) {
    let mut group = c.benchmark_group("motor * motor");

    let (a, b) = (
        Motor::<f64, Pga3>::new_translation_motor(1.0, 2.0, 3.0)
            * Motor::new_rotation_motor(0.0, 0.0, 1.1),
        Motor::<f64, Pga3>::new_rotation_motor(0.0, 0.0, 0.3)
            * Motor::new_translation_motor(0.5, -1.0, 2.0),
    );
    group.bench_function("seagrass f64", |bench| {
        bench.iter(|| black_box(&a) * black_box(&b))
    });
    let (a32, b32) = (
        Motor::<f32, Pga3>::new_translation_motor(1.0, 2.0, 3.0)
            * Motor::new_rotation_motor(0.0, 0.0, 1.1),
        Motor::<f32, Pga3>::new_rotation_motor(0.0, 0.0, 0.3)
            * Motor::new_translation_motor(0.5, -1.0, 2.0),
    );
    group.bench_function("seagrass f32", |bench| {
        bench.iter(|| black_box(&a32) * black_box(&b32))
    });

    {
        use clifford::specialized::projective::dim3::Motor;
        let a = Motor::<f64>::from_translation(1.0, 2.0, 3.0) * Motor::<f64>::from_rotation_z(1.1);
        let b = Motor::<f64>::from_rotation_z(0.3) * Motor::<f64>::from_translation(0.5, -1.0, 2.0);
        group.bench_function("clifford f64", |bench| {
            bench.iter(|| black_box(a) * black_box(b))
        });
    }
    {
        use geometric_algebra::ppga3d::Motor;
        let a = Motor::new(0.8, 0.1, -0.2, 0.3, 0.05, 1.0, 2.0, 3.0);
        let b = Motor::new(0.7, -0.1, 0.25, 0.2, -0.05, 0.5, -1.0, 2.0);
        group.bench_function("geometric_algebra f32", |bench| {
            bench.iter(|| black_box(a) * black_box(b))
        });
    }
    let (ia, ib) = isometries();
    group.bench_function("nalgebra Isometry3 f64", |bench| {
        bench.iter(|| black_box(ia) * black_box(ib))
    });
    let (qa, qb) = (
        UnitDualQuaternion::from_isometry(&ia),
        UnitDualQuaternion::from_isometry(&ib),
    );
    group.bench_function("nalgebra UnitDualQuaternion f64", |bench| {
        bench.iter(|| black_box(qa) * black_box(qb))
    });
    let (ma, mb) = (ia.to_homogeneous(), ib.to_homogeneous());
    group.bench_function("nalgebra Matrix4 f64", |bench| {
        bench.iter(|| black_box(ma) * black_box(mb))
    });
    group.finish();
}

fn motor_applied_to_point(c: &mut Criterion) {
    let mut group = c.benchmark_group("motor applied to a point");

    let motor = Motor::<f64, Pga3>::new_translation_motor(1.0, 2.0, 3.0)
        * Motor::new_rotation_motor(0.3, -0.2, 1.1);
    let point = Multivector::<f64, Pga3>::new_point(1.0, 2.0, 3.0);
    group.bench_function("seagrass f64", |bench| {
        bench.iter(|| black_box(&motor) / black_box(&point))
    });

    {
        use clifford::ops::Transform;
        use clifford::specialized::projective::dim3::{Motor, Point};
        let m = Motor::<f64>::from_translation(1.0, 2.0, 3.0) * Motor::<f64>::from_rotation_z(1.1);
        let p = Point::<f64>::from_cartesian(1.0, 2.0, 3.0);
        group.bench_function("clifford f64", |bench| {
            bench.iter(|| black_box(&m).transform(black_box(&p)))
        });
    }
    {
        use geometric_algebra::Transformation;
        use geometric_algebra::ppga3d::{Motor, Point};
        let m = Motor::new(0.8, 0.1, -0.2, 0.3, 0.05, 1.0, 2.0, 3.0);
        let p = Point::new(1.0, 1.0, 2.0, 3.0);
        group.bench_function("geometric_algebra f32", |bench| {
            bench.iter(|| black_box(m).transformation(black_box(p)))
        });
    }
    let (isometry, _) = isometries();
    let p = Point3::new(1.0, 2.0, 3.0);
    group.bench_function("nalgebra Isometry3 f64", |bench| {
        bench.iter(|| black_box(isometry).transform_point(black_box(&p)))
    });
    let dq = UnitDualQuaternion::from_isometry(&isometry);
    group.bench_function("nalgebra UnitDualQuaternion f64", |bench| {
        bench.iter(|| black_box(dq).transform_point(black_box(&p)))
    });
    let m = isometry.to_homogeneous();
    let v = Vector4::new(1.0, 2.0, 3.0, 1.0);
    group.bench_function("nalgebra Matrix4 f64", |bench| {
        bench.iter(|| black_box(m) * black_box(v))
    });
    group.finish();
}

/// Sparse values in larger algebras, where dense storage visits every coefficient pair.
fn conformal_and_6d(c: &mut Criterion) {
    let mut group = c.benchmark_group("larger algebras");

    let (a, b) = (
        Multivector::<f64, Cga3>::from_str("e1 + 2e2 + 3e3 + 0.5e4 + 1.5e5").unwrap(),
        Multivector::<f64, Cga3>::from_str("-e1 + e2 + 0.5e3 + 2e4 - e5").unwrap(),
    );
    group.bench_function("cga3 vector * vector, seagrass f64", |bench| {
        bench.iter(|| black_box(&a) * black_box(&b))
    });
    let rotor = Multivector::<f64, Cga3>::from_str("0.5 + 0.5e12 - 0.5e13 + 0.5e23").unwrap();
    group.bench_function("cga3 rotor applied to vector, seagrass f64", |bench| {
        bench.iter(|| black_box(&rotor).sandwich(black_box(&a)))
    });
    {
        use clifford::algebra::Multivector;
        use clifford::signature::Conformal3;
        let a = Multivector::<f64, Conformal3>::vector(&[1.0, 2.0, 3.0, 0.5, 1.5]);
        let b = Multivector::<f64, Conformal3>::vector(&[-1.0, 1.0, 0.5, 2.0, -1.0]);
        group.bench_function("cga3 vector * vector, clifford dense f64", |bench| {
            bench.iter(|| black_box(&a) * black_box(&b))
        });
        let mut coefficients = [0.0; 32];
        coefficients[0] = 0.5;
        (coefficients[3], coefficients[5], coefficients[6]) = (0.5, -0.5, 0.5);
        let rotor = Multivector::<f64, Conformal3>::from_coeffs(&coefficients);
        group.bench_function(
            "cga3 rotor applied to vector, clifford dense f64",
            |bench| bench.iter(|| black_box(&rotor).sandwich(black_box(&a))),
        );
    }

    let bivector = Multivector::<f64, Euclid6>::from_str(
        "e12 - 2e13 + e23 + 0.5e14 + e24 - e34 + 3e15 + e25 + e35 - e45 + e16 + 2e26 - e36 + e46 + e56",
    )
    .unwrap();
    let vector = Multivector::<f64, Euclid6>::from_str("e1 - e2 + 2e3 + 0.5e4 + e5 - 3e6").unwrap();
    group.bench_function("6d bivector * vector, seagrass f64", |bench| {
        bench.iter(|| black_box(&bivector) * black_box(&vector))
    });
    group.finish();
}

fn isometries() -> (Isometry3<f64>, Isometry3<f64>) {
    (
        Isometry3::from_parts(
            Translation3::new(1.0, 2.0, 3.0),
            UnitQuaternion::from_euler_angles(0.0, 0.0, 1.1),
        ),
        Isometry3::from_parts(
            Translation3::new(0.5, -1.0, 2.0),
            UnitQuaternion::from_euler_angles(0.0, 0.0, 0.3),
        ),
    )
}

criterion_group!(
    benches,
    motor_times_motor,
    motor_applied_to_point,
    conformal_and_6d
);
criterion_main!(benches);
