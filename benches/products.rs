//! Single-value throughput of products, unary operations, conversions and parsing.
//!
//! Results are recorded in `BENCHMARKS.md`.

use std::hint::black_box;

use criterion::{Criterion, criterion_group, criterion_main};
use seagrass_math::prelude::*;

/// 19 of the 64 blades of a 6D (3,3) algebra, mixed grades.
const SIX_19: &str = "3,3>1 + e1 - 2e2 + 0.5e12 + 3e3 - e13 + 2e23 + e123 + e4 - 0.5e14 + e24 + 2e34 \
                      + e5 - e15 + e25 + 3e6 - 2e16 + e36 + e123456";
/// 9 blades of the same algebra.
const SIX_9: &str = "3,3>2 - e1 + e2 + e12 - e34 + 0.5e5 + e56 - e246 + e12345";

fn six(s: &str) -> Multivector<f64, Dyn6> {
    Multivector::from_str(s).unwrap()
}

fn motor() -> Motor<f64, Pga3> {
    Motor::<f64, Pga3>::new_translation_motor(1.0, -2.0, 0.5)
        * Motor::<f64, Pga3>::new_rotation_motor(0.3, -0.2, 1.1)
}

fn products_6d(c: &mut Criterion) {
    let (a, b) = (six("6>e1 + 2e2"), six("6>3e5 + e6"));
    c.bench_function("6d vector * vector (2 x 2)", |bench| {
        bench.iter(|| black_box(&a) * black_box(&b))
    });

    let bivector = six(
        "6>e12 - 2e13 + e23 + 0.5e14 + e24 - e34 + 3e15 + e25 + e35 - e45 + e16 + 2e26 - e36 + e46 + e56",
    );
    let vector = six("6>e1 - e2 + 2e3 + 0.5e4 + e5 - 3e6");
    c.bench_function("6d bivector * vector (15 x 6)", |bench| {
        bench.iter(|| black_box(&bivector) * black_box(&vector))
    });

    let (x, y) = (six(SIX_19), six(SIX_9));
    c.bench_function("6d (3,3) 19 x 9 geometric", |bench| {
        bench.iter(|| black_box(&x) * black_box(&y))
    });
    c.bench_function("6d (3,3) 19 x 9 outer", |bench| {
        bench.iter(|| black_box(&x) ^ black_box(&y))
    });
    c.bench_function("6d (3,3) 19 x 9 inner", |bench| {
        bench.iter(|| black_box(&x) | black_box(&y))
    });
    c.bench_function("6d (3,3) 19 x 9 regressive", |bench| {
        bench.iter(|| black_box(&x) & black_box(&y))
    });
    c.bench_function("6d 19 dual", |bench| bench.iter(|| black_box(&x).dual()));
    c.bench_function("6d 19 reverse", |bench| {
        bench.iter(|| black_box(&x).reversed())
    });
    c.bench_function("6d 19 + 9 add", |bench| {
        bench.iter(|| black_box(&x) + black_box(&y))
    });
    c.bench_function("6d 19 scalar_norm", |bench| {
        bench.iter(|| black_box(&x).norm())
    });
}

fn pga3d(c: &mut Criterion) {
    let (m, n) = (
        motor(),
        Motor::<f64, Pga3>::new_rotation_motor(-0.4, 0.9, 0.1),
    );
    let point = Multivector::<f64, Pga3>::new_point(1.0, 2.0, 3.0);
    let other_point = Multivector::<f64, Pga3>::new_point(-2.0, 0.5, 4.0);
    let line = Multivector::<f64, Pga3>::new_line(0.3, -0.2, 1.1);
    let plane = Multivector::<f64, Pga3>::new_plane(1.0, 2.0, -1.0, 0.5);
    let matrix = Matrix4::new_translation(&Vector3::new(1.0, 2.0, 3.0))
        * Matrix4::from_euler_angles(0.3, -0.2, 1.1);

    c.bench_function("pga3d motor * motor", |bench| {
        bench.iter(|| black_box(&m) * black_box(&n))
    });
    c.bench_function("pga3d motor / point", |bench| {
        bench.iter(|| black_box(&m) / black_box(&point))
    });
    c.bench_function("pga3d point & point", |bench| {
        bench.iter(|| black_box(&point) & black_box(&other_point))
    });
    c.bench_function("pga3d line ^ plane", |bench| {
        bench.iter(|| black_box(&line) ^ black_box(&plane))
    });
    c.bench_function("pga3d exp(line)", |bench| {
        bench.iter(|| (black_box(&line) * 0.5).exp())
    });
    c.bench_function("pga3d motor log", |bench| {
        bench.iter(|| black_box(&m).log())
    });
    c.bench_function("pga3d motor to matrices", |bench| {
        bench.iter(|| black_box(&m).to_homogeneous_transformation_matrices())
    });
    c.bench_function("pga3d matrix to motor", |bench| {
        bench.iter(|| Motor::<f64, Pga3>::from(black_box(&matrix)))
    });
}

fn parsing(c: &mut Criterion) {
    let long = "3,0,1>1 + e1 + 2e2 + 3e12 + 4e3 + 5e13 + 6e23 + 7e123 + 8e4 + 9e14 + 10e24 + 11e124 + 12e34 \
                + 13e134 + 14e234 + 15e1234";
    c.bench_function("from_str 4d 16 terms", |bench| {
        bench.iter(|| Multivector::<f64, Dyn4>::from_str(black_box(long)).unwrap())
    });
    c.bench_function("from_str 3d 2 terms", |bench| {
        bench.iter(|| Multivector::<f64, Dyn3>::from_str(black_box("3>e1 + 2e12")).unwrap())
    });
}

criterion_group!(benches, products_6d, pga3d, parsing);
criterion_main!(benches);
