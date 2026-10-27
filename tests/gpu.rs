//! GPU batch products agree with the CPU. Run with `cargo test --features gpu`; without a GPU
//! (or a software adapter such as lavapipe or WARP) the tests report that they skipped.
#![cfg(feature = "gpu")]

use seagrass_math::prelude::*;
use seagrass_math::{BladeSet, Dynamic};

/// The GPU, or `None` with a message saying why the test is skipped. Every batch goes to the
/// device, however small.
fn gpu(mode: ShaderMode) -> Option<GpuBackend> {
    let options = GpuOptions {
        mode,
        min_batch: 0,
        ..GpuOptions::default()
    };
    match GpuBackend::with_options(options) {
        Ok(gpu) => Some(gpu),
        Err(err) => {
            // A skipped test must say so; there is no other channel for it.
            eprintln!("skipping GPU test: {err}");
            None
        }
    }
}

/// A deterministic generator of values in [-2, 2).
struct Values(u64);

impl Values {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }

    fn batch<M: Metric>(&mut self, metric: M, blades: M::Bits, len: usize) -> Batch<f32, M> {
        let planes = (0..blades.count_ones() as usize * len)
            .map(|_| (self.next() % 4001) as f32 / 1000.0 - 2.0)
            .collect();
        Batch::from_planes(metric, blades, len, planes).unwrap()
    }
}

/// The sum of the absolute coefficients of element `i`.
fn l1<M: Metric>(batch: &Batch<f32, M>, i: usize) -> f32 {
    let planes = batch.blades().count_ones() as usize;
    (0..planes)
        .map(|k| batch.planes()[k * batch.len() + i].abs())
        .sum()
}

/// Each coefficient agrees to `1e-6` relative to `scale(i)`, a bound on the sum of the absolute
/// values of the terms it adds up. Summing in a different order (or with fused multiply-adds, as
/// GPUs may) changes an `f32` sum by a few rounding steps of that size, however much the terms
/// cancel.
#[track_caller]
fn assert_close<M: Metric>(
    gpu: &Batch<f32, M>,
    cpu: &Batch<f32, M>,
    scale: impl Fn(usize) -> f32,
    what: &str,
) {
    assert_eq!(gpu.blades(), cpu.blades(), "{what}: stored blades");
    let planes = cpu.blades().count_ones() as usize;
    for i in 0..cpu.len() {
        let scale = scale(i).max(1.0);
        for k in 0..planes {
            let (g, c) = (
                gpu.planes()[k * cpu.len() + i],
                cpu.planes()[k * cpu.len() + i],
            );
            assert!(
                (g - c).abs() <= 1e-6 * scale,
                "{what}: element {i}, slot {k}: {g} != {c}"
            );
        }
    }
}

const PRODUCTS: [Product; 5] = [
    Product::Geometric,
    Product::Outer,
    Product::Inner,
    Product::LeftContraction,
    Product::RightContraction,
];

#[test]
fn pga3d_products_and_sandwiches_match_the_cpu() {
    for mode in [ShaderMode::Specialized, ShaderMode::Interpreted] {
        let Some(gpu) = gpu(mode) else { return };
        let mut values = Values(7);
        let pga = Pga3::default();
        let (motor, point, plane, line) = (0x9669_u16, 0x6880_u16, 0x0116_u16, 0x1668_u16);
        let len = 3_000; // more than one workgroup, not a multiple of one
        for (lhs, rhs) in [
            (motor, motor),
            (motor, point),
            (line, plane),
            (point, point),
            (plane, line),
        ] {
            let (a, b) = (values.batch(pga, lhs, len), values.batch(pga, rhs, len));
            let (ga, gb) = (gpu.upload(&a), gpu.upload(&b));
            for op in PRODUCTS {
                let what = format!("{mode:?} {op:?} {lhs:#x} x {rhs:#x}");
                let scale = |i| l1(&a, i) * l1(&b, i);
                assert_close(
                    &gpu.download(&gpu.product(op, &ga, &gb)),
                    &a.product(op, &b),
                    scale,
                    &what,
                );
            }
            let what = format!("{mode:?} sandwich {lhs:#x} / {rhs:#x}");
            let scale = |i| l1(&a, i).powi(2) * l1(&b, i);
            assert_close(
                &gpu.download(&gpu.sandwich(&ga, &gb)),
                &a.sandwich(&b),
                scale,
                &what,
            );
        }
    }
}

fn random_structures<B: BladeSet>(gpu: &GpuBackend, seed: u64) {
    let mut values = Values(seed);
    let dim = B::DIM;
    for p in (0..=dim).step_by(2) {
        let metric = Dynamic::<B>::new(p, dim - p - u8::from(p < dim), u8::from(p < dim)).unwrap();
        for _ in 0..3 {
            let (lhs, rhs) = (B::from_u64(values.next()), B::from_u64(values.next()));
            let (a, b) = (
                values.batch(metric, lhs, 700),
                values.batch(metric, rhs, 700),
            );
            let (ga, gb) = (gpu.upload(&a), gpu.upload(&b));
            for op in PRODUCTS {
                let what = format!("{op:?} in {metric} {lhs:?} x {rhs:?}");
                let scale = |i| l1(&a, i) * l1(&b, i);
                assert_close(
                    &gpu.download(&gpu.product(op, &ga, &gb)),
                    &a.product(op, &b),
                    scale,
                    &what,
                );
            }
            let what = format!("sandwich in {metric} {lhs:?} / {rhs:?}");
            let scale = |i| l1(&a, i).powi(2) * l1(&b, i);
            assert_close(
                &gpu.download(&gpu.sandwich(&ga, &gb)),
                &a.sandwich(&b),
                scale,
                &what,
            );
        }
    }
}

#[test]
fn random_structures_in_3d_to_6d_match_the_cpu() {
    for mode in [ShaderMode::Specialized, ShaderMode::Interpreted] {
        let Some(gpu) = gpu(mode) else { return };
        random_structures::<u8>(&gpu, 1);
        random_structures::<u16>(&gpu, 2);
        random_structures::<u32>(&gpu, 3);
        random_structures::<u64>(&gpu, 4);
    }
}

/// A frame of a simulation: compose every motor with a small motion, then move every point,
/// without reading anything back in between.
#[test]
fn chained_products_stay_on_the_device_and_match_the_cpu() {
    let Some(gpu) = gpu(ShaderMode::Specialized) else {
        return;
    };
    let len = 50_000;
    let pga = Pga3::default();
    let step = Motor::<f32, Pga3>::new_rotation_motor(0.0, 0.0, 0.01)
        * Motor::new_translation_motor(0.0, 0.0, 0.01);
    let steps = Batch::from_values(pga, &vec![step; len]);
    let motors = Batch::from_values(pga, &vec![Motor::<f32, Pga3>::identity(); len]);
    let points: Vec<Vector3<f32>> = (0..len)
        .map(|i| Vector3::new(i as f32 / len as f32, 1.0, 0.0))
        .collect();
    let points = Batch::<f32, Pga3>::from_points(&points);

    let (gsteps, gpoints) = (gpu.upload(&steps), gpu.upload(&points));
    let mut gmotors = gpu.upload(&motors);
    let (mut cmotors, mut cpoints) = (motors, points.clone());
    for _ in 0..10 {
        gmotors = gpu.product(Product::Geometric, &gsteps, &gmotors);
        cmotors = steps.geometric_product(&cmotors);
    }
    let gmoved = gpu.sandwich(&gmotors, &gpoints);
    cpoints = cmotors.sandwich(&cpoints);
    assert!(gmoved.is_on_device());
    // Unit motors and points with coordinates up to 1: the terms stay small.
    assert_close(
        &gpu.download(&gmoved),
        &cpoints,
        |_| 16.0,
        "ten composed steps, then a sandwich",
    );
}

#[test]
fn small_batches_stay_on_the_cpu() {
    let Ok(gpu) = GpuBackend::new() else { return };
    let points = Batch::<f32, Pga3>::from_points(&[Vector3::new(1.0, 2.0, 3.0)]);
    let uploaded = gpu.upload(&points);
    assert!(!uploaded.is_on_device());
    assert_eq!(gpu.download(&uploaded), points);
}
