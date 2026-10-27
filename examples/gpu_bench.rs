//! Batch throughput on the GPU and the CPU, against other Rust crates and a hand-written shader.
//!
//! ```text
//! cargo run --release --example gpu_bench --features gpu,rayon,simd            # every size
//! cargo run --release --example gpu_bench --features gpu,rayon,simd -- --quick  # up to 100k elements
//! cargo run --release --example gpu_bench --features gpu,rayon,simd -- --cpu    # CPU rows only
//! cargo run --release --example gpu_bench --features gpu,rayon,simd -- --chart  # redraw the chart from the CSV
//! ```
//!
//! Every row reports elements per second with the data already where it is computed (GPU rows:
//! resident in device memory, 50 dispatches queued per submission). GPU rows also report the
//! GPU's own timestamps for one dispatch, the time for one operation including upload and
//! readback, and the largest difference from the f32 CPU result. Results go to stdout as markdown
//! tables, to `target/gpu_bench.csv`, and as a chart of the "motor applied to a point" rows to
//! `target/gpu_bench.svg` (the readme's `docs/gpu-motor-point.svg`). With `--cpu`, only the CPU
//! rows run, without a GPU, and go to `target/cpu_bench.csv`.

use std::fmt::Write as _;
use std::hint::black_box;
use std::time::Instant;

use rayon::prelude::*;
use seagrass_math::prelude::*;
use seagrass_math::{BladeSet, Dyn5, GpuBatch};
use wgpu::util::DeviceExt;

const SIZES: [usize; 5] = [1_000, 10_000, 100_000, 1_000_000, 4_000_000];

/// Dispatches queued per submission when timing resident GPU work.
const QUEUED: usize = 50;

/// Seconds per call of `f`: two warm-up calls, then the median of five samples, each repeating
/// `f` until it runs for at least 20 ms.
fn seconds_per_call(mut f: impl FnMut()) -> f64 {
    f();
    f();
    let mut reps = 1;
    loop {
        let start = Instant::now();
        for _ in 0..reps {
            f();
        }
        if start.elapsed().as_secs_f64() > 0.02 || reps >= 1 << 22 {
            break;
        }
        reps *= 2;
    }
    let mut samples: Vec<f64> = (0..5)
        .map(|_| {
            let start = Instant::now();
            for _ in 0..reps {
                f();
            }
            start.elapsed().as_secs_f64() / reps as f64
        })
        .collect();
    samples.sort_by(f64::total_cmp);
    samples[2]
}

/// Keeps a GPU busy for about 200 ms, so it leaves the low clocks it drops to while the CPU rows
/// run; without this, the first GPU samples after CPU work measure the clock ramp.
fn warm_up(mut f: impl FnMut()) {
    let start = Instant::now();
    while start.elapsed().as_secs_f64() < 0.2 {
        f();
    }
}

/// One measured row at one batch size.
#[derive(Debug, Clone)]
struct Measurement {
    workload: &'static str,
    implementation: String,
    n: usize,
    /// Nanoseconds per element with the data in place.
    resident: f64,
    /// Nanoseconds per element for one operation including upload and readback.
    end_to_end: Option<f64>,
    /// Nanoseconds per element of one dispatch, from GPU timestamps.
    gpu_timestamp: Option<f64>,
    /// Largest difference from the f32 CPU result, relative to the largest CPU coefficient.
    max_difference: Option<f64>,
}

/// A small deterministic generator.
struct Random(u64);

impl Random {
    fn next(&mut self) -> f32 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        (self.0 % 2_000_001) as f32 / 1_000_000.0 - 1.0
    }
}

/// Rigid motions as rotation (Euler angles) and translation.
fn motions(n: usize, seed: u64) -> Vec<([f32; 3], [f32; 3])> {
    let mut random = Random(seed);
    (0..n)
        .map(|_| {
            (
                [random.next(), random.next(), random.next()],
                [random.next(), random.next(), random.next()],
            )
        })
        .collect()
}

fn points(n: usize, seed: u64) -> Vec<Vector3<f32>> {
    let mut random = Random(seed);
    (0..n)
        .map(|_| Vector3::new(random.next(), random.next(), random.next()))
        .collect()
}

fn motor_batch(motions: &[([f32; 3], [f32; 3])]) -> Batch<f32, Pga3> {
    let motors: Vec<Motor<f32, Pga3>> = motions
        .iter()
        .map(|([rx, ry, rz], [x, y, z])| {
            Motor::new_translation_motor(*x, *y, *z) * Motor::new_rotation_motor(*rx, *ry, *rz)
        })
        .collect();
    let mut batch = Batch::zeros(Pga3::default(), 0x9669_u16, motors.len());
    for (i, motor) in motors.iter().enumerate() {
        batch
            .set(i, motor.as_multivector())
            .expect("motors are even");
    }
    batch
}

/// The same batch storing all 16 blades of PGA3D, zeros included.
fn dense(batch: &Batch<f32, Pga3>) -> Batch<f32, Pga3> {
    let mut out = Batch::zeros(Pga3::default(), u16::MAX, batch.len());
    for i in 0..batch.len() {
        out.set(i, &batch.get(i)).expect("every blade is stored");
    }
    out
}

/// The largest difference between two batches with the same blades, relative to the largest
/// coefficient of `cpu`.
fn max_difference<M: Metric>(gpu: &Batch<f32, M>, cpu: &Batch<f32, M>) -> f64 {
    let scale = cpu.planes().iter().fold(1.0_f32, |m, c| m.max(c.abs()));
    let diff = gpu
        .planes()
        .iter()
        .zip(cpu.planes())
        .fold(0.0_f32, |m, (g, c)| m.max((g - c).abs()));
    f64::from(diff / scale)
}

/// Times one queued GPU operation: resident throughput, GPU timestamps, and one operation with
/// upload and readback. `run` queues the operation from uploaded operands into `out`.
#[allow(clippy::too_many_arguments)] // The pieces of one GPU measurement.
fn time_gpu<M: Metric>(
    workload: &'static str,
    implementation: &str,
    gpu: &GpuBackend,
    lhs: &Batch<f32, M>,
    rhs: &Batch<f32, M>,
    cpu: &Batch<f32, M>,
    run: impl Fn(&GpuBackend, &GpuBatch<M>, &GpuBatch<M>, &mut GpuBatch<M>),
) -> Measurement {
    let n = lhs.len();
    let (a, b) = (gpu.upload(lhs), gpu.upload(rhs));
    let mut out = gpu.upload(&Batch::zeros(lhs.metric(), cpu.blades(), n));
    warm_up(|| {
        for _ in 0..QUEUED {
            run(gpu, &a, &b, &mut out);
        }
        gpu.finish();
    });
    let resident = seconds_per_call(|| {
        for _ in 0..QUEUED {
            run(gpu, &a, &b, &mut out);
        }
        gpu.finish();
    }) / QUEUED as f64;
    let gpu_timestamp = gpu.time_on_gpu(|gpu| run(gpu, &a, &b, &mut out));
    let max_difference = Some(max_difference(&gpu.download(&out), cpu));
    let end_to_end = seconds_per_call(|| {
        let (a, b) = (gpu.upload(lhs), gpu.upload(rhs));
        let mut out = gpu.upload(&Batch::zeros(lhs.metric(), cpu.blades(), n));
        run(gpu, &a, &b, &mut out);
        black_box(gpu.download(&out));
    });
    Measurement {
        workload,
        implementation: implementation.to_owned(),
        n,
        resident: resident * 1e9 / n as f64,
        end_to_end: Some(end_to_end * 1e9 / n as f64),
        gpu_timestamp: gpu_timestamp.map(|ns| ns / n as f64),
        max_difference,
    }
}

fn cpu_row(workload: &'static str, implementation: &str, n: usize, f: impl FnMut()) -> Measurement {
    Measurement {
        workload,
        implementation: implementation.to_owned(),
        n,
        resident: seconds_per_call(f) * 1e9 / n as f64,
        end_to_end: None,
        gpu_timestamp: None,
        max_difference: None,
    }
}

/// A Seagrass backend that keeps its batches in host memory.
#[derive(Clone, Copy)]
enum Host {
    Cpu(CpuBackend),
    Simd(SimdBackend),
}

impl Host {
    /// The host backends, with their row names.
    fn all() -> [(&'static str, Host); 4] {
        [
            (
                "Seagrass CPU, 1 thread",
                Host::Cpu(CpuBackend::single_threaded()),
            ),
            ("Seagrass CPU, rayon", Host::Cpu(CpuBackend::parallel())),
            (
                "Seagrass SIMD, 1 thread",
                Host::Simd(SimdBackend::single_threaded()),
            ),
            ("Seagrass SIMD, rayon", Host::Simd(SimdBackend::parallel())),
        ]
    }

    fn product_into<M: Metric>(
        self,
        op: Product,
        lhs: &Batch<f32, M>,
        rhs: &Batch<f32, M>,
        out: &mut Batch<f32, M>,
    ) {
        match self {
            Host::Cpu(cpu) => cpu.product_into(op, lhs, rhs, out),
            Host::Simd(simd) => simd.product_into(op, lhs, rhs, out),
        }
    }

    fn sandwich_into<M: Metric>(
        self,
        versors: &Batch<f32, M>,
        operands: &Batch<f32, M>,
        out: &mut Batch<f32, M>,
    ) {
        match self {
            Host::Cpu(cpu) => cpu.sandwich_into(versors, operands, out),
            Host::Simd(simd) => simd.sandwich_into(versors, operands, out),
        }
    }
}

/// The Seagrass rows on the CPU, with and without explicit SIMD, on one thread and with rayon,
/// each writing into its own output batch like the competitors write into preallocated vectors.
fn seagrass_cpu_rows<M: Metric>(
    workload: &'static str,
    n: usize,
    output: &Batch<f32, M>,
    run: impl Fn(Host, &mut Batch<f32, M>),
    rows: &mut Vec<Measurement>,
) {
    for (name, backend) in Host::all() {
        let mut out = output.clone();
        rows.push(cpu_row(workload, name, n, || {
            run(backend, &mut out);
            black_box(&out);
        }));
    }
}

/// A plain wgpu device for the hand-written shader, the ceiling an engine would write by hand.
struct HandWritten {
    device: wgpu::Device,
    queue: wgpu::Queue,
    pipeline: wgpu::ComputePipeline,
}

const MATRIX_TIMES_POINT: &str = "
@group(0) @binding(0) var<storage, read> m: array<mat4x4<f32>>;
@group(0) @binding(1) var<storage, read> p: array<vec4<f32>>;
@group(0) @binding(2) var<storage, read_write> o: array<vec4<f32>>;
@compute @workgroup_size(256)
fn main(@builtin(global_invocation_id) gid: vec3<u32>, @builtin(num_workgroups) groups: vec3<u32>) {
    let i = gid.x + gid.y * groups.x * 256u;
    if (i >= arrayLength(&o)) { return; }
    o[i] = m[i] * p[i];
}
";

impl HandWritten {
    fn new() -> Option<Self> {
        let mut descriptor = wgpu::InstanceDescriptor::new_without_display_handle_from_env();
        descriptor.backends = wgpu::Backends::from_env().unwrap_or(wgpu::Backends::VULKAN);
        let instance = wgpu::Instance::new(descriptor);
        let adapter = pollster::block_on(instance.request_adapter(&wgpu::RequestAdapterOptions {
            power_preference: wgpu::PowerPreference::HighPerformance,
            force_fallback_adapter: false,
            compatible_surface: None,
            apply_limit_buckets: false,
        }))
        .ok()?;
        let (device, queue) = pollster::block_on(adapter.request_device(&wgpu::DeviceDescriptor {
            label: Some("hand-written"),
            required_features: wgpu::Features::empty(),
            required_limits: adapter.limits(),
            experimental_features: wgpu::ExperimentalFeatures::disabled(),
            memory_hints: wgpu::MemoryHints::Performance,
            trace: wgpu::Trace::Off,
        }))
        .ok()?;
        let module = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("matrix times point"),
            source: wgpu::ShaderSource::Wgsl(MATRIX_TIMES_POINT.into()),
        });
        let pipeline = device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
            label: Some("matrix times point"),
            layout: None,
            module: &module,
            entry_point: Some("main"),
            compilation_options: Default::default(),
            cache: None,
        });
        Some(Self {
            device,
            queue,
            pipeline,
        })
    }

    fn buffer(&self, data: &[f32], extra: wgpu::BufferUsages) -> wgpu::Buffer {
        self.device
            .create_buffer_init(&wgpu::util::BufferInitDescriptor {
                label: None,
                contents: bytemuck::cast_slice(data),
                usage: wgpu::BufferUsages::STORAGE | extra,
            })
    }

    fn dispatch(&self, buffers: [&wgpu::Buffer; 3], n: usize, times: usize) {
        let layout = self.pipeline.get_bind_group_layout(0);
        let bind_group = self.device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: None,
            layout: &layout,
            entries: &[
                wgpu::BindGroupEntry {
                    binding: 0,
                    resource: buffers[0].as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 1,
                    resource: buffers[1].as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 2,
                    resource: buffers[2].as_entire_binding(),
                },
            ],
        });
        let groups = (n as u32).div_ceil(256);
        let mut encoder = self
            .device
            .create_command_encoder(&wgpu::CommandEncoderDescriptor { label: None });
        for _ in 0..times {
            let mut pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
                label: None,
                timestamp_writes: None,
            });
            pass.set_pipeline(&self.pipeline);
            pass.set_bind_group(0, &bind_group, &[]);
            pass.dispatch_workgroups(groups.min(65_535), groups.div_ceil(65_535), 1);
        }
        self.queue.submit([encoder.finish()]);
        self.device
            .poll(wgpu::PollType::wait_indefinitely())
            .expect("device");
    }

    fn time(&self, matrices: &[f32], points: &[f32], n: usize) -> Measurement {
        let (m, p) = (
            self.buffer(matrices, wgpu::BufferUsages::empty()),
            self.buffer(points, wgpu::BufferUsages::empty()),
        );
        let o = self.buffer(&vec![0.0; 4 * n], wgpu::BufferUsages::COPY_SRC);
        warm_up(|| self.dispatch([&m, &p, &o], n, QUEUED));
        let resident = seconds_per_call(|| self.dispatch([&m, &p, &o], n, QUEUED)) / QUEUED as f64;
        let end_to_end = seconds_per_call(|| {
            let (m, p) = (
                self.buffer(matrices, wgpu::BufferUsages::empty()),
                self.buffer(points, wgpu::BufferUsages::empty()),
            );
            let o = self.buffer(&vec![0.0; 4 * n], wgpu::BufferUsages::COPY_SRC);
            self.dispatch([&m, &p, &o], n, 1);
            let staging = self.device.create_buffer(&wgpu::BufferDescriptor {
                label: None,
                size: o.size(),
                usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST,
                mapped_at_creation: false,
            });
            let mut encoder = self
                .device
                .create_command_encoder(&wgpu::CommandEncoderDescriptor { label: None });
            encoder.copy_buffer_to_buffer(&o, 0, &staging, 0, o.size());
            self.queue.submit([encoder.finish()]);
            staging.map_async(wgpu::MapMode::Read, .., |_| {});
            self.device
                .poll(wgpu::PollType::wait_indefinitely())
                .expect("device");
            black_box(staging.get_mapped_range(..).expect("mapped").len());
        });
        Measurement {
            workload: "motor applied to a point",
            implementation: "hand-written WGSL Mat4 × Vec4".to_owned(),
            n,
            resident: resident * 1e9 / n as f64,
            end_to_end: Some(end_to_end * 1e9 / n as f64),
            gpu_timestamp: None,
            max_difference: None,
        }
    }
}

fn pga_workloads(
    n: usize,
    gpus: &[(&str, &GpuBackend)],
    hand_written: Option<&HandWritten>,
    rows: &mut Vec<Measurement>,
) {
    let (motions_a, motions_b, xyz) = (motions(n, 1), motions(n, 2), points(n, 3));
    let (motors_a, motors_b, points) = (
        motor_batch(&motions_a),
        motor_batch(&motions_b),
        Batch::from_points(&xyz),
    );
    let (dense_a, dense_b, dense_points) = (dense(&motors_a), dense(&motors_b), dense(&points));

    // Motor × motor.
    let product = motors_a.geometric_product(&motors_b);
    let dense_product = dense_a.geometric_product(&dense_b);
    let w = "motor × motor";
    for (name, gpu) in gpus {
        rows.push(time_gpu(
            w,
            &format!("Seagrass GPU, {name}"),
            gpu,
            &motors_a,
            &motors_b,
            &product,
            |g, a, b, o| g.product_into(Product::Geometric, a, b, o),
        ));
    }
    if let Some(&(_, gpu)) = gpus.first() {
        rows.push(time_gpu(
            w,
            "Seagrass GPU, dense storage",
            gpu,
            &dense_a,
            &dense_b,
            &dense_product,
            |g, a, b, o| g.product_into(Product::Geometric, a, b, o),
        ));
    }
    seagrass_cpu_rows(
        w,
        n,
        &product,
        |cpu, out| cpu.product_into(Product::Geometric, &motors_a, &motors_b, out),
        rows,
    );
    competitors_motor_times_motor(w, n, &motions_a, &motions_b, rows);

    // Motor applied to a point.
    let moved = motors_a.sandwich(&points);
    let dense_moved = dense_a.sandwich(&dense_points);
    let w = "motor applied to a point";
    for (name, gpu) in gpus {
        rows.push(time_gpu(
            w,
            &format!("Seagrass GPU, {name}"),
            gpu,
            &motors_a,
            &points,
            &moved,
            |g, m, x, o| g.sandwich_into(m, x, o),
        ));
    }
    if let Some(&(_, gpu)) = gpus.first() {
        rows.push(time_gpu(
            w,
            "Seagrass GPU, dense storage",
            gpu,
            &dense_a,
            &dense_points,
            &dense_moved,
            |g, m, x, o| g.sandwich_into(m, x, o),
        ));
    }
    seagrass_cpu_rows(
        w,
        n,
        &moved,
        |cpu, out| cpu.sandwich_into(&motors_a, &points, out),
        rows,
    );
    let matrices: Vec<glam::Mat4> = motions_a.iter().map(glam_matrix).collect();
    if let Some(hand_written) = hand_written {
        let flat: Vec<f32> = matrices.iter().flat_map(|m| m.to_cols_array()).collect();
        let homogeneous: Vec<f32> = xyz.iter().flat_map(|p| [p.x, p.y, p.z, 1.0]).collect();
        rows.push(hand_written.time(&flat, &homogeneous, n));
    }
    competitors_motor_on_point(w, n, &motions_a, &xyz, rows);

    // A frame loop: every frame composes each motor with a small motion and moves each point.
    let w = "frame loop (per frame)";
    let frames = 100;
    let step = Batch::from_values(
        Pga3::default(),
        &vec![
            Motor::<f32, Pga3>::new_rotation_motor(0.0, 0.0, 0.001)
                * Motor::new_translation_motor(0.001, 0.0, 0.0);
            n
        ],
    );
    let step = {
        let mut padded = Batch::zeros(Pga3::default(), 0x9669_u16, n);
        for i in 0..n {
            padded.set(i, &step.get(i)).expect("even");
        }
        padded
    };
    if let Some(&(_, gpu)) = gpus.first() {
        let (g_step, g_points) = (gpu.upload(&step), gpu.upload(&points));
        let mut g_motors = gpu.upload(&motors_a);
        let mut g_next = gpu.upload(&motors_a);
        let mut g_moved = gpu.upload(&moved);
        warm_up(|| {
            gpu.sandwich_into(&g_motors, &g_points, &mut g_moved);
            gpu.finish();
        });
        let resident = seconds_per_call(|| {
            for _ in 0..frames {
                gpu.product_into(Product::Geometric, &g_step, &g_motors, &mut g_next);
                std::mem::swap(&mut g_motors, &mut g_next);
                gpu.sandwich_into(&g_motors, &g_points, &mut g_moved);
            }
            gpu.finish();
        }) / f64::from(frames);
        rows.push(Measurement {
            workload: w,
            implementation: "Seagrass GPU, specialised".to_owned(),
            n,
            resident: resident * 1e9 / n as f64,
            end_to_end: None,
            gpu_timestamp: None,
            max_difference: None,
        });
    }
    for (name, cpu) in Host::all() {
        let (mut motors, mut next, mut out) = (motors_a.clone(), motors_a.clone(), moved.clone());
        rows.push(cpu_row(w, name, n, || {
            cpu.product_into(Product::Geometric, &step, &motors, &mut next);
            std::mem::swap(&mut motors, &mut next);
            cpu.sandwich_into(&motors, &points, &mut out);
            black_box(&out);
        }));
    }
    let glam_step = glam::Mat4::from_rotation_translation(
        glam::Quat::from_rotation_z(0.001),
        glam::Vec3::new(0.001, 0.0, 0.0),
    );
    let glam_points: Vec<glam::Vec4> = xyz
        .iter()
        .map(|p| glam::Vec4::new(p.x, p.y, p.z, 1.0))
        .collect();
    let mut glam_matrices = matrices.clone();
    let mut glam_out = vec![glam::Vec4::ZERO; n];
    rows.push(cpu_row(w, "glam Mat4, 1 thread", n, || {
        for ((m, p), o) in glam_matrices
            .iter_mut()
            .zip(&glam_points)
            .zip(&mut glam_out)
        {
            *m = glam_step * *m;
            *o = *m * *p;
        }
        black_box(&glam_out);
    }));
    rows.push(cpu_row(w, "glam Mat4, rayon", n, || {
        glam_matrices
            .par_iter_mut()
            .zip(&glam_points)
            .zip(&mut glam_out)
            .for_each(|((m, p), o)| {
                *m = glam_step * *m;
                *o = *m * *p;
            });
        black_box(&glam_out);
    }));
}

fn glam_matrix(([rx, ry, rz], [x, y, z]): &([f32; 3], [f32; 3])) -> glam::Mat4 {
    glam::Mat4::from_rotation_translation(
        glam::Quat::from_euler(glam::EulerRot::XYZ, *rx, *ry, *rz),
        glam::Vec3::new(*x, *y, *z),
    )
}

fn nalgebra_isometry(([rx, ry, rz], [x, y, z]): &([f32; 3], [f32; 3])) -> nalgebra::Isometry3<f32> {
    nalgebra::Isometry3::from_parts(
        nalgebra::Translation3::new(*x, *y, *z),
        nalgebra::UnitQuaternion::from_euler_angles(*rx, *ry, *rz),
    )
}

fn ga_motor(([rx, ry, rz], [x, y, z]): &([f32; 3], [f32; 3])) -> geometric_algebra::ppga3d::Motor {
    geometric_algebra::ppga3d::Motor::new(1.0, *rx, *ry, *rz, 0.0, *x, *y, *z)
}

/// Runs `f` over paired inputs into outputs, on one thread and with rayon.
fn pairwise<A: Sync, B: Sync, O: Send + Clone>(
    workload: &'static str,
    name: &str,
    a: &[A],
    b: &[B],
    f: impl Fn(&A, &B) -> O + Sync,
    rows: &mut Vec<Measurement>,
) {
    let n = a.len();
    let mut out = vec![f(&a[0], &b[0]); n];
    rows.push(cpu_row(workload, &format!("{name}, 1 thread"), n, || {
        for ((o, a), b) in out.iter_mut().zip(a).zip(b) {
            *o = f(a, b);
        }
        black_box(&out);
    }));
    rows.push(cpu_row(workload, &format!("{name}, rayon"), n, || {
        out.par_iter_mut()
            .zip(a)
            .zip(b)
            .for_each(|((o, a), b)| *o = f(a, b));
        black_box(&out);
    }));
}

fn competitors_motor_times_motor(
    w: &'static str,
    n: usize,
    a: &[([f32; 3], [f32; 3])],
    b: &[([f32; 3], [f32; 3])],
    rows: &mut Vec<Measurement>,
) {
    let (ga, gb): (Vec<_>, Vec<_>) = (
        a.iter().map(glam_matrix).collect(),
        b.iter().map(glam_matrix).collect(),
    );
    pairwise(w, "glam Mat4", &ga, &gb, |a, b| *a * *b, rows);
    let (na, nb): (Vec<_>, Vec<_>) = (
        a.iter().map(nalgebra_isometry).collect(),
        b.iter().map(nalgebra_isometry).collect(),
    );
    let (ma, mb): (Vec<_>, Vec<_>) = (
        na.iter().map(|i| i.to_homogeneous()).collect(),
        nb.iter().map(|i| i.to_homogeneous()).collect(),
    );
    pairwise(w, "nalgebra Matrix4", &ma, &mb, |a, b| a * b, rows);
    let (qa, qb): (Vec<_>, Vec<_>) = (
        na.iter()
            .map(nalgebra::UnitDualQuaternion::from_isometry)
            .collect(),
        nb.iter()
            .map(nalgebra::UnitDualQuaternion::from_isometry)
            .collect(),
    );
    pairwise(
        w,
        "nalgebra UnitDualQuaternion",
        &qa,
        &qb,
        |a, b| a * b,
        rows,
    );
    let (gma, gmb): (Vec<_>, Vec<_>) = (
        a.iter().map(ga_motor).collect(),
        b.iter().map(ga_motor).collect(),
    );
    pairwise(
        w,
        "geometric_algebra Motor (SSE)",
        &gma,
        &gmb,
        |a, b| *a * *b,
        rows,
    );
    let lanes = |m: &[([f32; 3], [f32; 3])]| -> Vec<ultraviolet::Rotor3x8> {
        m.chunks(8)
            .map(|c| {
                let angle = |axis: usize| {
                    ultraviolet::f32x8::from(std::array::from_fn::<f32, 8, _>(|l| {
                        c.get(l).map_or(0.0, |m| m.0[axis])
                    }))
                };
                ultraviolet::Rotor3x8::from_euler_angles(angle(0), angle(1), angle(2))
            })
            .collect()
    };
    let (ua, ub) = (lanes(a), lanes(b));
    let mut rotor_rows = Vec::new();
    pairwise(
        w,
        "ultraviolet Rotor3x8 (AVX2, rotation only)",
        &ua,
        &ub,
        |a, b| *a * *b,
        &mut rotor_rows,
    );
    // Each call handled 8 elements.
    rows.extend(rotor_rows.into_iter().map(|mut m| {
        m.resident /= 8.0;
        m.n = n;
        m
    }));
}

fn competitors_motor_on_point(
    w: &'static str,
    n: usize,
    motions: &[([f32; 3], [f32; 3])],
    xyz: &[Vector3<f32>],
    rows: &mut Vec<Measurement>,
) {
    let matrices: Vec<_> = motions.iter().map(glam_matrix).collect();
    let points: Vec<_> = xyz
        .iter()
        .map(|p| glam::Vec4::new(p.x, p.y, p.z, 1.0))
        .collect();
    pairwise(
        w,
        "glam Mat4 × Vec4",
        &matrices,
        &points,
        |m, p| *m * *p,
        rows,
    );
    let isometries: Vec<_> = motions.iter().map(nalgebra_isometry).collect();
    let matrices: Vec<_> = isometries.iter().map(|i| i.to_homogeneous()).collect();
    let points4: Vec<_> = xyz.iter().map(|p| p.push(1.0)).collect();
    pairwise(
        w,
        "nalgebra Matrix4 × Vector4",
        &matrices,
        &points4,
        |m, p| m * p,
        rows,
    );
    let quaternions: Vec<_> = isometries
        .iter()
        .map(nalgebra::UnitDualQuaternion::from_isometry)
        .collect();
    let points3: Vec<_> = xyz.iter().map(|p| nalgebra::Point3::from(*p)).collect();
    pairwise(
        w,
        "nalgebra UnitDualQuaternion",
        &quaternions,
        &points3,
        |q, p| q.transform_point(p),
        rows,
    );
    let motors: Vec<_> = motions.iter().map(ga_motor).collect();
    let points: Vec<_> = xyz
        .iter()
        .map(|p| geometric_algebra::ppga3d::Point::new(1.0, p.x, p.y, p.z))
        .collect();
    pairwise(
        w,
        "geometric_algebra Motor (SSE)",
        &motors,
        &points,
        |m, p| geometric_algebra::Transformation::transformation(*m, *p),
        rows,
    );
    let rotors: Vec<_> = motions
        .chunks(8)
        .map(|c| {
            let angle = |axis: usize| {
                ultraviolet::f32x8::from(std::array::from_fn::<f32, 8, _>(|l| {
                    c.get(l).map_or(0.0, |m| m.0[axis])
                }))
            };
            ultraviolet::Rotor3x8::from_euler_angles(angle(0), angle(1), angle(2))
        })
        .collect();
    let vectors: Vec<_> = xyz
        .chunks(8)
        .map(|c| {
            let lane = |axis: usize| {
                ultraviolet::f32x8::from(std::array::from_fn::<f32, 8, _>(|l| {
                    c.get(l).map_or(0.0, |p| p[axis])
                }))
            };
            ultraviolet::Vec3x8::new(lane(0), lane(1), lane(2))
        })
        .collect();
    let mut rotor_rows = Vec::new();
    pairwise(
        w,
        "ultraviolet Rotor3x8 (AVX2, rotation only)",
        &rotors,
        &vectors,
        |r, v| {
            let mut v = *v;
            r.rotate_vec(&mut v);
            v
        },
        &mut rotor_rows,
    );
    rows.extend(rotor_rows.into_iter().map(|mut m| {
        m.resident /= 8.0;
        m.n = n;
        m
    }));
}

/// Products in other algebras through the same API and generated shaders.
fn other_algebras(n: usize, gpus: &[(&str, &GpuBackend)], rows: &mut Vec<Measurement>) {
    let mut random = Random(9);
    let mut batch = |metric_blades: u64, dim: u8| -> Vec<f32> {
        (0..metric_blades.count_ones() as usize * n * usize::from(dim > 0))
            .map(|_| random.next())
            .collect()
    };

    let cga = Dyn5::new(4, 1, 0).unwrap();
    let vectors = 0b0001_0000_0000_0001_0000_0001_0001_0110_u32; // e1..e5
    let (a, b) = (
        Batch::from_planes(cga, vectors, n, batch(u64::from(vectors), 5)).unwrap(),
        Batch::from_planes(cga, vectors, n, batch(u64::from(vectors), 5)).unwrap(),
    );
    let w = "CGA (4,1) vector × vector";
    let cpu = a.geometric_product(&b);
    for (name, gpu) in gpus {
        rows.push(time_gpu(
            w,
            &format!("Seagrass GPU, {name}"),
            gpu,
            &a,
            &b,
            &cpu,
            |g, a, b, o| g.product_into(Product::Geometric, a, b, o),
        ));
    }
    seagrass_cpu_rows(
        w,
        n,
        &cpu,
        |backend, out| backend.product_into(Product::Geometric, &a, &b, out),
        rows,
    );

    let sta = Dyn4::new(1, 3, 0).unwrap();
    let (rotors, vectors) = (0x9669_u16, 0x0116_u16);
    let (m, x) = (
        Batch::from_planes(sta, rotors, n, batch(u64::from(rotors), 4)).unwrap(),
        Batch::from_planes(sta, vectors, n, batch(u64::from(vectors), 4)).unwrap(),
    );
    let w = "STA (1,3) Lorentz transformation of a vector";
    let cpu = m.sandwich(&x);
    for (name, gpu) in gpus {
        rows.push(time_gpu(
            w,
            &format!("Seagrass GPU, {name}"),
            gpu,
            &m,
            &x,
            &cpu,
            |g, m, x, o| g.sandwich_into(m, x, o),
        ));
    }
    seagrass_cpu_rows(
        w,
        n,
        &cpu,
        |backend, out| backend.sandwich_into(&m, &x, out),
        rows,
    );

    let six = Dyn6::new(6, 0, 0).unwrap();
    let (bivectors, vectors) = (<u64 as BladeSet>::GRADE[2], <u64 as BladeSet>::GRADE[1]);
    let (a, b) = (
        Batch::from_planes(six, bivectors, n, batch(bivectors, 6)).unwrap(),
        Batch::from_planes(six, vectors, n, batch(vectors, 6)).unwrap(),
    );
    let w = "6D bivector × vector";
    let cpu = a.geometric_product(&b);
    for (name, gpu) in gpus {
        rows.push(time_gpu(
            w,
            &format!("Seagrass GPU, {name}"),
            gpu,
            &a,
            &b,
            &cpu,
            |g, a, b, o| g.product_into(Product::Geometric, a, b, o),
        ));
    }
    seagrass_cpu_rows(
        w,
        n,
        &cpu,
        |backend, out| backend.product_into(Product::Geometric, &a, &b, out),
        rows,
    );
}

/// How one series looks in the chart: its row name, legend label, palette slot, whether it runs
/// on the CPU (dashed) and its marker shape.
struct Series {
    implementation: &'static str,
    label: &'static str,
    slot: usize,
    cpu: bool,
    marker: Marker,
}

#[derive(Clone, Copy)]
enum Marker {
    Circle,
    Square,
    Diamond,
    Triangle,
}

/// The "motor applied to a point" rows as an SVG line chart of elements per second against batch
/// size, both axes logarithmic, for light and dark page themes.
fn chart(rows: &[Measurement], sizes: &[usize]) -> String {
    let series = [
        Series {
            implementation: "Seagrass GPU, specialised",
            label: "Seagrass GPU, specialised shader",
            slot: 1,
            cpu: false,
            marker: Marker::Circle,
        },
        Series {
            implementation: "hand-written WGSL Mat4 × Vec4",
            label: "Hand-written WGSL Mat4 × Vec4",
            slot: 2,
            cpu: false,
            marker: Marker::Square,
        },
        Series {
            implementation: "Seagrass GPU, interpreted",
            label: "Seagrass GPU, interpreter shader",
            slot: 3,
            cpu: false,
            marker: Marker::Diamond,
        },
        Series {
            implementation: "Seagrass CPU, rayon",
            label: "Seagrass CPU, 24 threads",
            slot: 4,
            cpu: true,
            marker: Marker::Circle,
        },
        Series {
            implementation: "glam Mat4 × Vec4, rayon",
            label: "glam Mat4 × Vec4, 24 threads",
            slot: 5,
            cpu: true,
            marker: Marker::Square,
        },
        Series {
            implementation: "glam Mat4 × Vec4, 1 thread",
            label: "glam Mat4 × Vec4, 1 thread",
            slot: 6,
            cpu: true,
            marker: Marker::Triangle,
        },
    ];
    let workload = "motor applied to a point";
    let points: Vec<Vec<(f64, f64)>> = series
        .iter()
        .map(|s| {
            sizes
                .iter()
                .filter_map(|&n| {
                    rows.iter()
                        .find(|r| {
                            r.workload == workload
                                && r.implementation == s.implementation
                                && r.n == n
                        })
                        .map(|r| ((n as f64).log10(), (1e9 / r.resident).log10()))
                })
                .collect()
        })
        .collect();
    let all: Vec<f64> = points.iter().flatten().map(|p| p.1).collect();
    // From the decade below the slowest point to a fifth of a decade above the fastest.
    let (y_min, y_max) = (
        all.iter().copied().fold(f64::MAX, f64::min).floor(),
        all.iter().copied().fold(f64::MIN, f64::max) + 0.2,
    );
    let (x_min, x_max) = (
        (sizes[0] as f64).log10(),
        (sizes[sizes.len() - 1] as f64).log10(),
    );

    let (width, height) = (760.0, 460.0);
    let (left, right, top, bottom) = (72.0, 24.0, 118.0, 52.0);
    let x = |v: f64| left + (v - x_min) / (x_max - x_min) * (width - left - right);
    let y = |v: f64| height - bottom - (v - y_min) / (y_max - y_min) * (height - top - bottom);

    let mut svg = String::new();
    let _ = write!(
        svg,
        r##"<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" width="{width}" height="{height}" role="img" aria-labelledby="title desc" class="viz">
<title id="title">Motor applied to a point: elements per second against batch size</title>
<desc id="desc">Line chart, both axes logarithmic. Series: {}.</desc>
<style>
.viz {{ --surface: #fcfcfb; --text: #0b0b0b; --text-2: #52514e; --muted: #898781; --grid: #e1e0d9; --axis: #c3c2b7;
  --s1: #2a78d6; --s2: #eb6834; --s3: #1baf7a; --s4: #eda100; --s5: #e87ba4; --s6: #008300; }}
@media (prefers-color-scheme: dark) {{ .viz {{ --surface: #1a1a19; --text: #f0efec; --text-2: #c3c2b7; --muted: #898781; --grid: #2c2c2a; --axis: #383835;
  --s1: #3987e5; --s2: #d95926; --s3: #199e70; --s4: #c98500; --s5: #d55181; --s6: #008300; }} }}
text {{ font-family: system-ui, -apple-system, "Segoe UI", sans-serif; fill: var(--text-2); font-size: 12px; }}
.title {{ fill: var(--text); font-size: 16px; font-weight: 600; }}
.tick {{ fill: var(--muted); font-size: 11px; }}
.grid {{ stroke: var(--grid); stroke-width: 1; }}
.axis {{ stroke: var(--axis); stroke-width: 1; }}
.line {{ fill: none; stroke-width: 2; stroke-linejoin: round; stroke-linecap: round; }}
.cpu {{ stroke-dasharray: 6 4; }}
.mark {{ stroke: var(--surface); stroke-width: 2; }}
{}
</style>
<rect width="{width}" height="{height}" fill="var(--surface)"/>
<text class="title" x="{left}" y="28">Motor applied to a point: elements per second</text>
<text x="{left}" y="48">Data resident where it is computed. RTX 3090 (Vulkan) and Ryzen 9 5900X, f32. Higher is better.</text>
"##,
        series
            .iter()
            .map(|s| s.label)
            .collect::<Vec<_>>()
            .join(", "),
        (1..=6)
            .map(|k| format!(".s{k} {{ stroke: var(--s{k}); }} .f{k} {{ fill: var(--s{k}); }}"))
            .collect::<Vec<_>>()
            .join("\n"),
    );

    // Legend: three columns of two rows.
    for (i, s) in series.iter().enumerate() {
        let (lx, ly) = (left + (i / 2) as f64 * 230.0, 72.0 + (i % 2) as f64 * 20.0);
        let dash = if s.cpu { " cpu" } else { "" };
        let _ = writeln!(
            svg,
            r#"<line class="line s{}{dash}" x1="{lx}" y1="{ly}" x2="{}" y2="{ly}"/>"#,
            s.slot,
            lx + 28.0
        );
        svg.push_str(&marker(s.marker, lx + 14.0, ly, s.slot));
        let _ = writeln!(
            svg,
            r#"<text x="{}" y="{}">{}</text>"#,
            lx + 36.0,
            ly + 4.0,
            s.label
        );
    }

    // Grid and ticks: one line per decade of throughput, one tick per batch size.
    let mut decade = y_min;
    while decade <= y_max {
        let label = match decade as i32 {
            6 => "1 M".to_owned(),
            7 => "10 M".to_owned(),
            8 => "100 M".to_owned(),
            9 => "1 G".to_owned(),
            10 => "10 G".to_owned(),
            11 => "100 G".to_owned(),
            d => format!("1e{d}"),
        };
        let _ = writeln!(
            svg,
            r#"<line class="grid" x1="{left}" y1="{0}" x2="{1}" y2="{0}"/>"#,
            y(decade),
            width - right
        );
        let _ = writeln!(
            svg,
            r#"<text class="tick" x="{}" y="{}" text-anchor="end">{label}</text>"#,
            left - 8.0,
            y(decade) + 4.0
        );
        decade += 1.0;
    }
    let _ = writeln!(
        svg,
        r#"<line class="axis" x1="{left}" y1="{0}" x2="{1}" y2="{0}"/>"#,
        height - bottom,
        width - right
    );
    for &n in sizes {
        let label = match n {
            n if n >= 1_000_000 => format!("{}M", n / 1_000_000),
            n => format!("{}k", n / 1_000),
        };
        let px = x((n as f64).log10());
        let _ = writeln!(
            svg,
            r#"<text class="tick" x="{px}" y="{}" text-anchor="middle">{label}</text>"#,
            height - bottom + 18.0
        );
    }
    let _ = writeln!(
        svg,
        r#"<text x="{}" y="{}" text-anchor="middle">Batch size (elements)</text>"#,
        (left + width - right) / 2.0,
        height - 10.0
    );
    let _ = writeln!(
        svg,
        r#"<text x="16" y="{0}" text-anchor="middle" transform="rotate(-90 16 {0})">Elements per second</text>"#,
        (top + height - bottom) / 2.0
    );

    // Lines first, then markers on top with a surface ring.
    for (s, points) in series.iter().zip(&points) {
        let path: Vec<String> = points
            .iter()
            .map(|&(px, py)| format!("{:.1},{:.1}", x(px), y(py)))
            .collect();
        let dash = if s.cpu { " cpu" } else { "" };
        let _ = writeln!(
            svg,
            r#"<polyline class="line s{}{dash}" points="{}"/>"#,
            s.slot,
            path.join(" ")
        );
    }
    for (s, points) in series.iter().zip(&points) {
        for &(px, py) in points {
            svg.push_str(&marker(s.marker, x(px), y(py), s.slot));
        }
    }
    svg.push_str("</svg>\n");
    svg
}

/// An 8 px marker with a 2 px ring in the surface color.
fn marker(shape: Marker, cx: f64, cy: f64, slot: usize) -> String {
    match shape {
        Marker::Circle => {
            format!(r#"<circle class="mark f{slot}" cx="{cx:.1}" cy="{cy:.1}" r="4.5"/>"#) + "\n"
        }
        Marker::Square => {
            format!(
                r#"<rect class="mark f{slot}" x="{:.1}" y="{:.1}" width="8" height="8" rx="1"/>"#,
                cx - 4.0,
                cy - 4.0
            ) + "\n"
        }
        Marker::Diamond => {
            format!(
                r#"<polygon class="mark f{slot}" points="{:.1},{:.1} {:.1},{:.1} {:.1},{:.1} {:.1},{:.1}"/>"#,
                cx,
                cy - 5.5,
                cx + 5.5,
                cy,
                cx,
                cy + 5.5,
                cx - 5.5,
                cy
            ) + "\n"
        }
        Marker::Triangle => {
            format!(
                r#"<polygon class="mark f{slot}" points="{:.1},{:.1} {:.1},{:.1} {:.1},{:.1}"/>"#,
                cx,
                cy - 5.0,
                cx + 5.0,
                cy + 4.0,
                cx - 5.0,
                cy + 4.0
            ) + "\n"
        }
    }
}

/// The fields of one CSV line, where a field in double quotes may contain commas.
fn csv_fields(line: &str) -> Vec<String> {
    let (mut fields, mut field, mut quoted) = (Vec::new(), String::new(), false);
    for c in line.chars() {
        match c {
            '"' => quoted = !quoted,
            ',' if !quoted => fields.push(std::mem::take(&mut field)),
            c => field.push(c),
        }
    }
    fields.push(field);
    fields
}

/// The rows of a CSV written by an earlier run.
fn read_csv(path: &std::path::Path) -> Vec<Measurement> {
    let text = std::fs::read_to_string(path).expect("run the benchmark first, without --chart");
    let workloads = [
        "motor × motor",
        "motor applied to a point",
        "frame loop (per frame)",
        "CGA (4,1) vector × vector",
        "STA (1,3) Lorentz transformation of a vector",
        "6D bivector × vector",
    ];
    text.lines()
        .skip(1)
        .filter_map(|line| {
            let fields = csv_fields(line);
            let number = |f: &String| f.parse::<f64>().ok();
            let [
                workload,
                implementation,
                n,
                resident,
                end_to_end,
                gpu_timestamp,
                max_difference,
            ] = &fields[..]
            else {
                return None;
            };
            Some(Measurement {
                workload: workloads.into_iter().find(|w| w == workload)?,
                implementation: implementation.clone(),
                n: n.parse().ok()?,
                resident: number(resident)?,
                end_to_end: number(end_to_end),
                gpu_timestamp: number(gpu_timestamp),
                max_difference: number(max_difference),
            })
        })
        .collect()
}

fn main() {
    let csv_path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("target")
        .join("gpu_bench.csv");
    if std::env::args().any(|a| a == "--chart") {
        let rows = read_csv(&csv_path);
        let mut sizes: Vec<usize> = rows.iter().map(|r| r.n).collect();
        sizes.sort_unstable();
        sizes.dedup();
        std::fs::write(csv_path.with_extension("svg"), chart(&rows, &sizes))
            .expect("target/ is writable");
        return;
    }
    let quick = std::env::args().any(|a| a == "--quick");
    let cpu_only = std::env::args().any(|a| a == "--cpu");
    let sizes: Vec<usize> = SIZES
        .iter()
        .copied()
        .filter(|&n| !quick || n <= 100_000)
        .collect();

    let cpu = format!(
        "CPU: {} ({} threads, SIMD backend on {})",
        std::env::var("PROCESSOR_IDENTIFIER").unwrap_or_else(|_| "unknown".into()),
        rayon::current_num_threads(),
        SimdBackend::default().instruction_set(),
    );
    let gpu = |mode| {
        GpuBackend::with_options(GpuOptions {
            min_batch: 0,
            mode,
            ..GpuOptions::default()
        })
        .expect("a GPU")
    };
    let (specialised, interpreted, hand_written);
    let gpus: Vec<(&str, &GpuBackend)> = if cpu_only {
        hand_written = None;
        println!("# CPU benchmark\n\n{cpu}. Release build, f32.\n");
        Vec::new()
    } else {
        (specialised, interpreted) = (gpu(ShaderMode::Specialized), gpu(ShaderMode::Interpreted));
        hand_written = HandWritten::new();
        let info = specialised.adapter_info();
        println!("# GPU benchmark\n");
        println!(
            "{cpu}. GPU: {} via {:?}, driver {} {}. Release build, f32.\n",
            info.name, info.backend, info.driver, info.driver_info
        );
        vec![("specialised", &specialised), ("interpreted", &interpreted)]
    };
    let mut rows = Vec::new();
    for &n in &sizes {
        eprintln!("n = {n}");
        pga_workloads(n, &gpus, hand_written.as_ref(), &mut rows);
        other_algebras(n, &gpus, &mut rows);
    }

    // Tables: million elements per second with data in place, one column per batch size.
    let mut workloads: Vec<&str> = rows.iter().map(|r| r.workload).collect();
    workloads.dedup();
    let mut seen = std::collections::HashSet::new();
    workloads.retain(|w| seen.insert(*w));
    for workload in workloads {
        println!("## {workload}: million elements per second, data in place\n");
        let header: String = sizes.iter().map(|n| format!(" {n} |")).collect();
        println!("| implementation |{header}");
        println!("|---|{}", "---:|".repeat(sizes.len()));
        let mut names: Vec<&str> = Vec::new();
        for r in rows.iter().filter(|r| r.workload == workload) {
            if !names.contains(&r.implementation.as_str()) {
                names.push(&r.implementation);
            }
        }
        for name in names {
            let mut line = format!("| {name} |");
            for &n in &sizes {
                match rows
                    .iter()
                    .find(|r| r.workload == workload && r.implementation == name && r.n == n)
                {
                    Some(r) => {
                        let _ = write!(line, " {:.0} |", 1e3 / r.resident);
                    }
                    None => line.push_str(" |"),
                }
            }
            println!("{line}");
        }
        println!();
    }

    let csv_path = if cpu_only {
        csv_path.with_file_name("cpu_bench.csv")
    } else {
        csv_path
    };
    if !cpu_only {
        println!("## GPU details at the largest size\n");
        println!(
            "| workload | implementation | GPU timestamp (ns/element) | with upload and readback (ns/element) | max difference from CPU |"
        );
        println!("|---|---|---:|---:|---:|");
        let largest = *sizes.last().expect("a size");
        for r in rows
            .iter()
            .filter(|r| r.n == largest && r.end_to_end.is_some())
        {
            println!(
                "| {} | {} | {} | {:.2} | {} |",
                r.workload,
                r.implementation,
                r.gpu_timestamp.map_or("-".into(), |t| format!("{t:.3}")),
                r.end_to_end.unwrap_or(f64::NAN),
                r.max_difference.map_or("-".into(), |d| format!("{d:.1e}"))
            );
        }
    }

    let mut csv = String::from(
        "workload,implementation,n,resident_ns_per_element,end_to_end_ns_per_element,gpu_timestamp_ns_per_element,max_difference\n",
    );
    for r in &rows {
        let _ = writeln!(
            csv,
            "\"{}\",\"{}\",{},{},{},{},{}",
            r.workload,
            r.implementation,
            r.n,
            r.resident,
            r.end_to_end.map_or(String::new(), |v| v.to_string()),
            r.gpu_timestamp.map_or(String::new(), |v| v.to_string()),
            r.max_difference.map_or(String::new(), |v| v.to_string()),
        );
    }
    let path = csv_path;
    std::fs::write(&path, csv).expect("target/ is writable");
    eprintln!("wrote {}", path.display());
    if cpu_only {
        return;
    }
    let path = path.with_extension("svg");
    std::fs::write(&path, chart(&rows, &sizes)).expect("target/ is writable");
    eprintln!("wrote {}", path.display());
}
