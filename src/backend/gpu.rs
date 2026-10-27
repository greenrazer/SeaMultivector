//! Batch products on a GPU through wgpu, as compute shaders generated from the cached plans.
//!
//! Every element of a batch is one GPU thread, and every thread runs the same plan: one shared
//! blade structure per batch means no divergence. Each output coefficient is a fixed sum that no
//! other thread writes, so nothing needs atomics.

use std::cell::{Cell, RefCell};
use std::collections::HashMap;
use std::error::Error;
use std::fmt::{self, Write as _};
use std::rc::Rc;
use std::sync::mpsc;

use wgpu::util::DeviceExt;

use super::Backend;
use crate::batch::{Batch, Product};
use crate::blade::BladeSet;
use crate::metric::Metric;
use crate::plan::{Op, ProductPlan, SandwichPlan, with_product_plan, with_sandwich_plan};

/// Threads per workgroup.
const WORKGROUP: u32 = 256;

/// The most workgroups along one dispatch dimension.
const MAX_GROUPS: u32 = 65_535;

/// How the GPU backend turns a plan into a shader.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub enum ShaderMode {
    /// One shader per plan, with every multiply-add written out. Compiling one takes about 10 ms
    /// the first time a structure pair is seen; it is then cached.
    #[default]
    Specialized,
    /// One shader for every plan, which reads the plan from a buffer. Slower (several times, in
    /// the measurements), but it never compiles anything new.
    Interpreted,
}

/// Options for [`GpuBackend::with_options`].
#[derive(Debug, Clone)]
pub struct GpuOptions {
    /// How plans become shaders.
    pub mode: ShaderMode,
    /// Batches with fewer elements stay on the CPU, where launching a GPU dispatch (about 15 µs)
    /// would cost more than the work.
    pub min_batch: usize,
    /// The wgpu backends to try, in order. `None` reads `WGPU_BACKEND` if set and otherwise tries
    /// Vulkan, Metal, DX12 and OpenGL.
    pub backends: Option<Vec<wgpu::Backends>>,
    /// Ask for the software adapter (lavapipe, WARP) instead of a GPU, for testing.
    pub force_fallback_adapter: bool,
}

impl Default for GpuOptions {
    fn default() -> Self {
        Self {
            mode: ShaderMode::Specialized,
            min_batch: 20_000,
            backends: None,
            force_fallback_adapter: false,
        }
    }
}

/// Why no GPU backend could be created.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum GpuError {
    /// No backend offered an adapter and a device; the messages say why each failed.
    Unavailable(Vec<String>),
}

impl fmt::Display for GpuError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Unavailable(reasons) => write!(f, "no GPU available: {}", reasons.join("; ")),
        }
    }
}

impl Error for GpuError {}

/// A compiled shader and the layout of its bindings.
#[derive(Debug)]
struct Kernel {
    pipeline: wgpu::ComputePipeline,
    layout: wgpu::BindGroupLayout,
}

/// What a specialized shader or an interpreter plan buffer depends on.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
struct ShaderKey {
    op: Op,
    signature: [u8; 3],
    lhs: u64,
    rhs: u64,
}

/// Timestamp queries around the passes recorded while timing.
#[derive(Debug)]
struct Timing {
    queries: wgpu::QuerySet,
    resolve: wgpu::Buffer,
    staging: wgpu::Buffer,
    passes: Cell<u32>,
}

/// Runs batch products on a GPU through wgpu, in `f32`.
///
/// Buffers live in device memory, products are queued into one command encoder, and only
/// [`download`](Backend::download) and [`finish`](Backend::finish) submit and wait, so a chain of
/// products stays on the GPU. Batches smaller than [`GpuOptions::min_batch`] are kept and
/// computed on the CPU.
///
/// # Examples
///
/// ```no_run
/// use seagrass_math::prelude::*;
///
/// let gpu = GpuBackend::new().expect("a GPU");
/// let motors = Batch::from_values(Pga3::default(), &vec![Motor::<f32, Pga3>::new_translation_motor(1.0, 0.0, 0.0); 100_000]);
/// let points = Batch::<f32, Pga3>::from_points(&vec![Vector3::new(0.0, 2.0, 0.0); 100_000]);
/// let moved = gpu.sandwich(&gpu.upload(&motors), &gpu.upload(&points));
/// assert_eq!(gpu.download(&moved).points_to_euclidean()[0], Vector3::new(1.0, 2.0, 0.0));
/// ```
pub struct GpuBackend {
    device: wgpu::Device,
    queue: wgpu::Queue,
    info: wgpu::AdapterInfo,
    options: GpuOptions,
    timestamp_period: f32,
    kernels: RefCell<HashMap<(ShaderKey, ShaderMode), Rc<Kernel>>>,
    plan_buffers: RefCell<HashMap<ShaderKey, Rc<wgpu::Buffer>>>,
    encoder: RefCell<Option<wgpu::CommandEncoder>>,
    timing: Option<Timing>,
    timing_active: Cell<bool>,
}

impl fmt::Debug for GpuBackend {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("GpuBackend")
            .field("adapter", &self.info.name)
            .field("backend", &self.info.backend)
            .finish()
    }
}

impl GpuBackend {
    /// A backend on the first GPU found, with default options.
    ///
    /// # Errors
    ///
    /// [`GpuError::Unavailable`] if no backend offers an adapter and a device.
    pub fn new() -> Result<Self, GpuError> {
        Self::with_options(GpuOptions::default())
    }

    /// A backend with the given options.
    ///
    /// # Errors
    ///
    /// [`GpuError::Unavailable`] if no backend offers an adapter and a device.
    pub fn with_options(options: GpuOptions) -> Result<Self, GpuError> {
        let candidates = match (&options.backends, wgpu::Backends::from_env()) {
            (Some(list), _) => list.clone(),
            (None, Some(from_env)) => vec![from_env],
            (None, None) => {
                vec![
                    wgpu::Backends::VULKAN,
                    wgpu::Backends::METAL,
                    wgpu::Backends::DX12,
                    wgpu::Backends::GL,
                ]
            }
        };
        let mut reasons = Vec::new();
        for backends in candidates {
            match Self::open(backends, options.force_fallback_adapter) {
                Ok((device, queue, info, timestamps)) => {
                    let timing = timestamps.then(|| Timing::new(&device));
                    let timestamp_period = queue.get_timestamp_period();
                    return Ok(Self {
                        device,
                        queue,
                        info,
                        options,
                        timestamp_period,
                        kernels: RefCell::default(),
                        plan_buffers: RefCell::default(),
                        encoder: RefCell::default(),
                        timing,
                        timing_active: Cell::new(false),
                    });
                }
                Err(reason) => reasons.push(format!("{backends:?}: {reason}")),
            }
        }
        Err(GpuError::Unavailable(reasons))
    }

    fn open(
        backends: wgpu::Backends,
        force_fallback_adapter: bool,
    ) -> Result<(wgpu::Device, wgpu::Queue, wgpu::AdapterInfo, bool), String> {
        let mut descriptor = wgpu::InstanceDescriptor::new_without_display_handle_from_env();
        descriptor.backends = backends;
        let instance = wgpu::Instance::new(descriptor);
        let adapter = pollster::block_on(instance.request_adapter(&wgpu::RequestAdapterOptions {
            power_preference: wgpu::PowerPreference::HighPerformance,
            force_fallback_adapter,
            compatible_surface: None,
            apply_limit_buckets: false,
        }))
        .map_err(|e| e.to_string())?;
        let timestamps = adapter.features().contains(wgpu::Features::TIMESTAMP_QUERY);
        let required_features = if timestamps {
            wgpu::Features::TIMESTAMP_QUERY
        } else {
            wgpu::Features::empty()
        };
        let (device, queue) = pollster::block_on(adapter.request_device(&wgpu::DeviceDescriptor {
            label: Some("seagrass"),
            required_features,
            required_limits: adapter.limits(),
            experimental_features: wgpu::ExperimentalFeatures::disabled(),
            memory_hints: wgpu::MemoryHints::Performance,
            trace: wgpu::Trace::Off,
        }))
        .map_err(|e| e.to_string())?;
        Ok((device, queue, adapter.get_info(), timestamps))
    }

    /// The adapter in use: its name, backend and driver.
    pub fn adapter_info(&self) -> &wgpu::AdapterInfo {
        &self.info
    }

    /// The options in use.
    pub fn options(&self) -> &GpuOptions {
        &self.options
    }

    /// Changes how plans become shaders, for products queued from now on.
    pub fn set_shader_mode(&mut self, mode: ShaderMode) {
        self.options.mode = mode;
    }

    /// Changes the size below which batches stay on the CPU, for uploads from now on.
    pub fn set_min_batch(&mut self, min_batch: usize) {
        self.options.min_batch = min_batch;
    }

    /// Runs `queue` (which queues products on this backend) and returns how long the GPU spent on
    /// the products it queued, in nanoseconds, from the GPU's own timestamps. `None` if the adapter
    /// has no timestamp queries.
    pub fn time_on_gpu(&self, queue: impl FnOnce(&Self)) -> Option<f64> {
        let Some(timing) = &self.timing else {
            queue(self);
            self.finish();
            return None;
        };
        self.finish();
        timing.passes.set(0);
        self.timing_active.set(true);
        queue(self);
        self.timing_active.set(false);
        if timing.passes.get() == 0 {
            return Some(0.0);
        }
        let mut encoder = self.take_encoder();
        encoder.resolve_query_set(&timing.queries, 0..2, &timing.resolve, 0);
        encoder.copy_buffer_to_buffer(&timing.resolve, 0, &timing.staging, 0, 16);
        self.queue.submit([encoder.finish()]);
        let bytes = self.read(&timing.staging);
        let ticks: &[u64] = bytemuck::cast_slice(&bytes);
        Some(ticks[1].wrapping_sub(ticks[0]) as f64 * f64::from(self.timestamp_period))
    }

    fn take_encoder(&self) -> wgpu::CommandEncoder {
        self.encoder.borrow_mut().take().unwrap_or_else(|| {
            self.device
                .create_command_encoder(&wgpu::CommandEncoderDescriptor {
                    label: Some("seagrass"),
                })
        })
    }

    fn submit_pending(&self) {
        if let Some(encoder) = self.encoder.borrow_mut().take() {
            self.queue.submit([encoder.finish()]);
        }
    }

    /// Maps a `MAP_READ` buffer that a submitted command copied into, and copies it out.
    fn read(&self, buffer: &wgpu::Buffer) -> Vec<u8> {
        let (sender, receiver) = mpsc::channel();
        buffer.map_async(wgpu::MapMode::Read, .., move |result| {
            let _ = sender.send(result);
        });
        self.device
            .poll(wgpu::PollType::wait_indefinitely())
            .expect("the device is lost");
        receiver
            .recv()
            .expect("the map callback runs during the poll")
            .expect("mapping a staging buffer succeeds");
        let bytes = buffer
            .get_mapped_range(..)
            .expect("the buffer is mapped")
            .to_vec();
        buffer.unmap();
        bytes
    }

    fn device_buffer(&self, contents: &[f32]) -> wgpu::Buffer {
        self.device
            .create_buffer_init(&wgpu::util::BufferInitDescriptor {
                label: Some("batch"),
                contents: bytemuck::cast_slice(contents),
                usage: wgpu::BufferUsages::STORAGE
                    | wgpu::BufferUsages::COPY_SRC
                    | wgpu::BufferUsages::COPY_DST,
            })
    }

    fn kernel(
        &self,
        key: ShaderKey,
        mode: ShaderMode,
        source: impl FnOnce() -> String,
        with_plan: bool,
    ) -> Rc<Kernel> {
        // All interpreted products share one shader, as do all interpreted sandwiches.
        let cache_key = match mode {
            ShaderMode::Specialized => (key, mode),
            ShaderMode::Interpreted => {
                let op = if key.op == Op::Sandwich {
                    Op::Sandwich
                } else {
                    Op::Geometric
                };
                (
                    ShaderKey {
                        op,
                        signature: [0; 3],
                        lhs: 0,
                        rhs: 0,
                    },
                    mode,
                )
            }
        };
        if let Some(kernel) = self.kernels.borrow().get(&cache_key) {
            return Rc::clone(kernel);
        }
        let kernel = Rc::new(self.compile(&source(), with_plan));
        self.kernels
            .borrow_mut()
            .insert(cache_key, Rc::clone(&kernel));
        kernel
    }

    fn compile(&self, source: &str, with_plan: bool) -> Kernel {
        let module = self
            .device
            .create_shader_module(wgpu::ShaderModuleDescriptor {
                label: Some("seagrass plan"),
                source: wgpu::ShaderSource::Wgsl(source.into()),
            });
        let storage = |binding, read_only| wgpu::BindGroupLayoutEntry {
            binding,
            visibility: wgpu::ShaderStages::COMPUTE,
            ty: wgpu::BindingType::Buffer {
                ty: wgpu::BufferBindingType::Storage { read_only },
                has_dynamic_offset: false,
                min_binding_size: None,
            },
            count: None,
        };
        let mut entries = vec![storage(0, true), storage(1, true), storage(2, false)];
        if with_plan {
            entries.push(storage(3, true));
        }
        let layout = self
            .device
            .create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
                label: Some("seagrass plan"),
                entries: &entries,
            });
        let pipeline_layout = self
            .device
            .create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
                label: Some("seagrass plan"),
                bind_group_layouts: &[Some(&layout)],
                immediate_size: 0,
            });
        let pipeline = self
            .device
            .create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
                label: Some("seagrass plan"),
                layout: Some(&pipeline_layout),
                module: &module,
                entry_point: Some("main"),
                compilation_options: Default::default(),
                cache: None,
            });
        Kernel { pipeline, layout }
    }

    fn plan_buffer(
        &self,
        key: ShaderKey,
        words: impl FnOnce() -> Vec<[u32; 4]>,
    ) -> Rc<wgpu::Buffer> {
        if let Some(buffer) = self.plan_buffers.borrow().get(&key) {
            return Rc::clone(buffer);
        }
        let buffer = Rc::new(
            self.device
                .create_buffer_init(&wgpu::util::BufferInitDescriptor {
                    label: Some("seagrass plan entries"),
                    contents: bytemuck::cast_slice(&words()),
                    usage: wgpu::BufferUsages::STORAGE,
                }),
        );
        self.plan_buffers
            .borrow_mut()
            .insert(key, Rc::clone(&buffer));
        buffer
    }

    /// Queues one dispatch of `kernel` over `len` elements.
    fn dispatch(&self, kernel: &Kernel, buffers: &[&wgpu::Buffer], len: usize) {
        let entries: Vec<wgpu::BindGroupEntry<'_>> = buffers
            .iter()
            .enumerate()
            .map(|(binding, buffer)| wgpu::BindGroupEntry {
                binding: binding as u32,
                resource: buffer.as_entire_binding(),
            })
            .collect();
        let bind_group = self.device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("seagrass plan"),
            layout: &kernel.layout,
            entries: &entries,
        });
        let groups = (len as u32).div_ceil(WORKGROUP);
        let (x, y) = (groups.min(MAX_GROUPS), groups.div_ceil(MAX_GROUPS));

        let mut encoder = self.take_encoder();
        {
            let timestamp_writes = self.timing_active.get().then(|| {
                let timing = self
                    .timing
                    .as_ref()
                    .expect("timing is only active with timestamps");
                let first = timing.passes.get() == 0;
                timing.passes.set(timing.passes.get() + 1);
                wgpu::ComputePassTimestampWrites {
                    query_set: &timing.queries,
                    beginning_of_pass_write_index: first.then_some(0),
                    end_of_pass_write_index: Some(1),
                }
            });
            let mut pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
                label: None,
                timestamp_writes,
            });
            pass.set_pipeline(&kernel.pipeline);
            pass.set_bind_group(0, &bind_group, &[]);
            pass.dispatch_workgroups(x, y, 1);
        }
        *self.encoder.borrow_mut() = Some(encoder);
    }

    /// Moves a host-side buffer to the device, for products that mix the two.
    fn resident<M: Metric>(&self, buffer: &GpuBatch<M>) -> wgpu::Buffer {
        match &buffer.data {
            Data::Device(device) => device.clone(),
            Data::Host(batch) => self.device_buffer(batch.planes()),
        }
    }

    fn output<B: BladeSet>(&self, blades: B, len: usize) -> wgpu::Buffer {
        self.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("batch"),
            size: (blades.count_ones() as usize * len * 4) as u64,
            usage: wgpu::BufferUsages::STORAGE
                | wgpu::BufferUsages::COPY_SRC
                | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        })
    }

    fn key<M: Metric>(op: Op, metric: &M, lhs: M::Bits, rhs: M::Bits) -> ShaderKey {
        ShaderKey {
            op,
            signature: [metric.positive(), metric.negative(), metric.zero()],
            lhs: lhs.to_u64(),
            rhs: rhs.to_u64(),
        }
    }
}

impl Timing {
    fn new(device: &wgpu::Device) -> Self {
        Self {
            queries: device.create_query_set(&wgpu::QuerySetDescriptor {
                label: Some("seagrass timing"),
                ty: wgpu::QueryType::Timestamp,
                count: 2,
            }),
            resolve: device.create_buffer(&wgpu::BufferDescriptor {
                label: Some("seagrass timing"),
                size: 16,
                usage: wgpu::BufferUsages::QUERY_RESOLVE | wgpu::BufferUsages::COPY_SRC,
                mapped_at_creation: false,
            }),
            staging: device.create_buffer(&wgpu::BufferDescriptor {
                label: Some("seagrass timing"),
                size: 16,
                usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST,
                mapped_at_creation: false,
            }),
            passes: Cell::new(0),
        }
    }
}

/// Where a [`GpuBatch`] keeps its coefficients.
#[derive(Debug)]
enum Data<M: Metric> {
    /// Small batches stay in host memory and run on the CPU.
    Host(Batch<f32, M>),
    /// Planes in device memory, laid out like [`Batch::planes`].
    Device(wgpu::Buffer),
}

/// A batch held by a [`GpuBackend`]: in device memory, or in host memory when it is smaller than
/// [`GpuOptions::min_batch`].
#[derive(Debug)]
pub struct GpuBatch<M: Metric> {
    metric: M,
    blades: M::Bits,
    len: usize,
    data: Data<M>,
}

impl<M: Metric> GpuBatch<M> {
    /// The number of elements.
    pub fn len(&self) -> usize {
        self.len
    }

    /// Whether the batch has no elements.
    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// The blades every element stores.
    pub fn blades(&self) -> M::Bits {
        self.blades
    }

    /// The metric.
    pub fn metric(&self) -> M {
        self.metric
    }

    /// Whether the coefficients live in device memory, rather than on the host because the batch
    /// is small.
    pub fn is_on_device(&self) -> bool {
        matches!(self.data, Data::Device(_))
    }

    fn assert_matches(&self, rhs: &Self) {
        assert_eq!(self.metric, rhs.metric, "batches have different metrics");
        assert_eq!(self.len, rhs.len, "batches have different lengths");
    }
}

impl GpuBackend {
    fn queue_product<M: Metric>(
        &self,
        op: Product,
        lhs: &GpuBatch<M>,
        rhs: &GpuBatch<M>,
        reuse: Option<&GpuBatch<M>>,
    ) -> GpuBatch<M> {
        lhs.assert_matches(rhs);
        let (metric, len) = (lhs.metric, lhs.len);
        if let (Data::Host(a), Data::Host(b)) = (&lhs.data, &rhs.data) {
            let out = a.product_with(op, b, true);
            return GpuBatch {
                metric,
                blades: out.blades(),
                len,
                data: Data::Host(out),
            };
        }
        let key = Self::key(op.op(), &metric, lhs.blades, rhs.blades);
        with_product_plan(op.op(), &metric, lhs.blades, rhs.blades, |plan| {
            let blades: M::Bits = plan.output();
            if blades.is_empty() || lhs.blades.is_empty() || rhs.blades.is_empty() {
                let out = Batch::zeros(metric, blades, len);
                return GpuBatch {
                    metric,
                    blades,
                    len,
                    data: Data::Host(out),
                };
            }
            let (a, b) = (self.resident(lhs), self.resident(rhs));
            let out = self.output_reusing(reuse, blades, len, &[&a, &b]);
            match self.options.mode {
                ShaderMode::Specialized => {
                    let kernel =
                        self.kernel(key, ShaderMode::Specialized, || product_shader(plan), false);
                    self.dispatch(&kernel, &[&a, &b, &out], len);
                }
                ShaderMode::Interpreted => {
                    let kernel =
                        self.kernel(key, ShaderMode::Interpreted, product_interpreter, true);
                    let words = self.plan_buffer(key, || product_words(plan));
                    self.dispatch(&kernel, &[&a, &b, &out, &words], len);
                }
            }
            GpuBatch {
                metric,
                blades,
                len,
                data: Data::Device(out),
            }
        })
    }

    fn queue_sandwich<M: Metric>(
        &self,
        versors: &GpuBatch<M>,
        operands: &GpuBatch<M>,
        reuse: Option<&GpuBatch<M>>,
    ) -> GpuBatch<M> {
        versors.assert_matches(operands);
        let (metric, len) = (versors.metric, versors.len);
        if let (Data::Host(m), Data::Host(x)) = (&versors.data, &operands.data) {
            let out = m.sandwich_with(x, true);
            return GpuBatch {
                metric,
                blades: out.blades(),
                len,
                data: Data::Host(out),
            };
        }
        let key = Self::key(Op::Sandwich, &metric, versors.blades, operands.blades);
        with_sandwich_plan(&metric, versors.blades, operands.blades, |plan| {
            let blades: M::Bits = plan.output();
            if blades.is_empty() {
                let out = Batch::zeros(metric, blades, len);
                return GpuBatch {
                    metric,
                    blades,
                    len,
                    data: Data::Host(out),
                };
            }
            let (m, x) = (self.resident(versors), self.resident(operands));
            let out = self.output_reusing(reuse, blades, len, &[&m, &x]);
            match self.options.mode {
                ShaderMode::Specialized => {
                    let kernel = self.kernel(
                        key,
                        ShaderMode::Specialized,
                        || sandwich_shader(plan),
                        false,
                    );
                    self.dispatch(&kernel, &[&m, &x, &out], len);
                }
                ShaderMode::Interpreted => {
                    let kernel =
                        self.kernel(key, ShaderMode::Interpreted, sandwich_interpreter, true);
                    let words = self.plan_buffer(key, || sandwich_words(plan));
                    self.dispatch(&kernel, &[&m, &x, &out, &words], len);
                }
            }
            GpuBatch {
                metric,
                blades,
                len,
                data: Data::Device(out),
            }
        })
    }

    /// The buffer of `reuse` if it has the output's structure and is not also an input of the
    /// same dispatch (a buffer cannot be read and written by one dispatch); otherwise a new one.
    fn output_reusing<M: Metric>(
        &self,
        reuse: Option<&GpuBatch<M>>,
        blades: M::Bits,
        len: usize,
        inputs: &[&wgpu::Buffer],
    ) -> wgpu::Buffer {
        match reuse {
            Some(GpuBatch {
                blades: b,
                len: l,
                data: Data::Device(buffer),
                ..
            }) if *b == blades && *l == len && !inputs.contains(&buffer) => buffer.clone(),
            _ => self.output(blades, len),
        }
    }
}

impl Backend<f32> for GpuBackend {
    type Buffer<M: Metric> = GpuBatch<M>;

    fn upload<M: Metric>(&self, batch: &Batch<f32, M>) -> GpuBatch<M> {
        let data = if batch.len() < self.options.min_batch || batch.planes().is_empty() {
            Data::Host(batch.clone())
        } else {
            Data::Device(self.device_buffer(batch.planes()))
        };
        GpuBatch {
            metric: batch.metric(),
            blades: batch.blades(),
            len: batch.len(),
            data,
        }
    }

    fn download<M: Metric>(&self, buffer: &GpuBatch<M>) -> Batch<f32, M> {
        let Data::Device(device) = &buffer.data else {
            let Data::Host(batch) = &buffer.data else {
                unreachable!()
            };
            return batch.clone();
        };
        let staging = self.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("seagrass download"),
            size: device.size(),
            usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        let mut encoder = self.take_encoder();
        encoder.copy_buffer_to_buffer(device, 0, &staging, 0, device.size());
        self.queue.submit([encoder.finish()]);
        let planes = bytemuck::cast_slice(&self.read(&staging)).to_vec();
        Batch::from_planes(buffer.metric, buffer.blades, buffer.len, planes)
            .expect("the device buffer holds every plane")
    }

    fn product<M: Metric>(&self, op: Product, lhs: &GpuBatch<M>, rhs: &GpuBatch<M>) -> GpuBatch<M> {
        self.queue_product(op, lhs, rhs, None)
    }

    fn sandwich<M: Metric>(&self, versors: &GpuBatch<M>, operands: &GpuBatch<M>) -> GpuBatch<M> {
        self.queue_sandwich(versors, operands, None)
    }

    fn product_into<M: Metric>(
        &self,
        op: Product,
        lhs: &GpuBatch<M>,
        rhs: &GpuBatch<M>,
        out: &mut GpuBatch<M>,
    ) {
        let result = self.queue_product(op, lhs, rhs, Some(&*out));
        *out = result;
    }

    fn sandwich_into<M: Metric>(
        &self,
        versors: &GpuBatch<M>,
        operands: &GpuBatch<M>,
        out: &mut GpuBatch<M>,
    ) {
        let result = self.queue_sandwich(versors, operands, Some(&*out));
        *out = result;
    }

    fn finish(&self) {
        self.submit_pending();
        self.device
            .poll(wgpu::PollType::wait_indefinitely())
            .expect("the device is lost");
    }
}

/// The bindings and thread index shared by every shader. `a`, `b` and `o` hold planes of `n`
/// elements; `n` comes from the output's length and the number of output planes.
fn shader_header(out_planes: usize, with_plan: bool) -> String {
    let mut s = String::new();
    let _ = writeln!(s, "@group(0) @binding(0) var<storage, read> a: array<f32>;");
    let _ = writeln!(s, "@group(0) @binding(1) var<storage, read> b: array<f32>;");
    let _ = writeln!(
        s,
        "@group(0) @binding(2) var<storage, read_write> o: array<f32>;"
    );
    if with_plan {
        let _ = writeln!(
            s,
            "@group(0) @binding(3) var<storage, read> plan: array<vec4<u32>>;"
        );
    }
    let _ = writeln!(s, "@compute @workgroup_size({WORKGROUP})");
    let _ = writeln!(
        s,
        "fn main(@builtin(global_invocation_id) gid: vec3<u32>, @builtin(num_workgroups) groups: vec3<u32>) {{"
    );
    let _ = writeln!(s, "    let i = gid.x + gid.y * groups.x * {WORKGROUP}u;");
    if with_plan {
        let _ = writeln!(s, "    let n = arrayLength(&o) / plan[0].y;");
    } else {
        let _ = writeln!(s, "    let n = arrayLength(&o) / {out_planes}u;");
    }
    let _ = writeln!(s, "    if (i >= n) {{ return; }}");
    s
}

/// Appends `± factors` to an expression, starting it if it is empty.
fn push_term(expression: &mut String, coefficient: i32, factors: &str) {
    let magnitude = coefficient.unsigned_abs();
    let scaled = if magnitude == 1 {
        factors.to_owned()
    } else {
        format!("{magnitude}.0 * {factors}")
    };
    match (expression.is_empty(), coefficient < 0) {
        (true, false) => expression.push_str(&scaled),
        (true, true) => expression.push_str(&format!("-{scaled}")),
        (false, false) => expression.push_str(&format!(" + {scaled}")),
        (false, true) => expression.push_str(&format!(" - {scaled}")),
    }
}

/// A shader for one product plan, with every multiply-add written out.
fn product_shader(plan: &ProductPlan) -> String {
    let outputs = plan.outputs();
    let mut s = shader_header(outputs.len(), false);
    let mut loaded = (
        std::collections::BTreeSet::new(),
        std::collections::BTreeSet::new(),
    );
    let mut body = String::new();
    for (k, entries) in outputs.enumerate() {
        let mut expression = String::new();
        for entry in entries {
            loaded.0.insert(entry.lhs);
            loaded.1.insert(entry.rhs);
            push_term(
                &mut expression,
                if entry.negative { -1 } else { 1 },
                &format!("a{} * b{}", entry.lhs, entry.rhs),
            );
        }
        if expression.is_empty() {
            expression.push_str("0.0");
        }
        let _ = writeln!(body, "    o[{k}u * n + i] = {expression};");
    }
    for slot in loaded.0 {
        let _ = writeln!(s, "    let a{slot} = a[{slot}u * n + i];");
    }
    for slot in loaded.1 {
        let _ = writeln!(s, "    let b{slot} = b[{slot}u * n + i];");
    }
    s.push_str(&body);
    s.push_str("}\n");
    s
}

/// A shader for one sandwich plan: `a` holds the versors and `b` the operands.
fn sandwich_shader(plan: &SandwichPlan) -> String {
    let outputs = plan.outputs();
    let mut s = shader_header(outputs.len(), false);
    let mut loaded = (
        std::collections::BTreeSet::new(),
        std::collections::BTreeSet::new(),
    );
    let mut body = String::new();
    for (k, terms) in outputs.enumerate() {
        let mut expression = String::new();
        for term in terms {
            loaded.0.insert(term.first);
            loaded.0.insert(term.second);
            loaded.1.insert(term.operand);
            let factors = format!("m{} * x{} * m{}", term.first, term.operand, term.second);
            push_term(&mut expression, i32::from(term.coefficient), &factors);
        }
        if expression.is_empty() {
            expression.push_str("0.0");
        }
        let _ = writeln!(body, "    o[{k}u * n + i] = {expression};");
    }
    for slot in loaded.0 {
        let _ = writeln!(s, "    let m{slot} = a[{slot}u * n + i];");
    }
    for slot in loaded.1 {
        let _ = writeln!(s, "    let x{slot} = b[{slot}u * n + i];");
    }
    s.push_str(&body);
    s.push_str("}\n");
    s
}

/// The shared shader that runs any product plan from a buffer. `plan[0]` is
/// `(entries, output planes, 0, 0)`; each entry is `(lhs slot, rhs slot, output slot, negative)`.
fn product_interpreter() -> String {
    let mut s = shader_header(0, true);
    s.push_str(
        "    var sum: array<f32, 64>;
    for (var k = 0u; k < plan[0].y; k = k + 1u) { sum[k] = 0.0; }
    for (var e = 1u; e <= plan[0].x; e = e + 1u) {
        let p = plan[e];
        let product = a[p.x * n + i] * b[p.y * n + i];
        sum[p.z] = sum[p.z] + select(product, -product, p.w != 0u);
    }
    for (var k = 0u; k < plan[0].y; k = k + 1u) { o[k * n + i] = sum[k]; }
}
",
    );
    s
}

/// The shared shader that runs any sandwich plan from a buffer. `plan[0]` is
/// `(terms, output planes, 0, 0)`; each term is
/// `(first | operand << 8 | second << 16, coefficient as i32 bits, output slot, 0)`.
fn sandwich_interpreter() -> String {
    let mut s = shader_header(0, true);
    s.push_str(
        "    var sum: array<f32, 64>;
    for (var k = 0u; k < plan[0].y; k = k + 1u) { sum[k] = 0.0; }
    for (var e = 1u; e <= plan[0].x; e = e + 1u) {
        let p = plan[e];
        let first = p.x & 255u;
        let operand = (p.x >> 8u) & 255u;
        let second = p.x >> 16u;
        let c = f32(bitcast<i32>(p.y));
        sum[p.z] = sum[p.z] + c * a[first * n + i] * b[operand * n + i] * a[second * n + i];
    }
    for (var k = 0u; k < plan[0].y; k = k + 1u) { o[k * n + i] = sum[k]; }
}
",
    );
    s
}

fn product_words(plan: &ProductPlan) -> Vec<[u32; 4]> {
    let outputs: Vec<_> = plan.outputs().collect();
    let entries: Vec<[u32; 4]> = outputs
        .iter()
        .enumerate()
        .flat_map(|(k, entries)| {
            entries.iter().map(move |e| {
                [
                    u32::from(e.lhs),
                    u32::from(e.rhs),
                    k as u32,
                    u32::from(e.negative),
                ]
            })
        })
        .collect();
    let mut words = vec![[entries.len() as u32, outputs.len() as u32, 0, 0]];
    words.extend(entries);
    words
}

fn sandwich_words(plan: &SandwichPlan) -> Vec<[u32; 4]> {
    let outputs: Vec<_> = plan.outputs().collect();
    let terms: Vec<[u32; 4]> = outputs
        .iter()
        .enumerate()
        .flat_map(|(k, terms)| {
            terms.iter().map(move |t| {
                let slots =
                    u32::from(t.first) | u32::from(t.operand) << 8 | u32::from(t.second) << 16;
                [slots, i32::from(t.coefficient) as u32, k as u32, 0]
            })
        })
        .collect();
    let mut words = vec![[terms.len() as u32, outputs.len() as u32, 0, 0]];
    words.extend(terms);
    words
}
