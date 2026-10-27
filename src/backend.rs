//! Where batch products run: on the CPU, with explicit SIMD instructions with the `simd` feature,
//! or on a GPU with the `gpu` feature.

use crate::batch::{Batch, Product};
use crate::coefficient::Coefficient;
use crate::metric::Metric;

#[cfg(feature = "gpu")]
mod gpu;
#[cfg(feature = "simd")]
mod simd;

#[cfg(feature = "gpu")]
pub use gpu::{GpuBackend, GpuBatch, GpuError, GpuOptions, ShaderMode};
#[cfg(feature = "simd")]
pub use simd::SimdBackend;

/// A place to run elementwise batch products: the CPU, or a GPU.
///
/// A backend keeps batches in its own [`Buffer`](Backend::Buffer) type: host memory for the CPU,
/// device memory for a GPU. [`upload`](Backend::upload) moves a [`Batch`] in, the products
/// work buffer to buffer, so a chain of operations never leaves the device, and
/// [`download`](Backend::download) brings a result back. Products may be queued; `download` and
/// [`finish`](Backend::finish) wait for them.
///
/// # Examples
///
/// ```
/// use seagrass_math::prelude::*;
///
/// fn move_points<B: Backend<f64>>(backend: &B, motors: &Batch<f64, Pga3>, points: &Batch<f64, Pga3>) -> Batch<f64, Pga3> {
///     let (motors, points) = (backend.upload(motors), backend.upload(points));
///     let moved = backend.sandwich(&motors, &points);
///     backend.download(&moved)
/// }
///
/// let motors = Batch::from_values(Pga3::default(), &[Motor::<f64, Pga3>::new_translation_motor(1.0, 0.0, 0.0)]);
/// let points = Batch::<f64, Pga3>::from_points(&[Vector3::new(0.0, 2.0, 0.0)]);
/// let moved = move_points(&CpuBackend::default(), &motors, &points);
/// assert_eq!(moved.points_to_euclidean(), [Vector3::new(1.0, 2.0, 0.0)]);
/// ```
pub trait Backend<T: Coefficient> {
    /// A batch held by this backend.
    type Buffer<M: Metric>;

    /// Copies a batch into the backend.
    fn upload<M: Metric>(&self, batch: &Batch<T, M>) -> Self::Buffer<M>;

    /// Copies a buffer back, waiting for the work that produces it.
    fn download<M: Metric>(&self, buffer: &Self::Buffer<M>) -> Batch<T, M>;

    /// The elementwise product `op(lhs[i], rhs[i])`.
    ///
    /// # Panics
    ///
    /// Panics if the buffers differ in length or runtime metric.
    fn product<M: Metric>(
        &self,
        op: Product,
        lhs: &Self::Buffer<M>,
        rhs: &Self::Buffer<M>,
    ) -> Self::Buffer<M>;

    /// The elementwise sandwich `versors[i] * operands[i] * reverse(versors[i])`.
    ///
    /// # Panics
    ///
    /// Panics if the buffers differ in length or runtime metric.
    fn sandwich<M: Metric>(
        &self,
        versors: &Self::Buffer<M>,
        operands: &Self::Buffer<M>,
    ) -> Self::Buffer<M>;

    /// Like [`product`](Backend::product), writing into `out` and reusing its storage when it
    /// already has the result's structure, as in a loop that recomputes the same batch.
    ///
    /// # Panics
    ///
    /// Panics if `lhs` and `rhs` differ in length or runtime metric.
    fn product_into<M: Metric>(
        &self,
        op: Product,
        lhs: &Self::Buffer<M>,
        rhs: &Self::Buffer<M>,
        out: &mut Self::Buffer<M>,
    ) {
        *out = self.product(op, lhs, rhs);
    }

    /// Like [`sandwich`](Backend::sandwich), writing into `out` and reusing its storage when it
    /// already has the result's structure.
    ///
    /// # Panics
    ///
    /// Panics if `versors` and `operands` differ in length or runtime metric.
    fn sandwich_into<M: Metric>(
        &self,
        versors: &Self::Buffer<M>,
        operands: &Self::Buffer<M>,
        out: &mut Self::Buffer<M>,
    ) {
        *out = self.sandwich(versors, operands);
    }

    /// Waits until all queued products have run.
    fn finish(&self);
}

/// Runs batch products on the CPU, optionally on all cores with the `rayon` feature.
///
/// Its buffers are plain [`Batch`] values, and its products equal the single-value products bit
/// for bit.
#[derive(Debug, Clone, Copy, Default)]
pub struct CpuBackend {
    parallel: bool,
}

impl CpuBackend {
    /// A backend that runs on the calling thread.
    pub fn single_threaded() -> Self {
        Self { parallel: false }
    }

    /// A backend that splits batches across all cores. Without the `rayon` feature it runs on the
    /// calling thread.
    pub fn parallel() -> Self {
        Self { parallel: true }
    }
}

impl<T: Coefficient> Backend<T> for CpuBackend {
    type Buffer<M: Metric> = Batch<T, M>;

    fn upload<M: Metric>(&self, batch: &Batch<T, M>) -> Batch<T, M> {
        batch.clone()
    }

    fn download<M: Metric>(&self, buffer: &Batch<T, M>) -> Batch<T, M> {
        buffer.clone()
    }

    fn product<M: Metric>(&self, op: Product, lhs: &Batch<T, M>, rhs: &Batch<T, M>) -> Batch<T, M> {
        lhs.product_with(op, rhs, self.parallel)
    }

    fn sandwich<M: Metric>(&self, versors: &Batch<T, M>, operands: &Batch<T, M>) -> Batch<T, M> {
        versors.sandwich_with(operands, self.parallel)
    }

    fn product_into<M: Metric>(
        &self,
        op: Product,
        lhs: &Batch<T, M>,
        rhs: &Batch<T, M>,
        out: &mut Batch<T, M>,
    ) {
        lhs.product_into(op, rhs, out, self.parallel);
    }

    fn sandwich_into<M: Metric>(
        &self,
        versors: &Batch<T, M>,
        operands: &Batch<T, M>,
        out: &mut Batch<T, M>,
    ) {
        versors.sandwich_into(operands, out, self.parallel);
    }

    fn finish(&self) {}
}
