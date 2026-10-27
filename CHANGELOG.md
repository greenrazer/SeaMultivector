# Changelog

All notable changes to this crate are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the crate follows
[Semantic Versioning](https://semver.org/) (before 1.0, a breaking change bumps the minor version).

## [Unreleased]

The 0.2.0 release: one generic engine instead of six copies of the operator surface, cached
product plans, and an API that names every type it uses. Code written against 0.1 keeps compiling
through deprecated aliases, except where noted under "Removed" and "Changed".

### Added

- `Multivector<T, M>`, one type for every algebra, generic over the metric `M`:
  `Signature<P, Q, R>` fixes the signature at compile time (zero-sized; mixing two signatures does
  not compile) and `Dynamic<B>` chooses it at runtime. Aliases `Euclid3`..`Euclid6`, `Pga2`,
  `Pga3`, `Pga4`, `Cga2`, `Cga3`, `Sta` and `Dyn3`..`Dyn6`.
- `Typed<T, M, K>` for values of known grade parity, with the kinds `Even`, `Odd` and `Any` and
  the aliases `Motor`, `Flector` and `Versor`. Product kinds follow from the operands' kinds
  (`KindMul`), checked at compile time, and work in every metric, not only PGA3D.
- Numeric constructors that skip the parser: `zero`, `from_scalar`, `basis`, `from_blades` and
  `from_parts` (which checks that the coefficient count matches the bitset).
- `Display` for `Multivector` and `Typed`, in the syntax `FromStr` reads, including the metric
  prefix; parsing the output gives back an equal value. A precision such as `{:.3}` applies to
  every coefficient.
- Accessors `metric`, `blades`, `coefficients`, `coefficients_mut`, `into_parts`, `get` and
  `pseudoscalar`, and `IntoIterator` for `&Multivector` and `&mut Multivector`.
- `Typed::new`/`TryFrom<Multivector>` with a parity check, `into_inner`, `as_multivector`,
  `apply`, and `Multivector::sandwich`.
- Exports that were unreachable: `BladeIndex` (with a constant for every blade, `E1` to
  `E123456`), `BladeSet`, `MotorCreationOrder`, `Motor::to_dual_quaternion`, and the error types.
- Error types `ParseError`, `MetricError`, `KindError` and `BladeError` implement `Display` and
  `std::error::Error`. Parse errors report `UnexpectedToken { position }` and `MismatchedMetric`.
- `nalgebra` feature, on by default, for the PGA3D matrix, vector and dual quaternion
  conversions.
- `Batch<T, M>`: many multivectors that share one blade set, stored as one array per blade.
  Elementwise `product` (any `Product`), `geometric_product` and `sandwich` run each cached plan
  once for the whole batch. Conversions `from_values`/`to_values`, `get`/`set`, `plane`, and for
  PGA3D `from_points`, `points_to_euclidean`, `from_matrices` and `motors_to_matrices`.
- The `Backend` trait (`upload`, `product`, `sandwich`, `product_into`, `sandwich_into`,
  `download`, `finish`) and `CpuBackend`, single-threaded or, with the `rayon` feature, across
  all cores.
- `gpu` feature: `GpuBackend` runs batch products in `f32` through wgpu, as WGSL generated from
  each plan and cached per structure pair, or through one plan-interpreter shader
  (`ShaderMode`). Products are queued into one submission and results stay on the device
  (`GpuBatch`) until downloaded. Batches below `GpuOptions::min_batch` (20,000 by default) run on
  the CPU. Vulkan is tried first, then Metal, DX12 and OpenGL, or the `WGPU_BACKEND` variable's
  choice.
- `simd` feature: `SimdBackend` runs batch products for `f32` and `f64` with explicit SIMD
  instructions through fearless_simd, on the widest instruction set the CPU supports when the
  program runs (AVX-512, AVX2 with FMA, SSE4.2 or SSE2 on x86, NEON on ARM, SIMD128 on
  WebAssembly), on one thread or, with `rayon`, across all cores. Results agree with
  `CpuBackend` to within rounding rather than bit for bit; other coefficient types run the
  `CpuBackend` loops.
- `examples/four_algebras.rs`: one function rotating, translating, boosting and dilating in four
  algebras.
- `examples/gpu_bench.rs`: batch throughput against glam, nalgebra, geometric_algebra,
  ultraviolet and a hand-written WGSL shader, from 1,000 to 4 million elements, with `--cpu` for
  the CPU and SIMD rows alone.
- Benchmarks (`benches/products.rs`, `benches/compare.rs`) and `BENCHMARKS.md`, which opens with
  tables comparing Seagrass with each library as ratios.
- `ARCHITECTURE.md`: how products find their pairs with bit operations, how plans are cached, and
  how batches run on each backend. The readme's "How it works" section moved there, and its
  performance sections shrank to one summary table, with the detail in `BENCHMARKS.md`.

### Changed

- Products run cached plans: the bitset structure work runs once per pair of operand blade sets,
  per thread. PGA3D motor × motor takes 51 ns instead of 239 ns, a motor applied to a point 57 ns
  instead of 511 ns, a 6D bivector × vector 123 ns instead of 697 ns (see `BENCHMARKS.md`).
- Coefficients are stored inline for values with up to 16 blades, so 3D and 4D values never
  allocate.
- The parser is hand-written: 10 to 20 times faster than the regex it replaces. Input the regex
  skipped silently, such as `e1 + * e2` or a stray character, is now an error.
- Method names follow the crate's code standards: `reversed`/`reverse`,
  `involuted`/`involute`, `conjugated`/`conjugate`, `dualize`/`undualize`,
  `poincare_dualize`, `left_complemented`/`left_complement`,
  `right_complemented`/`right_complement`, `norm`, `normalized`/`normalize`. The 0.1 names remain
  as deprecated methods.
- `grade` and `grade_iter` take the grade as a number, `grade(2)`, instead of `&SeaGrade(2)`.
- `iter()` yields `(BladeIndex, &T)` pairs.
- `pseudoscalar()` returns the pseudoscalar coefficient; use `grade(n)` for the part.
- Applying a motor to a flector (`Motor / Flector`) gives a `Flector`, and `Versor / Motor` a
  `Motor`: a sandwich keeps the kind of the value it is applied to.
- `inverse()` on typed values is the exact versor inverse `reverse / (v reverse(v))`, not the
  reverse; the two agree for unit versors.
- The product traits (`GeometricProduct` and the rest) take the right operand by reference.
- Edition 2024, minimum Rust 1.89 (required by nalgebra 0.35), nalgebra 0.35.

### Deprecated

- `SeaMultivector3D`..`SeaMultivector6D`, `SeaMultivectorPGA3D`, `SeaMotorPGA3D`,
  `SeaFlectorPGA3D`, `SeaVersorPGA3D`, `SeaBasisBladeIndexLow`, `ParseMultivectorError` and
  `ParseMetricError`, now aliases of the new names.
- The 0.1 method names listed under "Changed", and `basis_blade_coefficient`.

### Removed

- The traits `Multivector`, `MultivectorMut`, `MultivectorConstMetric`, `Versor`, `Motor` and
  `Flector`, the five iterator traits and the unary operation traits; their methods are inherent.
  `Multivector`, `Motor`, `Flector` and `Versor` now name types.
- The public fields `metric`, `multivector`, `basis_blade_set`, `coefficients` and the `.0` of the
  PGA3D wrappers; use the accessors.
- `SeaGrade`, `SeaMetric`, `SeaBasisBladeLow` and `IntoVersorError`.
- `Motor += Flector` and `Motor += Versor`, which always panicked; mixed-kind sums are
  multivectors.
- The dependencies `regex`, `lazy_static` and `funty`.

### Fixed

- `exp` of arguments whose square is not a scalar: it squared once too often above a norm of 10
  (`exp(10 + 10e1)` returned `exp(2A)`), and its 14-term series was up to 1.7% off below that. It
  now scales by the sum of absolute coefficients until that sum is at most 1/2, which bounds every
  power of the argument, and is exact to `f64` precision.
- `exp` of arguments whose square is zero, such as the ideal lines behind translation motors,
  takes the closed form `1 + x` instead of the series.
- The 6D left complement had four wrong signs, so `undual(dual(x))` was not `x` in 6D algebras
  with a null vector. The complement and dual sign tables are now derived from the reordering sign
  at compile time.
- `Motor::log` lost accuracy near half turns, where it divided by a scalar part close to zero, and
  returned a scalar instead of a bivector for exact half turns. It now uses a closed form that
  never divides by the scalar part.
- `to_dual_quaternion` (previously private) passed its quaternion parts to nalgebra in the wrong
  order; it now agrees with nalgebra's own conversions.

## [0.1.0]

First version, used as a git dependency.
