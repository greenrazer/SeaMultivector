# Seagrass Math

Sparse geometric algebra for Rust. Pick any signature up to 6 dimensions, at compile time or at runtime, build any multivector, and pay only for the basis blades you use.

```rust
use seagrass_math::prelude::*;

// 6D Euclidean geometric algebra: 64 basis blades.
let a = Multivector::<f64, Euclid6>::from_str("e1 + 2e2").unwrap();
let b = Multivector::<f64, Euclid6>::from_str("3e5 + e6").unwrap();

// The product visits 2 × 2 = 4 coefficient pairs instead of 64 × 64 = 4,096,
// and the result stores only the 4 blades that can be non-zero.
let ab = &a * &b;
assert_eq!(ab.to_string(), "6,0,0>3e15 + 6e25 + e16 + 2e26");
assert_eq!(ab.iter().count(), 4);
```

* Any signature (p, q, r) from 3 to 6 dimensions, fixed at compile time or chosen at runtime, including null basis vectors for projective geometry
* Sparse values: a bitset of the blades present plus only their coefficients, stored inline up to 16 blades
* Products whose structure comes from [multiplying a whole bitset of blades at once](ARCHITECTURE.md#the-trick-multiply-a-whole-set-of-blades-at-once), worked out once per pair of operand structures and cached
* Geometric, outer, regressive, inner, contraction, and sandwich products as operators
* Motors, flectors, and versors whose grade parity is checked at compile time, in any algebra
* A PGA3D toolkit with points, planes, lines, join and meet, `exp` and `log`, and nalgebra conversions
* Batches of millions of values, run from the same plans on every core, with explicit SIMD instructions, or on a GPU through wgpu
* Generic over `f32` and `f64`

## Contents

* [Why Seagrass](#why-seagrass)
* [Getting started](#getting-started)
* [Examples](#examples)
* [API overview](#api-overview)
* [Performance](#performance)
* [How it compares](#how-it-compares)
* [Status and limitations](#status-and-limitations)
* [Running the tests](#running-the-tests)

## Why Seagrass

Geometric algebra libraries usually make you choose between two designs.

* Dense, general libraries accept any algebra and any multivector. They store all 2ⁿ coefficients of every value, and a product of two values visits all 2ⁿ × 2ⁿ coefficient pairs, most of which are zero.
* Specialized libraries generate code for one algebra and a fixed menu of object types such as points, planes, and motors. They are fast, but a new signature, a higher dimension, or a mixed-grade intermediate result means regenerating code or reaching for a different tool.

Seagrass takes a third route. Every multivector carries a bitset recording which basis blades it holds, and every product works out the structure of its result with [one bit trick](ARCHITECTURE.md#the-trick-multiply-a-whole-set-of-blades-at-once) before it touches a single float. That structure depends only on which blades the operands hold, so it is computed once per pair of operand structures, cached, and replayed as a list of multiply-adds. [ARCHITECTURE.md](ARCHITECTURE.md) explains each step.

| Property          | Dense, general      | Specialized, generated                | Seagrass                                            |
|-------------------|---------------------|---------------------------------------|-----------------------------------------------------|
| Signature         | Any                 | Fixed when the code is generated      | Any (p, q, r), at compile time or at runtime        |
| Values            | Any multivector     | Predefined object types               | Any multivector                                     |
| Storage per value | All 2ⁿ coefficients | Fixed per type                        | One integer plus the coefficients that are stored   |
| Work per product  | All 2ⁿ × 2ⁿ pairs   | Non-zero pairs, decided ahead of time | Non-zero pairs, decided once per pair of structures |

The result is a library that is several times faster than a dense general library in the larger algebras and much slower than code generated for one algebra. [Performance](#performance) has the measurements.

## Getting started

### Installation

```toml
[dependencies]
seagrass-math = "0.2"
```

The default `nalgebra` feature adds the PGA3D conversions to and from nalgebra matrices, vectors, and dual quaternions; turn it off with `default-features = false` if you do not need them. Three more features are off by default: `rayon` splits batch products across every core, `simd` runs them with explicit SIMD instructions chosen when the program runs (AVX-512, AVX2, SSE, or NEON), and `gpu` runs them on a GPU through wgpu.

Everything lives in the prelude, which also re-exports `FromStr` and the nalgebra `Matrix4` and `Vector3` types used by the PGA3D conversions:

```rust
use seagrass_math::prelude::*;
```

### Your first multivectors

A multivector's second type parameter is its metric. A runtime metric reads its signature from the text before `>`: the number of positive, negative, and null basis vectors.

```rust
use seagrass_math::prelude::*;

let mv = |s: &str| Multivector::<f32, Dyn4>::from_str(s).unwrap();

// Signature (3, 1, 0): e1, e2, and e3 square to +1, and e4 squares to -1.
let a = mv("3,1,0>e1 + e12");
let b = mv("3,1,0>e3 + e123 + e4");

assert_eq!(&a * &b, mv("3,1,0>-e3 + e13 + e14 + e23 + e123 + e124")); // geometric product
assert_eq!(&a ^ &b, mv("3,1,0>e13 + e14 + e123 + e124"));            // outer (wedge) product
assert_eq!(&a | &b, mv("3,1,0>-e3 + e23"));                           // inner product
assert_eq!(&a << &b, mv("3,1,0>-e3 + e23"));                          // left contraction
assert_eq!(&a >> &b, mv("3,1,0>"));                                   // right contraction (zero)
```

Building values from numbers skips the parser, which matters in hot loops:

```rust
use seagrass_math::prelude::*;

let sta = Sta::default(); // spacetime algebra, signature (1, 3, 0)
let event = Multivector::<f64, Sta>::from_blades(sta, [(BladeIndex::E1, 2.0), (BladeIndex::E2, 1.0)]).unwrap();
assert_eq!((&event * &event).scalar(), 3.0); // the interval t² - x²
assert_eq!(event.to_string(), "1,3,0>2e1 + e2");
```

### Choosing an algebra

A compile-time signature costs nothing at runtime, and mixing two of them is a type error. A runtime signature is 3 bytes in every value and is checked when two values meet.

| Metric type                | Signature fixed  | Examples                                                                        |
|----------------------------|------------------|---------------------------------------------------------------------------------|
| `Signature<P, Q, R>`       | At compile time  | `Euclid3`, `Pga3` (3, 0, 1), `Cga3` (4, 1, 0), `Sta` (1, 3, 0), `Signature<3, 3, 0>` |
| `Dynamic<B>`, `Dyn3`..`Dyn6` | At runtime       | `"3>"` Euclidean 3D, `"2,0,1>"` PGA 2D, `"5,1>"` CGA 4D                          |

The dimension p + q + r picks the bitset type: `u8` for 3D (8 blades), `u16` for 4D, `u32` for 5D, and `u64` for 6D. Basis vectors are numbered from 1 and ordered positive, then negative, then null, so in (3, 0, 1) e4 is the null projective vector that many PGA texts call e0. `T` is `f32` or `f64`.

### String syntax

* Terms look like `2.5e13` or `-e124`, and a bare number such as `3` is the scalar part. There is no exponent notation: `1e5` is the blade e5.
* Basis vectors inside a blade must be ascending, so `e13` is valid and `e31` is rejected.
* An empty right-hand side, as in `"3>"`, is the zero multivector.
* The `p,q,r>` prefix is required for runtime metrics and optional for compile-time ones, where it must match.
* `Display` writes the same syntax, so `value.to_string().parse()` gives back an equal value.
* Parsing returns a `Result` whose error names the problem, such as `MissingMetric`, `InvalidBasisVectorOrdering`, `DuplicateBasisBlade`, or `UnexpectedToken { position }`.

## Examples

### One engine, many algebras

The same types and operators cover Euclidean, spacetime, and conformal algebras. `exp()` switches to a closed form whenever its argument squares to a scalar.

```rust
use std::f64::consts::FRAC_PI_2;
use seagrass_math::prelude::*;

// Euclidean 3D: rotate e1 by 90° in the e12 plane with the rotor R = exp(-θ/2 · e12).
let vga = |s: &str| Multivector::<f64, Euclid3>::from_str(s).unwrap();
let rotor = (vga("e12") * (-FRAC_PI_2 / 2.0)).exp();
let mut v = rotor.sandwich(&vga("e1")); // R e1 R̃
v.simplify();
assert_eq!(v, vga("e2"));

// Spacetime algebra, signature (1, 3): a boost is the exponential of a timelike bivector.
let sta = |s: &str| Multivector::<f64, Sta>::from_str(s).unwrap();
let rapidity = 2f64.ln(); // velocity 0.6c, Lorentz factor 1.25
let boost = (sta("e12") * (rapidity / 2.0)).exp(); // e12² = +1, so this is cosh + sinh·e12
let mut t = boost.sandwich(&sta("e1"));
t.round_to_decimals(12);
assert_eq!(t, sta("1.25e1 - 0.75e2"));

// Conformal 3D, signature (4, 1): the inner product of two points is -½·distance².
let cga = |s: &str| Multivector::<f64, Cga3>::from_str(s).unwrap();
let (n_inf, n_0) = (cga("e4 + e5"), cga("-0.5e4 + 0.5e5"));
let point = |x: f64, y: f64, z: f64| {
    cga("e1") * x + cga("e2") * y + cga("e3") * z + &n_inf * (0.5 * (x * x + y * y + z * z)) + &n_0
};
let distance_squared = -2.0 * (&point(1.0, 2.0, 3.0) | &point(4.0, 6.0, 3.0)).scalar();
assert_eq!(distance_squared, 25.0);
```

[`examples/four_algebras.rs`](examples/four_algebras.rs) runs one generic function in four algebras chosen by metric strings: it rotates a vector, translates a point, boosts an event, and dilates a conformal point.

### Projective geometry for graphics

In PGA3D, points and planes are flectors (odd), and lines and rigid motions are motors (even). `&` joins and `^` meets, and the kind of every result is known at compile time.

```rust
use std::f64::consts::{FRAC_1_SQRT_2, FRAC_PI_2};
use seagrass_math::prelude::*;

let point = Multivector::<f64, Pga3>::new_point;
let plane = Multivector::<f64, Pga3>::new_plane;

// Join three points into the plane through them: -7x - 3y + 4z + 1 = 0.
let (a, b, c) = (point(1.0, 2.0, 3.0), point(2.0, 1.0, 4.0), point(2.0, 5.0, 7.0));
let through: Flector<f64, Pga3> = &a & &b & &c;
assert_eq!(through, plane(-7.0, -3.0, 4.0, 1.0));

// Rigid motions are motors. `*` composes them right to left, like matrices,
// and `/` applies one with the sandwich product M x M̃.
let rotate = Motor::<f64, Pga3>::new_rotation_motor(0.0, 0.0, FRAC_PI_2); // 90° about the z axis
let translate = Motor::<f64, Pga3>::new_translation_motor(0.0, 0.0, 5.0);
let motor = &translate * &rotate;

let p = &motor / point(1.0, 0.0, 0.0);
assert!((p.point_to_euclidean() - Vector3::new(0.0, 1.0, 5.0)).norm() < 1e-12);

// Screw interpolation: exp(½·log M) moves exactly halfway along the same motion.
let halfway = (motor.log() * 0.5).exp();
let q = &halfway / point(1.0, 0.0, 0.0);
assert!((q.point_to_euclidean() - Vector3::new(FRAC_1_SQRT_2, FRAC_1_SQRT_2, 2.5)).norm() < 1e-12);
```

### A pinhole camera in three operations

Projection is a motor, a join, and a meet, with no projection matrix involved.

```rust
use seagrass_math::prelude::*;

let eye = Multivector::<f64, Pga3>::new_point(0.0, 0.0, -1.0);
let image_plane = Multivector::<f64, Pga3>::new_plane(0.0, 0.0, 1.0, 0.0); // z = 0
let view = Motor::<f64, Pga3>::new_translation_motor(-3.0, -3.0, 0.0);      // world to camera space

let point = &view / Multivector::<f64, Pga3>::new_point(3.0, 6.0, 2.0); // 1. move into camera space
let ray = &eye & &point;                                                // 2. join: the line through eye and point
let pixel = ray ^ &image_plane;                                         // 3. meet: where that line hits the image

assert!((pixel.point_to_euclidean() - Vector3::new(0.0, 1.0, 0.0)).norm() < 1e-12);
```

### Interop with nalgebra

Motors convert to and from homogeneous matrices, and the inverse matrix falls out of the same closed-form expression ([how](ARCHITECTURE.md#sandwiches-and-matrices)), so no matrix inversion runs. Motors also convert from quaternions and to nalgebra dual quaternions.

```rust
use seagrass_math::prelude::*;

let m = Matrix4::new_translation(&Vector3::new(1.0, 2.0, 3.0)) * Matrix4::from_euler_angles(0.3, -0.2, 1.1);
let motor = Motor::<f64, Pga3>::from(m);

let (forward, inverse) = motor.to_homogeneous_transformation_matrices();
assert!((forward - m).norm() < 1e-12);
assert!((forward * inverse - Matrix4::identity()).norm() < 1e-12);
```

### Millions of values at once

A `Batch` holds many values that share one set of blades, stored as one array per blade, and runs each product's plan once for the whole batch. A `Backend` decides where that happens: `CpuBackend` on one thread or, with the `rayon` feature, on every core; `SimdBackend` (the `simd` feature) the same way with explicit SIMD instructions; and `GpuBackend` (the `gpu` feature) on a GPU in `f32`. The same code runs on any of them.

```rust
use seagrass_math::prelude::*;

// Moves every point by its own motor, on any backend.
fn move_points<B: Backend<f32>>(backend: &B, motors: &Batch<f32, Pga3>, points: &Batch<f32, Pga3>) -> Batch<f32, Pga3> {
    let (motors, points) = (backend.upload(motors), backend.upload(points));
    backend.download(&backend.sandwich(&motors, &points))
}

let motors: Vec<Motor<f32, Pga3>> = (0..100_000).map(|i| Motor::new_translation_motor(i as f32, 0.0, 0.0)).collect();
let motors = Batch::from_values(Pga3::default(), &motors);
let points = Batch::<f32, Pga3>::from_points(&vec![Vector3::new(0.0, 1.0, 0.0); 100_000]);

let moved = move_points(&CpuBackend::parallel(), &motors, &points);
assert_eq!(moved.points_to_euclidean()[3], Vector3::new(3.0, 1.0, 0.0));
// With the `simd` feature: move_points(&SimdBackend::parallel(), &motors, &points)
// With the `gpu` feature: move_points(&GpuBackend::new()?, &motors, &points)
```

On the GPU, products queue into one submission and their results stay in device memory until downloaded, so a chain of products never crosses the bus. Batches below 20,000 elements stay on the CPU, where they are faster ([benchmarks](BENCHMARKS.md#batches-on-the-gpu-and-the-cpu)). The GPU backend tries Vulkan first, then Metal, DX12, and OpenGL; set `WGPU_BACKEND` to choose.

## API overview

### Types

| Type                      | Represents                                                                 |
|---------------------------|----------------------------------------------------------------------------|
| `Multivector<T, M>`       | Any multivector of the algebra with metric `M`                             |
| `Motor<T, M>`             | An even multivector: a rotor, motor, or line (`Typed<T, M, Even>`)         |
| `Flector<T, M>`           | An odd multivector: a reflection, point, or plane (`Typed<T, M, Odd>`)     |
| `Versor<T, M>`            | A multivector whose parity is not tracked (`Typed<T, M, Any>`)             |
| `BladeIndex`              | A basis blade, with constants `BladeIndex::E1` to `BladeIndex::E123456`    |
| `Signature<P, Q, R>`, `Dynamic<B>` | Compile-time and runtime metrics                                  |

Products of typed values are typed: a motor times a motor is a motor, a motor times a flector is a flector, and a motor applied to a flector with `/` is a flector. Mixing a typed value with a plain multivector gives a plain multivector, and `Typed::new` (or `try_from`) checks the parity of a plain one.

### Operators

| Expression                         | Meaning                                                                  |
|------------------------------------|--------------------------------------------------------------------------|
| `a * b`                            | Geometric product                                                        |
| `a ^ b`                            | Outer (wedge) product, the meet in PGA                                   |
| `a & b`                            | Regressive (vee) product, the join in PGA                                |
| `a \| b`                           | Inner product: left plus right contraction, scalar part counted once     |
| `a << b`                           | Left contraction                                                         |
| `a >> b`                           | Right contraction                                                        |
| `!a`                               | Dual                                                                     |
| `-a`                               | Negation                                                                 |
| `a + b`, `a - b`                   | Sum and difference                                                       |
| `a + s`, `a - s`, `a * s`, `a / s` | Arithmetic with a scalar `s: T`                                          |
| `m / x`                            | Sandwich product: applies the motor, flector, or versor `m` to `x`       |

Binary operators accept owned values or references, so `&a * &b` multiplies without moving either operand. Each operator forwards to a named trait method, such as `GeometricProduct::geometric_product`, for generic code. Operands with different runtime metrics panic.

### Methods

Most unary operations come in two forms: one returns a new value and the other works in place.

| Returns a new value                       | In place                              | Operation                   |
|-------------------------------------------|---------------------------------------|-----------------------------|
| `reversed()`                              | `reverse()`                           | Reverse                     |
| `involuted()`                             | `involute()`                          | Grade involution            |
| `conjugated()`                            | `conjugate()`                         | Clifford conjugate          |
| `negated()`                               | `negate()`                            | Negation                    |
| `dual()`, `undual()`                      | `dualize()`, `undualize()`            | Metric dual and its inverse |
| `poincare_dual()`                         | `poincare_dualize()`                  | Metric-free dual            |
| `left_complemented()`, `right_complemented()` | `left_complement()`, `right_complement()` | Left and right complements  |
| `normalized()`                            | `normalize()`                         | Divide by `norm()`          |
| `exp()`                                   |                                       | Exponential                 |

Other helpers include `get(index)`, `scalar()`, and `pseudoscalar()` for coefficients, `grade(k)` for a grade part, `iter()` and `iter_mut()` over the stored blades, `sandwich(x)`, `inverse()` and `apply(x)` on typed values, and cleanup with `remove_zeros()`, `remove_approx_zeros(eps)`, `round_to_decimals(n)`, and `simplify()`. Values are built with `from_str`, `zero`, `from_scalar`, `basis`, `from_blades`, and `from_parts`.

### PGA3D

| Type                       | Highlights                                                                                                                                       |
|----------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------|
| `Multivector<T, Pga3>`     | `new_point`, `new_homogeneous_point`, `new_plane`, `new_line`, `new_ideal_line`, `new_motor`, `point_to_euclidean`, `point_to_homogeneous`, `raw` |
| `Motor<T, Pga3>`           | `new_rotation_motor`, `new_translation_motor`, `new_combined_motor`, `from_quaternion`, `log`, `From<Matrix4>`, `to_homogeneous_transformation_matrices`, `to_dual_quaternion` |
| `Flector<T, Pga3>`         | `from_point`, `from_plane`, `from_origin_plane`, `From<Matrix4>`, `to_homogeneous_transformation_matrices`                                         |

`new_point` and `new_plane` return flectors, while `new_line`, `new_ideal_line`, and `new_motor` return motors.

### Names from 0.1

The 0.1 type names (`SeaMultivector3D` to `SeaMultivector6D`, `SeaMultivectorPGA3D`, `SeaMotorPGA3D`, `SeaFlectorPGA3D`, `SeaVersorPGA3D`) and method names (`reversion_conjugated()` and the rest) still compile as deprecated aliases. [`CHANGELOG.md`](CHANGELOG.md) lists every change.

## Performance

Each row compares Seagrass with the fastest other library doing the same work, on an AMD Ryzen 9 5900X and an NVIDIA GeForce RTX 3090:

| Task                                                                    | Fastest other library                         | Seagrass is        |
|-------------------------------------------------------------------------|-----------------------------------------------|--------------------|
| One PGA3D product: motor × motor, or a motor applied to a point         | geometric_algebra (`f32`), nalgebra `Matrix4` | 18 to 22× slower   |
| One CGA (4, 1) product: vector × vector, or a rotor applied to a vector | clifford (dense)                              | 2.9 to 6.9× faster |
| Move 100k points on one CPU thread (`SimdBackend`)                      | glam `Mat4 × Vec4`                            | 3.2× slower        |
| Move 4M points on the GPU (`GpuBackend`)                                | a hand-written WGSL shader                    | 1.5× faster        |

* Libraries built for one 3D algebra win on the CPU, because they know each product's structure at compile time. Seagrass still beats some of them: on one thread it combines motions 1.4 to 2.5 times as fast as glam and nalgebra, and moves points 8 times as fast as geometric_algebra.
* The GPU wins only while data stays on it. Copying each batch there and back is slower than one CPU thread.

[`BENCHMARKS.md`](BENCHMARKS.md) compares Seagrass with each library at several batch sizes and thread counts, with every number and how it was measured.

## How it compares

Here is how Seagrass sits next to some well-known geometric algebra libraries, based on their documentation as of October 2026:

| Library                                                           | Language   | Algebras                                          | Algebra chosen                  | Values                                            | SIMD         |
|-------------------------------------------------------------------|------------|---------------------------------------------------|---------------------------------|---------------------------------------------------|--------------|
| **Seagrass**                                                      | Rust       | Any (p, q, r) in 3 to 6 dimensions                | At compile time or at runtime   | Any multivector, stored sparsely                  | Batches only |
| [geometric_algebra](https://github.com/Lichtso/geometric_algebra) | Rust       | Prebuilt PGAs up to 3D, others from its generator | By code generation              | Generated types and a general multivector         | Yes          |
| [clifford](https://github.com/DevonMorris/clifford)               | Rust       | 15 prebuilt algebras plus signature types         | At compile time                 | Specialized types and a dense generic multivector | No           |
| [g3](https://github.com/wrnrlr/g3)                                | Rust       | PGA 3D                                            | Built in                        | Predefined object types                           | Yes, nightly |
| [ultraviolet](https://github.com/termhn/ultraviolet)              | Rust       | Euclidean 2D and 3D rotors                        | Built in                        | Rotors and bivectors only                         | Yes          |
| [klein](https://github.com/jeremyong/klein)                       | C++        | PGA 3D                                            | Built in                        | Predefined object types                           | Yes, SSE     |
| [ganja.js](https://github.com/enkimute/ganja.js)                  | JavaScript | Any (p, q, r)                                     | At runtime, by generating code  | Any multivector, dense up to 6D                   | No           |
| [kingdon](https://github.com/tBuLi/kingdon)                       | Python     | Any (p, q, r)                                     | At runtime, by symbolic codegen | Any multivector, stored sparsely                  | No           |
| [clifford](https://github.com/pygae/clifford)                     | Python     | Any (p, q, r)                                     | At runtime                      | Any multivector, stored densely                   | No           |

Two things set Seagrass apart in that list.

* It is the only Rust library here that can pick the signature at runtime and stores arbitrary multivectors sparsely. The general-purpose Rust crates fix the algebra at compile time or through a generator, and the SIMD-focused ones offer a fixed menu of types.
* Its closest relative is kingdon, which also exploits sparsity at runtime but generates and compiles code for each new combination of input blades. Seagrass gets the same structural information from bit operations and caches it as a plan, so there is no code-generation step on first use, and it runs as native Rust.

If you need more than 6 dimensions, ganja.js, kingdon, or the Python clifford package handle larger algebras.

## Status and limitations

* Seagrass is pre-1.0, so expect API changes; deprecated names stay for at least one release.
* Dimensions run from 3 to 6, one bitset width per dimension.
* Plans are cached per thread. A thread that meets more than 4,096 distinct structure pairs clears its cache and rebuilds plans as needed, which costs a few microseconds each.
* Values with more than 16 stored blades (only possible in 5D and 6D) keep their coefficients on the heap.
* Operands with different runtime metrics panic; compile-time signatures make that a type error instead.
* The SIMD backend has kernels for `f32` and `f64`; other coefficient types run the `CpuBackend` loops. Its tests ran on SSE2, SSE4.2, and AVX2; the AVX-512 and NEON paths are compiled from the same code by fearless_simd but have not run on such hardware.
* The GPU backend computes in `f32` only, and each batch shares one set of blades across its values. It compiles a shader the first time it meets a pair of blade sets, which blocked for 1 to 10 ms on the benchmark machine; `ShaderMode::Interpreted` avoids that at a fifth to a sixth of the speed.
* On the Windows machine behind these benchmarks, wgpu 30 could not create a DX12 device, which is why the backend tries Vulkan first.

## Running the tests

```bash
cargo test
cargo bench --bench products   # single-value timings
cargo bench --bench compare    # single values against other crates
cargo test --release --features gpu,rayon,simd                              # batches on every core, with SIMD, and on the GPU
cargo run --release --example gpu_bench --features gpu,rayon,simd           # the batch tables in BENCHMARKS.md
cargo run --release --example gpu_bench --features gpu,rayon,simd -- --cpu  # the CPU and SIMD batch tables, without a GPU
```

Golden values cover products under Euclidean, projective, hyperbolic, and mixed signatures such as (3, 1, 1) in up to 6 dimensions, duals with two null vectors as in (3, 0, 2), parsing, `exp` and `log`, nalgebra conversions, and PGA camera projections. Property tests check associativity, reversion, duals, complements, quadratic forms, and the `Display` round trip on random values in twelve signatures. In all 74 signatures from 3 to 6 dimensions, at compile time and at runtime, every product, dual, complement, and reversion is checked twice: against a reference in the tests that multiplies blades by hand, and against golden values computed once by an independent implementation and stored in `tests/fixtures/golden`, one file per signature. Every cached plan is also checked against the per-call product it replaces.
