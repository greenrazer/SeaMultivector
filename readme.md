# Seagrass Math

Sparse geometric algebra for Rust. Pick any metric at runtime, build any multivector, and pay only for the basis blades you use.

```rust
use seagrass_math::prelude::*;

// 6D Euclidean geometric algebra: 64 basis blades.
let a = SeaMultivector6D::<f64>::from_str("6>e1 + 2e2").unwrap();
let b = SeaMultivector6D::<f64>::from_str("6>3e5 + e6").unwrap();

// The product visits 2 × 2 = 4 coefficient pairs instead of 64 × 64 = 4,096,
// and the result stores only the 4 blades that can be non-zero.
let ab = &a * &b;
assert_eq!(ab, SeaMultivector6D::<f64>::from_str("6>3e15 + e16 + 6e25 + 2e26").unwrap());
assert_eq!(ab.iter().count(), 4);
```

* Any signature (p, q, r) from 3 to 6 dimensions, including null basis vectors for projective geometry
* Sparse values: a bitset of the blades present plus only their coefficients
* Products that skip structural zeros, found at runtime with bitwise operations
* Geometric, outer, regressive, inner, and contraction products as operators
* A PGA3D toolkit with points, planes, lines, motors, join and meet, `exp` and `log`, and nalgebra conversions
* Generic over `f32` and `f64`

## Contents

* [Why Seagrass](#why-seagrass)
* [Getting started](#getting-started)
* [Examples](#examples)
* [How it works](#how-it-works)
* [API overview](#api-overview)
* [How it compares](#how-it-compares)
* [Status and limitations](#status-and-limitations)
* [Running the tests](#running-the-tests)

## Why Seagrass

Geometric algebra libraries usually make you choose between two designs.

* Dense, general libraries accept any algebra and any multivector. They store all 2ⁿ coefficients of every value, and a product of two values visits all 2ⁿ × 2ⁿ coefficient pairs, most of which are zero.
* Specialized libraries generate code for one algebra and a fixed menu of object types such as points, planes, and motors. They are fast, but a new signature, a higher dimension, or a mixed-grade intermediate result means regenerating code or reaching for a different tool.

Seagrass takes a third route. Every multivector carries a bitset recording which basis blades it holds, and every product works out the structure of its result with bitwise operations before it touches a single float. You keep the generality of the first design and do work in line with the second.

| Property          | Dense, general      | Specialized, generated                | Seagrass                                          |
|-------------------|---------------------|---------------------------------------|---------------------------------------------------|
| Signature         | Any                 | Fixed when the code is generated      | Any (p, q, r), chosen at runtime                  |
| Values            | Any multivector     | Predefined object types               | Any multivector                                   |
| Storage per value | All 2ⁿ coefficients | Fixed per type                        | One integer plus the coefficients that are stored |
| Work per product  | All 2ⁿ × 2ⁿ pairs   | Non-zero pairs, decided ahead of time | Non-zero pairs, decided at runtime                |

### The work you skip

The table counts floating-point multiply-adds for a single product. The dense column multiplies every coefficient pair of two full 2ⁿ-element arrays. The Seagrass column is the number of stored blade pairs that can contribute to the result, which is exactly what its inner loop visits. Pairs that share a null basis vector vanish, and Seagrass skips those too, which matters in PGA.

| Product                             | Algebra (basis blades) | Dense | Seagrass |
|-------------------------------------|------------------------|------:|---------:|
| Rotor × vector                      | Euclidean 3D (8)       |    64 |       12 |
| Motor × motor                       | PGA 3D (16)            |   256 |       48 |
| Motor applied to a point (sandwich) | PGA 3D (16)            |   512 |       68 |
| Vector × vector                     | CGA 3D (32)            | 1,024 |       25 |
| Vector × vector                     | Euclidean 6D (64)      | 4,096 |       36 |
| Bivector × vector                   | Euclidean 6D (64)      | 4,096 |       90 |

> [!NOTE]
> A dense implementation that prunes pairs which are always zero in the PGA product table would drop the two PGA rows to 192 and 384. Specialized libraries skip zeros too, and often go further with hand-derived formulas, but only for the type combinations they ship. Seagrass finds the non-zero pairs at runtime for any two multivectors.

## Getting started

### Installation

Seagrass is not on crates.io yet. Add it as a git dependency:

```toml
[dependencies]
seagrass-math = { git = "https://github.com/greenrazer/SeaMultivector" }
```

Everything lives in the prelude, which also re-exports `FromStr` and the nalgebra `Matrix4` and `Vector3` types used by the PGA conversions:

```rust
use seagrass_math::prelude::*;
```

### Your first multivectors

A general multivector carries its metric. The text before `>` picks the signature: the number of positive, negative, and null basis vectors.

```rust
use seagrass_math::prelude::*;

let mv = |s: &str| SeaMultivector4D::<f32>::from_str(s).unwrap();

// Signature (3, 1, 0): e1, e2, and e3 square to +1, and e4 squares to -1.
let a = mv("3,1,0>e1 + e12");
let b = mv("3,1,0>e3 + e123 + e4");

assert_eq!(&a * &b, mv("3,1,0>-e3 + e13 + e14 + e23 + e123 + e124")); // geometric product
assert_eq!(&a ^ &b, mv("3,1,0>e13 + e14 + e123 + e124"));            // outer (wedge) product
assert_eq!(&a | &b, mv("3,1,0>-e3 + e23"));                           // inner product
assert_eq!(&a << &b, mv("3,1,0>-e3 + e23"));                          // left contraction
assert_eq!(&a >> &b, mv("3,1,0>"));                                   // right contraction (zero)
```

### Choosing an algebra

Pick a type by dimension, then any signature whose counts add up to that dimension.

| Type                  | Dimension | Basis blades | Bitset | Example signatures                 |
|-----------------------|-----------|--------------|--------|------------------------------------|
| `SeaMultivector3D<T>` | 3         | 8            | `u8`   | `3>` Euclidean 3D, `2,0,1>` PGA 2D |
| `SeaMultivector4D<T>` | 4         | 16           | `u16`  | `3,0,1>` PGA 3D, `1,3>` spacetime  |
| `SeaMultivector5D<T>` | 5         | 32           | `u32`  | `4,1>` CGA 3D, `4,0,1>` PGA 4D     |
| `SeaMultivector6D<T>` | 6         | 64           | `u64`  | `6>` Euclidean 6D, `5,1>` CGA 4D   |

Basis vectors are numbered from 1 and ordered positive, then negative, then null. In `3,0,1>`, e4 is the null projective vector that many PGA texts call e0. `T` is `f32` or `f64`.

### String syntax

* Terms look like `2.5e13` or `-e124`, and a bare number such as `3` is the scalar part.
* Basis vectors inside a blade must be ascending, so `e13` is valid and `e31` is rejected.
* An empty right-hand side, as in `"3>"`, is the zero multivector.
* Parsing returns a `Result` whose error names the problem, such as `MissingMetric`, `InvalidBasisVectorOrdering`, `DuplicateBasisBlade`, or `InvalidMetricDimensions`.

## Examples

### One engine, many algebras

The same types and operators cover Euclidean, spacetime, and conformal algebras. `exp()` switches to a closed form whenever its argument squares to a scalar.

```rust
use std::f64::consts::FRAC_PI_2;
use seagrass_math::prelude::*;

// Euclidean 3D: rotate e1 by 90° in the e12 plane with the rotor R = exp(-θ/2 · e12).
let vga = |s: &str| SeaMultivector3D::<f64>::from_str(&format!("3>{s}")).unwrap();
let rotor = (vga("e12") * (-FRAC_PI_2 / 2.0)).exp();
let mut v = &(&rotor * &vga("e1")) * &rotor.reversion_conjugated();
v.remove_approx_zeros(1e-12);
v.round_to_decimals(12);
assert_eq!(v, vga("e2"));

// Spacetime algebra, signature (1, 3): a boost is the exponential of a timelike bivector.
let sta = |s: &str| SeaMultivector4D::<f64>::from_str(&format!("1,3>{s}")).unwrap();
let rapidity = 2f64.ln(); // velocity 0.6c, Lorentz factor 1.25
let boost = (sta("e12") * (rapidity / 2.0)).exp(); // e12² = +1, so this is cosh + sinh·e12
let mut t = &(&boost * &sta("e1")) * &boost.reversion_conjugated();
t.round_to_decimals(12);
assert_eq!(t, sta("1.25e1 - 0.75e2"));

// Conformal 3D, signature (4, 1): the inner product of two points is -½·distance².
let cga = |s: &str| SeaMultivector5D::<f64>::from_str(&format!("4,1>{s}")).unwrap();
let (n_inf, n_0) = (cga("e4 + e5"), cga("-0.5e4 + 0.5e5"));
let point = |x: f64, y: f64, z: f64| {
    cga("e1") * x + cga("e2") * y + cga("e3") * z + &n_inf * (0.5 * (x * x + y * y + z * z)) + &n_0
};
let distance_squared = -2.0 * (&point(1.0, 2.0, 3.0) | &point(4.0, 6.0, 3.0)).scalar();
assert_eq!(distance_squared, 25.0);
```

### Projective geometry for graphics

The PGA3D types fix the (3, 0, 1) metric at compile time, so nothing is parsed or checked at runtime. Points, planes, and lines have constructors, `&` joins, and `^` meets.

```rust
use std::f64::consts::{FRAC_1_SQRT_2, FRAC_PI_2};
use seagrass_math::prelude::*;

// Join three points into the plane through them: -7x - 3y + 4z + 1 = 0.
let a = SeaMultivectorPGA3D::new_point(1.0, 2.0, 3.0);
let b = SeaMultivectorPGA3D::new_point(2.0, 1.0, 4.0);
let c = SeaMultivectorPGA3D::new_point(2.0, 5.0, 7.0);
assert_eq!(&a & &b & &c, SeaMultivectorPGA3D::new_plane(-7.0, -3.0, 4.0, 1.0));

// Rigid motions are motors. `*` composes them right to left, like matrices,
// and `/` applies one with the sandwich product M x M̃.
let rotate = SeaMotorPGA3D::new_rotation_motor(0.0, 0.0, FRAC_PI_2); // 90° about the z axis
let translate = SeaMotorPGA3D::new_translation_motor(0.0, 0.0, 5.0);
let motor = &translate * &rotate;

let p = &motor / SeaMultivectorPGA3D::new_point(1.0, 0.0, 0.0);
assert!((p.point_to_euclidean() - Vector3::new(0.0, 1.0, 5.0)).norm() < 1e-12);

// Screw interpolation: exp(½·log M) moves exactly halfway along the same motion.
let halfway = (motor.log() * 0.5).exp();
let q = &halfway / SeaMultivectorPGA3D::new_point(1.0, 0.0, 0.0);
assert!((q.point_to_euclidean() - Vector3::new(FRAC_1_SQRT_2, FRAC_1_SQRT_2, 2.5)).norm() < 1e-12);
```

### A pinhole camera in three operations

Projection is a motor, a join, and a meet, with no projection matrix involved.

```rust
use seagrass_math::prelude::*;

let eye = SeaMultivectorPGA3D::new_point(0.0, 0.0, -1.0);
let image_plane = SeaMultivectorPGA3D::new_plane(0.0, 0.0, 1.0, 0.0); // z = 0
let view = SeaMotorPGA3D::new_translation_motor(-3.0, -3.0, 0.0);    // world to camera space

let point = &view / SeaMultivectorPGA3D::new_point(3.0, 6.0, 2.0); // 1. move into camera space
let ray = &eye & &point;                                           // 2. join: the line through eye and point
let pixel = ray ^ &image_plane;                                    // 3. meet: where that line hits the image

assert!((pixel.point_to_euclidean() - Vector3::new(0.0, 1.0, 0.0)).norm() < 1e-12);
```

### Interop with nalgebra

Motors convert to and from homogeneous matrices, and the inverse matrix falls out of the same closed-form expression, so no matrix inversion runs. Motors can also be built from quaternions and dual quaternions with `from_quaternion`, `from_dual_quaternion`, and `from_combined_quaternion`.

```rust
use seagrass_math::prelude::*;

let m = Matrix4::new_translation(&Vector3::new(1.0, 2.0, 3.0)) * Matrix4::from_euler_angles(0.3, -0.2, 1.1);
let motor = SeaMotorPGA3D::<f64>::from(m);

let (forward, inverse) = motor.to_homogeneous_transformation_matrices();
assert!((forward - m).norm() < 1e-12);
assert!((forward * inverse - Matrix4::identity()).norm() < 1e-12);
```

## How it works

### Blades are bit patterns

Each basis blade is indexed by a bitmask of its basis vectors, so the blade produced by multiplying two blades is the XOR of their indices.

```text
blade   1     e1    e2    e12   e3    e13   e23   e123
index   0     1     2     3     4     5     6     7
bits    000   001   010   011   100   101   110   111
```

### A multivector is one integer and a short vector

A second bitset with one bit per blade (8 bits in 3D, 64 in 6D) records which blades are present. Coefficients are packed in index order, and a blade's slot is the number of set bits below it, which is one `count_ones` call.

```text
3 + 2e1 - e23                   3D, so the bitset is a u8

bitset         0b0100_0011      blades 0 (1), 1 (e1), and 6 (e23)
coefficients   [3.0, 2.0, -1.0]
```

### Products run in two phases

1. Structure comes first. The output bitset is computed from the input bitsets with shifts and masks, because multiplying by a basis vector permutes an entire bitset at once. The output coefficient vector is then allocated once, at its final length.
2. Numbers come second. For each output blade, a few shifts and one AND find every stored left-hand blade that has a stored partner. The signs of all those pairs, from blade reordering and from the metric, come out of a single mask assembled from precomputed 64-bit constants, so the inner loop is a plain multiply-add.

### Unary operations are word-sized

* Reversion, grade involution, and Clifford conjugation flip signs under one constant mask.
* The dual is `reverse_bits()` on the bitset, a reversed coefficient vector, and a sign mask. There is no multiplication by the pseudoscalar, which is why the dual also works in degenerate metrics such as PGA, where the pseudoscalar has no inverse.
* `exp()` uses a closed form (cosine and sine, hyperbolic cosine and sine, or 1 + x) when its argument squares to a scalar, and otherwise a Taylor series, scaled and squared for large arguments. `log()` on PGA3D motors is closed form.

## API overview

### Operators

| Expression                         | Meaning                                                                  |
|------------------------------------|--------------------------------------------------------------------------|
| `a * b`                            | Geometric product                                                        |
| `a ^ b`                            | Outer (wedge) product, the meet in PGA                                   |
| `a & b`                            | Regressive (vee) product, the join in PGA                                |
| `a \| b`                           | Inner product, the sum of the left and right contractions                |
| `a << b`                           | Left contraction                                                         |
| `a >> b`                           | Right contraction                                                        |
| `!a`                               | Dual                                                                     |
| `-a`                               | Negation                                                                 |
| `a + b`, `a - b`                   | Sum and difference                                                       |
| `a + s`, `a - s`, `a * s`, `a / s` | Arithmetic with a scalar `s: T`                                          |
| `m / x`                            | Sandwich product: applies the PGA3D motor, flector, or versor `m` to `x` |

Binary operators accept owned values or references, so `&a * &b` multiplies without moving either operand. Both operands must share a metric.

### Methods

Most unary operations come in two forms: one returns a new value and the other works in place.

| Returns a new value                           | In place                                  | Operation                   |
|-----------------------------------------------|-------------------------------------------|-----------------------------|
| `reversion_conjugated()`                      | `reversion_conjugation()`                 | Reverse                     |
| `involution_conjugated()`                     | `involution_conjugation()`                | Grade involution            |
| `clifford_conjugated()`                       | `clifford_conjugation()`                  | Clifford conjugate          |
| `negated()`                                   | `negate()`                                | Negation                    |
| `dual()`, `undual()`                          | `dualization()`, `undualization()`        | Metric dual and its inverse |
| `poincare_dual()`                             | `poincare_dualization()`                  | Metric-free dual            |
| `left_complimented()`, `right_complimented()` | `left_compliment()`, `right_compliment()` | Left and right complements  |
| `scalar_normalized()`                         | `scalar_normalize()`                      | Divide by `scalar_norm()`   |
| `exp()`                                       |                                           | Exponential                 |

Other helpers include `scalar()` and `pseudoscalar()` for those parts of a value, `iter()` and `iter_mut()` over the stored blades, and cleanup with `remove_zeros()`, `remove_approx_zeros(eps)`, and `round_to_decimals(n)` on general types or `simplify()` and `simplify_epsilon(eps)` on PGA3D types.

### PGA3D types

| Type                     | Represents                                 | Highlights                                                                                                                                                 |
|--------------------------|--------------------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `SeaMultivectorPGA3D<T>` | Any PGA3D multivector                      | `new_point`, `new_homogeneous_point`, `new_plane`, `new_line`, `new_ideal_line`, `new_motor`, `point_to_euclidean`, `point_to_homogeneous`                 |
| `SeaMotorPGA3D<T>`       | Rotations, translations, and screw motions | `new_rotation_motor`, `new_translation_motor`, `from_quaternion`, `from_dual_quaternion`, `log`, `From<Matrix4>`, `to_homogeneous_transformation_matrices` |
| `SeaFlectorPGA3D<T>`     | Points, planes, and reflections            | `from_point`, `from_plane`, `from_origin_plane`, `From<Matrix4>`, `to_homogeneous_transformation_matrices`                                                 |
| `SeaVersorPGA3D<T>`      | General versors, such as motor × flector   | The shared operators and methods                                                                                                                           |

`new_point` and `new_plane` return flectors, while `new_line`, `new_ideal_line`, and `new_motor` return motors. The traits `Multivector<T>`, `Versor<T>`, `Motor<T>`, and `Flector<T>` are exported for code that should work across types.

## How it compares

Here is how Seagrass sits next to some well-known geometric algebra libraries, based on their documentation as of October 2026:

| Library                                                           | Language   | Algebras                                          | Algebra chosen                  | Values                                            | SIMD         |
|-------------------------------------------------------------------|------------|---------------------------------------------------|---------------------------------|---------------------------------------------------|--------------|
| **Seagrass**                                                      | Rust       | Any (p, q, r) in 3 to 6 dimensions                | At runtime                      | Any multivector, stored sparsely                  | No           |
| [geometric_algebra](https://github.com/Lichtso/geometric_algebra) | Rust       | Prebuilt PGAs up to 3D, others from its generator | By code generation              | Generated types and a general multivector         | Yes          |
| [clifford](https://github.com/DevonMorris/clifford)               | Rust       | 15 prebuilt algebras plus signature types         | At compile time                 | Specialized types and a dense generic multivector | No           |
| [g3](https://github.com/wrnrlr/g3)                                | Rust       | PGA 3D                                            | Built in                        | Predefined object types                           | Yes, nightly |
| [ultraviolet](https://github.com/termhn/ultraviolet)              | Rust       | Euclidean 2D and 3D rotors                        | Built in                        | Rotors and bivectors only                         | Yes          |
| [klein](https://github.com/jeremyong/klein)                       | C++        | PGA 3D                                            | Built in                        | Predefined object types                           | Yes, SSE     |
| [ganja.js](https://github.com/enkimute/ganja.js)                  | JavaScript | Any (p, q, r)                                     | At runtime, by generating code  | Any multivector, dense up to 6D                   | No           |
| [kingdon](https://github.com/tBuLi/kingdon)                       | Python     | Any (p, q, r)                                     | At runtime, by symbolic codegen | Any multivector, stored sparsely                  | No           |
| [clifford](https://github.com/pygae/clifford)                     | Python     | Any (p, q, r)                                     | At runtime                      | Any multivector, stored densely                   | No           |

Two things set Seagrass apart in that list.

* It is the only Rust library here that picks the signature at runtime and stores arbitrary multivectors sparsely. The general-purpose Rust crates fix the algebra at compile time or through a generator, and the SIMD-focused ones offer a fixed menu of types.
* Its closest relative is kingdon, which also exploits sparsity at runtime but generates and compiles code for each new combination of input blades. Seagrass skips the same structural zeros with bit operations instead, so there is no code-generation step on first use, and it runs as native Rust.

Seagrass is not the right tool for every job. If you work in exactly one algebra and need peak throughput, especially PGA 3D, a SIMD library such as geometric_algebra or g3 will be faster. If you need more than 6 dimensions, ganja.js, kingdon, or the Python clifford package handle larger algebras.

## Status and limitations

* Seagrass is pre-1.0 and not yet on crates.io, so expect API changes.
* Dimensions run from 3 to 6, one bitset width per dimension.
* Each value owns a small heap-allocated coefficient vector, so peak throughput trails the SIMD libraries listed in [How it compares](#how-it-compares).
* Mixing metrics in one operation panics.
* Multivectors do not implement `Display` yet. Compare against parsed strings, as the examples do, or walk `iter()`.

## Running the tests

```bash
cargo test
```

Products are checked under Euclidean, projective, and hyperbolic metrics and mixed signatures such as (3, 1, 1), in up to 6 dimensions. Duals are also tested with two null vectors, as in (3, 0, 2). Further tests cover parsing, `exp` and `log`, nalgebra conversions, and PGA camera projections.
