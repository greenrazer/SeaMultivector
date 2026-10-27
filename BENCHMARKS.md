# Benchmarks

How fast Seagrass is, and how it compares with other Rust libraries doing the same work.
Everything here ran on one machine: an AMD Ryzen 9 5900X (12 cores, 24 threads) and an NVIDIA
GeForce RTX 3090, on Windows 11.

* [At a glance](#at-a-glance): Seagrass against each library, as ratios
* [0.1 to 0.2](#01-to-02): single values before and after the refactor
* [Other crates](#other-crates): single values against clifford, geometric_algebra and nalgebra
* [Batches on the GPU and the CPU](#batches-on-the-gpu-and-the-cpu): millions of values at once
* [Batches on the CPU with SIMD](#batches-on-the-cpu-with-simd): the SIMD backend against the
  same libraries

## At a glance

Each table compares Seagrass with other libraries doing the same work. "Faster" and "slower"
describe Seagrass: "3.2× slower" means the other library did the work 3.2 times as fast. "About
even" means within 10%, about how far results moved between runs. Bold marks where Seagrass is
faster, and the sections after this one have the numbers behind every ratio.

### Summary

The fastest other library doing the same work, for four typical tasks. The table leaves out
ultraviolet's rotors, because they only rotate.

| Task | Fastest other library | Seagrass is |
|---|---|---|
| One PGA3D product: motor × motor, or a motor applied to a point | geometric_algebra (`f32`), nalgebra `Matrix4` | 18 to 22× slower |
| One CGA (4,1) product: vector × vector, or a rotor applied to a vector | clifford (dense) | **2.9 to 6.9× faster** |
| Move 100k points on one CPU thread (`SimdBackend`) | glam `Mat4 × Vec4` | 3.2× slower |
| Move 4M points on the GPU, data kept on the GPU (`GpuBackend`) | a hand-written WGSL shader | **1.5× faster** |

### One value at a time

One operation per call on one thread, from `benches/compare.rs` ([numbers](#other-crates)). Both
sides run in `f64`, except geometric_algebra, which only has `f32`. For points, nalgebra's
`Matrix4` multiplies a `Vector4`, and clifford's CGA rows use its dense multivector. A dash means
the library was not measured on that operation.

| Operation | clifford | geometric_algebra | nalgebra `Isometry3` | nalgebra `UnitDualQuaternion` | nalgebra `Matrix4` |
|---|---|---|---|---|---|
| PGA3D motor × motor | 3.1× slower | 18× slower (both `f32`) | 3.8× slower | 6.6× slower | 5.8× slower |
| PGA3D motor applied to a point | 5.0× slower | 3.1× slower (its `f32` against `f64`) | 6.6× slower | 13× slower | 22× slower |
| CGA (4,1) vector × vector | **2.9× faster** | — | — | — | — |
| CGA (4,1) rotor applied to a vector | **6.9× faster** | — | — | — | — |
| 6D bivector × vector | no 6D algebras | — | — | — | — |

### Batches on one CPU thread

Seagrass's `SimdBackend` against each library on one thread, in `f32`, at three batch sizes
([numbers](#batches-on-the-cpu-with-simd)). From 1M elements every library waits on memory, which
narrows most gaps.

| Task | Against | 10k | 100k | 1M |
|---|---|---|---|---|
| Motor applied to a point | glam `Mat4 × Vec4` | 3.5× slower | 3.2× slower | about even |
| | nalgebra `Matrix4 × Vector4` | 1.1× slower | about even | **1.2× faster** |
| | nalgebra `UnitDualQuaternion` | 1.2× slower | 1.2× slower | 1.3× slower |
| | geometric_algebra `Motor` (SSE) | **7.6× faster** | **8.3× faster** | **6.4× faster** |
| | ultraviolet `Rotor3x8` (rotation only) | 3.6× slower | 3.3× slower | 3.0× slower |
| Motor × motor | glam `Mat4` | **1.2× faster** | **1.4× faster** | **1.4× faster** |
| | nalgebra `Matrix4` | **2.1× faster** | **2.5× faster** | **1.6× faster** |
| | nalgebra `UnitDualQuaternion` | **1.3× faster** | **1.5× faster** | 1.2× slower |
| | geometric_algebra `Motor` (SSE) | 2.3× slower | 2.0× slower | 1.5× slower |
| | ultraviolet `Rotor3x8` (rotation only) | 4.0× slower | 3.3× slower | 3.8× slower |
| Frame loop | glam `Mat4` | 1.5× slower | 1.3× slower | 1.7× slower |

### Batches on 24 CPU threads

The faster of Seagrass's two CPU backends against each library, all split across 24 threads with
rayon. Results on 24 threads moved by up to 1.4× between runs, so ratios smaller than that can
flip.

| Task | Against | 10k | 100k | 1M |
|---|---|---|---|---|
| Motor applied to a point | glam `Mat4 × Vec4` | **1.5× faster** | **1.1× faster** | **1.6× faster** |
| | nalgebra `Matrix4 × Vector4` | **1.6× faster** | **1.2× faster** | **1.6× faster** |
| | nalgebra `UnitDualQuaternion` | **1.8× faster** | **1.4× faster** | 1.4× slower |
| | geometric_algebra `Motor` (SSE) | **2.3× faster** | **3.1× faster** | **1.2× faster** |
| | ultraviolet `Rotor3x8` (rotation only) | **1.2× faster** | 1.2× slower | 3.0× slower |
| Motor × motor | glam `Mat4` | **1.3× faster** | **1.2× faster** | **2.3× faster** |
| | nalgebra `Matrix4` | **1.4× faster** | **1.2× faster** | **2.4× faster** |
| | nalgebra `UnitDualQuaternion` | **1.3× faster** | 1.2× slower | about even |
| | geometric_algebra `Motor` (SSE) | **1.1× faster** | 1.3× slower | about even |
| | ultraviolet `Rotor3x8` (rotation only) | 1.1× slower | 1.9× slower | 3.6× slower |
| Frame loop | glam `Mat4` | 1.4× slower | 1.7× slower | 1.3× slower |

### Batches on the GPU

Seagrass's `GpuBackend` applying a motor to each point with the data kept in GPU memory, against
a hand-written WGSL shader that multiplies each point by a 4×4 matrix, and against the fastest CPU
result of any other library on one thread or 24 ([numbers](#batches-on-the-gpu-and-the-cpu)). At
every size that fastest CPU result is ultraviolet's rotors, which only rotate and so do less work.

| Elements | Against the hand-written shader | Against the fastest CPU library |
|---|---|---|
| 1k | 2.7× slower | 41× slower |
| 10k | 2.1× slower | 2.9× slower |
| 100k | 1.1× slower | **2.9× faster** |
| 1M | **1.4× faster** | **7.8× faster** |
| 4M | **1.5× faster** | **21× faster** |

Copying each batch to the GPU and back costs more than the work: with upload and readback,
Seagrass on the GPU managed at most 60 million points per second, against 130 to 172 for Seagrass
on one CPU thread.

## 0.1 to 0.2

Single-value timings from `benches/products.rs`, median of criterion's samples (rustc 1.96, release
profile, `f64`). Run them with:

```bash
cargo bench --bench products
```

Run-to-run noise on this machine is about ±5%; a busy core elsewhere on the machine moves results
by much more, so compare runs taken back to back on an idle machine.

"0.1" is the crate before the refactor, with the corrected `exp`. "0.2" is one generic engine with
cached product plans and inline coefficient storage.

| Case | 0.1 (ns) | 0.2 (ns) | Change |
|---|---:|---:|---:|
| PGA3D motor × motor | 239.1 | 50.7 | −79% |
| PGA3D motor applied to a point (`/`) | 511.2 | 56.7 | −89% |
| PGA3D point & point (join) | 188.4 | 85.0 | −55% |
| PGA3D line ^ plane (meet) | 87.0 | 35.7 | −59% |
| PGA3D exp(line) → motor | 272.4 | 86.8 | −68% |
| PGA3D motor log | 823.5 | 311.5 | −62% |
| PGA3D motor → matrix pair | 16.3 | 17.8 | +9% |
| PGA3D matrix → motor | 1,084.2 | 186.5 | −83% |
| 6D vector × vector (2 × 2 blades) | 101.8 | 37.1 | −64% |
| 6D bivector × vector (15 × 6 blades) | 697.0 | 123.2 | −82% |
| 6D (3,3) 19 × 9 blades, geometric | 1,500.3 | 191.9 | −87% |
| 6D (3,3) 19 × 9 blades, outer | 696.9 | 141.7 | −80% |
| 6D (3,3) 19 × 9 blades, inner | 1,129.3 | 123.6 | −89% |
| 6D (3,3) 19 × 9 blades, regressive | 361.3 | 160.1 | −56% |
| 6D 19 + 9 blades, add | 127.0 | 111.5 | −12% |
| 6D 19 blades, dual | 55.7 | 65.2 | +17% |
| 6D 19 blades, reverse | 66.5 | 60.0 | −10% |
| 6D 19 blades, norm | 49.0 | 33.1 | −32% |
| `from_str`, 4D, 16 terms | 5,841.7 | 289.8 | −95% |
| `from_str`, 3D, 2 terms | 891.6 | 52.9 | −94% |

Where the time went:

* Products look up a plan per (operation, signature, left blades, right blades) and run its list of
  multiply-adds; the bitset work, popcounts and mask building run once per plan instead of on
  every call. A loop that repeats one product finds its plan with a four-word compare.
* The sandwich plan is flattened into trilinear terms with symbolic cancellation: 40 terms for a
  motor applied to a point, where the two products have 68 pairs.
* Values with up to 16 coefficients (every 3D and 4D value) live inline, with no heap allocation.
  Larger 5D and 6D values spill to the heap, which is why the dual of a 19-blade 6D value is
  slower than before: copying spilled storage costs a little more than copying a `Vec`.
* `exp` no longer runs the Taylor series for arguments whose square stores no blades, such as the
  ideal lines behind every translation motor, which is most of the matrix → motor gain.
* The hand-written parser replaced a regex.

## Other crates

From `benches/compare.rs`, run with `cargo bench --bench compare`. Each crate uses its fastest
single-value form: `clifford` and `geometric_algebra` through their generated PGA3D types (and
`clifford`'s dense generic multivector for the conformal rows), nalgebra through isometries, dual
quaternions and matrices. `geometric_algebra` only offers `f32`.

| Operation | Implementation | Time (ns) |
|---|---|---:|
| PGA3D motor × motor | Seagrass `f64` | 44.3 |
| | Seagrass `f32` | 35.5 |
| | clifford 0.3 `f64` | 14.1 |
| | geometric_algebra 0.3 `f32` (SIMD) | 2.0 |
| | nalgebra 0.35 `Isometry3` | 11.6 |
| | nalgebra 0.35 `UnitDualQuaternion` | 6.7 |
| | nalgebra 0.35 `Matrix4` | 7.6 |
| PGA3D motor applied to a point | Seagrass `f64` | 72.9 |
| | clifford 0.3 `f64` | 14.5 |
| | geometric_algebra 0.3 `f32` (SIMD) | 23.2 |
| | nalgebra 0.35 `Isometry3` | 11.0 |
| | nalgebra 0.35 `UnitDualQuaternion` | 5.5 |
| | nalgebra 0.35 `Matrix4` × `Vector4` | 3.3 |
| CGA (4,1) vector × vector | Seagrass `f64` | 49.9 |
| | clifford 0.3 dense `f64` | 143.6 |
| CGA (4,1) rotor applied to a vector | Seagrass `f64` | 66.0 |
| | clifford 0.3 dense `f64` | 457.3 |
| 6D bivector × vector | Seagrass `f64` | 124.8 |
| | clifford 0.3 | no 6D algebras |

Where Seagrass wins:

* Any signature, at compile time or at runtime, up to 6D, from one code path. The conformal and 6D
  rows have no specialized competitor, and the dense general library is 3 to 7 times slower because
  it visits every coefficient pair.
* Sparse values. A product visits only the pairs that can be non-zero, so a 6D value with a handful
  of blades costs a handful of multiply-adds.
* Mixed-grade intermediate results, which a fixed menu of types cannot represent.

Where it does not:

* Peak PGA3D throughput. Code generated for one algebra, SIMD types, and plain matrices are 3 to 22
  times faster for motor operations, because they know the structure at compile time and need no
  plan lookup. If you only ever work in PGA3D and need every nanosecond, use one of them.

The motors here are screws about the z axis with 6 stored blades each. `benches/products.rs`
multiplies a general motor (8 blades) by a rotation (4 blades), which has fewer pairs (50.7 ns),
and applies the general motor to a point (56.7 ns).

## Batches on the GPU and the CPU

From `examples/gpu_bench.rs`, run with
`cargo run --release --example gpu_bench --features gpu,rayon,simd`, on an NVIDIA GeForce RTX 3090
(Vulkan, driver 610.88) and the same Ryzen 9 5900X (24 threads), `f32` throughout. The raw rows
are in [`docs/gpu_bench.csv`](docs/gpu_bench.csv).

Unless a table says otherwise, numbers are millions of elements per second with the data already
where it is computed: in device memory for GPU rows, with 50 operations queued per submission,
and in preallocated output buffers for CPU rows. Each is the median of five samples after a
warm-up. Every competitor runs in its documented fastest form on one thread and under rayon:
glam `Mat4`, nalgebra `Matrix4` and `UnitDualQuaternion`, geometric_algebra's SSE motors, and
ultraviolet's AVX2 `Rotor3x8`, which handles eight rotations per value and only rotates, so its
rows do less work than a rigid motion. The hand-written WGSL shader multiplies a 4×4 matrix by a
point, which is what an engine would write.

Compared with an earlier full run, 80 of the 96 rows at 1M and 4M elements agreed within 10%.
GPU rows from 1k to 100k did not: a dispatch there costs more than its work, so they depend on the
GPU's clock state and varied by up to 2.6× between runs. Read those rows as an order of magnitude.

![Motor applied to a point: elements per second against batch size, for Seagrass on the GPU with the specialised and interpreter shaders, a hand-written WGSL matrix shader, Seagrass on 24 CPU threads, and glam on 1 and 24 threads](docs/gpu-motor-point.svg)

Where the GPU wins:

* Large batches that stay on the GPU. From 100k elements, the specialised shader is faster than
  every CPU row: 2.9 times ultraviolet on 24 threads at 100k (3,938 against 1,342) and 7.8 times at
  1M (11,348 against 1,464).
* Against the hand-written shader, from 1M elements: 1.4 times as fast at 1M and 1.5 times at 4M.
  Both move about 850 GB/s, close to the card's 936 GB/s memory bandwidth, and a motor, a point, and
  the result take 64 bytes per element where a matrix, a point, and the result take 96.
* Sparse storage. Storing all 16 blades of every value takes 192 bytes per element and runs 3 times
  slower at 4M (4,373 against 13,226).
* Many operations in a row. The frame loop composes each of 1M motors with a small motion and moves
  each point, 100 frames with the data on the GPU throughout, at 21 times the speed of glam on 24
  threads (4,528 against 215).
* Any algebra. The conformal, spacetime, and 6D products run 29 to 38 times faster than Seagrass on
  24 CPU threads at 1M, with no algebra-specific code.
* The interpreter shader, which runs any plan without compiling a new shader, is 5.9 times slower
  than the specialised one at 1M, and still faster than every CPU row there (1,919 against 1,464).

Where it does not:

* Small batches. At 10k elements, glam on one thread (1,167) beats every GPU row, and the crossover
  lies between 10k and 100k. A GPU dispatch costs more than its work at these sizes, which is why
  `GpuBackend` keeps batches below 20,000 elements on the CPU by default.
* One-off operations. With upload and readback, Seagrass on the GPU reaches at most 60 million
  elements per second and the hand-written shader 65, against 130 to 172 for Seagrass on one CPU
  thread. The GPU pays only when data stays on it.
* Small batches against the hand-written shader. At 10k, Seagrass's specialised shader runs at half
  its speed (409 against 850).
* Seagrass on the CPU against glam. On one thread up to 100k, glam's matrices are 6 to 8 times
  faster than `CpuBackend` (1,167 against 150 at 10k), and 3.2 to 4.7 times faster than
  `SimdBackend` ([below](#batches-on-the-cpu-with-simd)).

GPU results agree with the CPU to within 3.6 × 10⁻⁷ of the largest coefficient.

### Every workload at 1M elements

The same API and generated shaders in every workload:

| At 1M elements | Seagrass GPU, specialised | Seagrass GPU, interpreted | Seagrass CPU, 24 threads | Fastest other CPU row |
|---|---:|---:|---:|---|
| PGA3D motor × motor | 7,809 | 1,459 | 257 | 987, ultraviolet (rotation only) |
| Frame loop: compose each motor, then move its point | 4,528 | — | 148 | 215, glam on 24 threads |
| CGA (4,1) vector × vector | 7,437 | 1,336 | 215 | none |
| STA (1,3) Lorentz transformation of a vector | 11,523 | 1,770 | 396 | none |
| 6D bivector × vector | 4,166 | 780 | 109 | none |

### Motor × motor

| Implementation | 1k | 10k | 100k | 1M | 4M |
|---|---:|---:|---:|---:|---:|
| Seagrass GPU, specialised | 79 | 781 | 3,917 | 7,809 | 8,372 |
| Seagrass GPU, interpreted | 43 | 413 | 1,297 | 1,459 | 1,543 |
| Seagrass GPU, dense storage | 74 | 727 | 2,736 | 4,072 | 4,414 |
| Seagrass CPU, 1 thread | 180 | 175 | 163 | 104 | 85 |
| Seagrass CPU, 24 threads | 20 | 209 | 446 | 257 | 239 |
| glam Mat4, 1 thread | 296 | 291 | 279 | 119 | 112 |
| glam Mat4, 24 threads | 29 | 164 | 517 | 112 | 117 |
| nalgebra Matrix4, 1 thread | 167 | 160 | 150 | 104 | 100 |
| nalgebra Matrix4, 24 threads | 30 | 170 | 521 | 125 | 115 |
| nalgebra UnitDualQuaternion, 1 thread | 264 | 256 | 254 | 209 | 196 |
| nalgebra UnitDualQuaternion, 24 threads | 30 | 174 | 713 | 291 | 238 |
| geometric_algebra Motor (SSE), 1 thread | 752 | 740 | 733 | 277 | 200 |
| geometric_algebra Motor (SSE), 24 threads | 32 | 197 | 791 | 288 | 242 |
| ultraviolet Rotor3x8 (AVX2, rotation only), 1 thread | 1,322 | 1,273 | 1,289 | 642 | 474 |
| ultraviolet Rotor3x8 (AVX2, rotation only), 24 threads | 52 | 266 | 1,212 | 987 | 501 |

### Motor applied to a point

| Implementation | 1k | 10k | 100k | 1M | 4M |
|---|---:|---:|---:|---:|---:|
| Seagrass GPU, specialised | 29 | 409 | 3,938 | 11,348 | 13,226 |
| Seagrass GPU, interpreted | 34 | 414 | 1,541 | 1,919 | 1,999 |
| Seagrass GPU, dense storage | 68 | 648 | 2,473 | 4,158 | 4,373 |
| Seagrass CPU, 1 thread | 162 | 150 | 172 | 139 | 130 |
| Seagrass CPU, 24 threads | 18 | 283 | 723 | 549 | 375 |
| hand-written WGSL Mat4 × Vec4 | 80 | 850 | 4,370 | 8,272 | 8,867 |
| glam Mat4 × Vec4, 1 thread | 1,172 | 1,167 | 1,047 | 291 | 254 |
| glam Mat4 × Vec4, 24 threads | 32 | 212 | 1,004 | 350 | 282 |
| nalgebra Matrix4 × Vector4, 1 thread | 366 | 356 | 360 | 198 | 220 |
| nalgebra Matrix4 × Vector4, 24 threads | 29 | 197 | 954 | 337 | 281 |
| nalgebra UnitDualQuaternion, 1 thread | 409 | 410 | 248 | 348 | 340 |
| nalgebra UnitDualQuaternion, 24 threads | 31 | 178 | 770 | 757 | 478 |
| geometric_algebra Motor (SSE), 1 thread | 41 | 42 | 43 | 43 | 41 |
| geometric_algebra Motor (SSE), 24 threads | 29 | 131 | 350 | 457 | 389 |
| ultraviolet Rotor3x8 (AVX2, rotation only), 1 thread | 1,194 | 1,191 | 1,035 | 762 | 515 |
| ultraviolet Rotor3x8 (AVX2, rotation only), 24 threads | 53 | 271 | 1,342 | 1,464 | 639 |

### Frame loop (per frame)

| Implementation | 1k | 10k | 100k | 1M | 4M |
|---|---:|---:|---:|---:|---:|
| Seagrass GPU, specialised | 23 | 209 | 1,332 | 4,528 | 5,222 |
| Seagrass CPU, 1 thread | 89 | 81 | 90 | 58 | 55 |
| Seagrass CPU, 24 threads | 9 | 100 | 250 | 148 | 144 |
| glam Mat4, 1 thread | 237 | 236 | 232 | 142 | 123 |
| glam Mat4, 24 threads | 29 | 155 | 443 | 215 | 165 |

### CGA (4,1) vector × vector

| Implementation | 1k | 10k | 100k | 1M | 4M |
|---|---:|---:|---:|---:|---:|
| Seagrass GPU, specialised | 37 | 500 | 2,268 | 7,437 | 8,168 |
| Seagrass GPU, interpreted | 39 | 430 | 1,162 | 1,336 | 1,426 |
| Seagrass CPU, 1 thread | 236 | 251 | 220 | 105 | 99 |
| Seagrass CPU, 24 threads | 18 | 143 | 373 | 215 | 184 |

### STA (1,3) Lorentz transformation of a vector

| Implementation | 1k | 10k | 100k | 1M | 4M |
|---|---:|---:|---:|---:|---:|
| Seagrass GPU, specialised | 75 | 754 | 6,254 | 11,523 | 12,944 |
| Seagrass GPU, interpreted | 36 | 353 | 1,255 | 1,770 | 1,929 |
| Seagrass CPU, 1 thread | 84 | 89 | 85 | 76 | 78 |
| Seagrass CPU, 24 threads | 18 | 229 | 541 | 396 | 397 |

### 6D bivector × vector

| Implementation | 1k | 10k | 100k | 1M | 4M |
|---|---:|---:|---:|---:|---:|
| Seagrass GPU, specialised | 75 | 700 | 2,680 | 4,166 | 4,462 |
| Seagrass GPU, interpreted | 28 | 279 | 683 | 780 | 835 |
| Seagrass CPU, 1 thread | 95 | 89 | 90 | 46 | 41 |
| Seagrass CPU, 24 threads | 16 | 55 | 240 | 109 | 99 |

### With upload and readback

| Workload | Implementation | 1k | 10k | 100k | 1M | 4M |
|---|---|---:|---:|---:|---:|---:|
| motor × motor | Seagrass GPU, specialised | 7 | 36 | 33 | 32 | 21 |
| motor × motor | Seagrass GPU, interpreted | 7 | 36 | 32 | 30 | 20 |
| motor × motor | Seagrass GPU, dense storage | 7 | 11 | 17 | 16 | 10 |
| motor applied to a point | Seagrass GPU, specialised | 5 | 24 | 47 | 60 | 46 |
| motor applied to a point | Seagrass GPU, interpreted | 7 | 45 | 55 | 59 | 47 |
| motor applied to a point | Seagrass GPU, dense storage | 6 | 22 | 17 | 16 | 10 |
| motor applied to a point | hand-written WGSL Mat4 × Vec4 | 5 | 32 | 65 | 56 | 35 |
| CGA (4,1) vector × vector | Seagrass GPU, specialised | 7 | 30 | 20 | 23 | 17 |
| CGA (4,1) vector × vector | Seagrass GPU, interpreted | 7 | 20 | 23 | 22 | 16 |
| STA (1,3) Lorentz transformation of a vector | Seagrass GPU, specialised | 8 | 46 | 55 | 58 | 41 |
| STA (1,3) Lorentz transformation of a vector | Seagrass GPU, interpreted | 7 | 44 | 52 | 57 | 47 |
| 6D bivector × vector | Seagrass GPU, specialised | 6 | 10 | 12 | 10 | 8 |
| 6D bivector × vector | Seagrass GPU, interpreted | 6 | 10 | 12 | 9 | 9 |

### GPU timestamps

| Workload | Implementation | 1k | 10k | 100k | 1M | 4M | Largest difference from the CPU |
|---|---|---:|---:|---:|---:|---:|---:|
| motor × motor | Seagrass GPU, specialised | 3.136 | 0.298 | 0.130 | 0.115 | 0.111 | 2.4e-07 |
| motor × motor | Seagrass GPU, interpreted | 12.480 | 1.341 | 0.640 | 0.663 | 0.641 | 0.0e+00 |
| motor × motor | Seagrass GPU, dense storage | 3.968 | 0.368 | 0.251 | 0.232 | 0.222 | 1.8e-07 |
| motor applied to a point | Seagrass GPU, specialised | 8.736 | 0.758 | 0.114 | 0.074 | 0.072 | 3.6e-07 |
| motor applied to a point | Seagrass GPU, interpreted | 17.408 | 1.270 | 0.525 | 0.497 | 0.493 | 2.8e-07 |
| motor applied to a point | Seagrass GPU, dense storage | 4.512 | 0.390 | 0.263 | 0.226 | 0.222 | 3.6e-07 |
| CGA (4,1) vector × vector | Seagrass GPU, specialised | 5.248 | 0.522 | 0.138 | 0.122 | 0.117 | 1.8e-07 |
| CGA (4,1) vector × vector | Seagrass GPU, interpreted | 10.944 | 1.158 | 0.722 | 0.729 | 0.695 | 0.0e+00 |
| STA (1,3) Lorentz transformation of a vector | Seagrass GPU, specialised | 3.040 | 0.269 | 0.063 | 0.075 | 0.072 | 2.8e-07 |
| STA (1,3) Lorentz transformation of a vector | Seagrass GPU, interpreted | 16.032 | 1.654 | 0.657 | 0.536 | 0.518 | 2.3e-07 |
| 6D bivector × vector | Seagrass GPU, specialised | 3.296 | 0.342 | 0.250 | 0.225 | 0.220 | 1.8e-07 |
| 6D bivector × vector | Seagrass GPU, interpreted | 22.528 | 2.250 | 1.347 | 1.247 | 1.194 | 0.0e+00 |

GPU timestamps are nanoseconds per element for one dispatch, from the GPU's own timer. The last
column is the largest difference from the `f32` CPU result over all sizes, relative to the
largest CPU coefficient.

## Batches on the CPU with SIMD

From `examples/gpu_bench.rs` without its GPU rows, run with
`cargo run --release --example gpu_bench --features gpu,rayon,simd -- --cpu`, on the same
Ryzen 9 5900X (24 threads), where `SimdBackend` picks AVX2 with FMA, `f32` throughout. The raw
rows are in [`docs/cpu_bench.csv`](docs/cpu_bench.csv). The method and the competitors are those
of the GPU section above, and every row in these tables comes from the same run, which is newer
than the GPU tables' run.

Compared with a second run, 191 of the 230 rows agreed within 10%, and 110 of the 115 one-thread
rows did; the largest one-thread change was nalgebra's `UnitDualQuaternion` motor × motor at 100k,
by 1.5×. Rows on 24 threads varied by up to 1.4×.

Where SIMD helps:

* One thread, any algebra. Every workload at every size ran 1.3 to 2.3 times as fast as with
  `CpuBackend`, in each of two runs, with no algebra-specific code.

Where it does not:

* Many threads. Spread over 24 threads, a batch soon waits on memory rather than arithmetic. At 100k
  elements SIMD ran 0.9 to 1.7 times as fast as `CpuBackend` across two runs (STA: 919 against 532),
  and from 1M elements 0.8 to 1.0 times.
* Code built for one algebra. On one thread up to 100k, glam's matrices still move points 3.2 to 4.7
  times as fast (1,151 against 331 at 10k), and ultraviolet's AVX2 rotors, which only rotate,
  compose 3.2 to 5.4 times as fast. The SIMD kernels interpret a plan, so every multiply-add loads
  both of its operands, where code built for one algebra keeps them in registers.
* Bit-for-bit agreement. Fused multiply-adds and a different summation order change the last bits.
  Relative to the sum of the absolute terms, the tests allow differences from `CpuBackend` of 4
  machine epsilons and see at most 0.4.

### One thread, every workload

| One thread | CPU at 10k | SIMD at 10k | CPU at 1M | SIMD at 1M |
|---|---:|---:|---:|---:|
| PGA3D motor × motor | 189 | 341 | 103 | 168 |
| Motor applied to a point | 174 | 331 | 145 | 277 |
| Frame loop: compose each motor, then move its point | 90 | 156 | 60 | 95 |
| CGA (4,1) vector × vector | 251 | 416 | 109 | 151 |
| STA (1,3) Lorentz transformation of a vector | 87 | 153 | 81 | 182 |
| 6D bivector × vector | 101 | 165 | 46 | 71 |

### Motor × motor

| Implementation | 1k | 10k | 100k | 1M | 4M |
|---|---:|---:|---:|---:|---:|
| Seagrass CPU, 1 thread | 185 | 189 | 190 | 103 | 99 |
| Seagrass CPU, 24 threads | 18 | 190 | 554 | 291 | 238 |
| Seagrass SIMD, 1 thread | 256 | 341 | 391 | 168 | 149 |
| Seagrass SIMD, 24 threads | 19 | 207 | 565 | 263 | 236 |
| glam Mat4, 1 thread | 296 | 295 | 283 | 123 | 108 |
| glam Mat4, 24 threads | 27 | 155 | 488 | 127 | 117 |
| nalgebra Matrix4, 1 thread | 168 | 161 | 155 | 104 | 107 |
| nalgebra Matrix4, 24 threads | 28 | 153 | 479 | 124 | 118 |
| nalgebra UnitDualQuaternion, 1 thread | 271 | 265 | 262 | 206 | 204 |
| nalgebra UnitDualQuaternion, 24 threads | 27 | 159 | 669 | 286 | 237 |
| geometric_algebra Motor (SSE), 1 thread | 815 | 793 | 764 | 250 | 213 |
| geometric_algebra Motor (SSE), 24 threads | 28 | 181 | 720 | 287 | 237 |
| ultraviolet Rotor3x8 (AVX2, rotation only), 1 thread | 1,375 | 1,347 | 1,304 | 637 | 463 |
| ultraviolet Rotor3x8 (AVX2, rotation only), 24 threads | 47 | 234 | 1,096 | 1,056 | 506 |

### Motor applied to a point

| Implementation | 1k | 10k | 100k | 1M | 4M |
|---|---:|---:|---:|---:|---:|
| Seagrass CPU, 1 thread | 172 | 174 | 173 | 145 | 123 |
| Seagrass CPU, 24 threads | 20 | 243 | 923 | 557 | 406 |
| Seagrass SIMD, 1 thread | 254 | 331 | 359 | 277 | 252 |
| Seagrass SIMD, 24 threads | 18 | 293 | 1,052 | 508 | 393 |
| glam Mat4 × Vec4, 1 thread | 1,206 | 1,151 | 1,152 | 293 | 261 |
| glam Mat4 × Vec4, 24 threads | 29 | 193 | 927 | 345 | 282 |
| nalgebra Matrix4 × Vector4, 1 thread | 374 | 366 | 367 | 225 | 224 |
| nalgebra Matrix4 × Vector4, 24 threads | 29 | 183 | 914 | 347 | 281 |
| nalgebra UnitDualQuaternion, 1 thread | 414 | 412 | 413 | 354 | 340 |
| nalgebra UnitDualQuaternion, 24 threads | 28 | 167 | 729 | 786 | 483 |
| geometric_algebra Motor (SSE), 1 thread | 44 | 43 | 43 | 43 | 43 |
| geometric_algebra Motor (SSE), 24 threads | 28 | 126 | 340 | 454 | 388 |
| ultraviolet Rotor3x8 (AVX2, rotation only), 1 thread | 1,199 | 1,201 | 1,178 | 824 | 562 |
| ultraviolet Rotor3x8 (AVX2, rotation only), 24 threads | 46 | 237 | 1,258 | 1,686 | 654 |

### Frame loop (per frame)

| Implementation | 1k | 10k | 100k | 1M | 4M |
|---|---:|---:|---:|---:|---:|
| Seagrass CPU, 1 thread | 87 | 90 | 90 | 60 | 57 |
| Seagrass CPU, 24 threads | 9 | 97 | 217 | 159 | 144 |
| Seagrass SIMD, 1 thread | 136 | 156 | 182 | 95 | 93 |
| Seagrass SIMD, 24 threads | 9 | 93 | 254 | 136 | 142 |
| glam Mat4, 1 thread | 234 | 236 | 236 | 164 | 150 |
| glam Mat4, 24 threads | 25 | 134 | 424 | 209 | 168 |

### CGA (4,1) vector × vector

| Implementation | 1k | 10k | 100k | 1M | 4M |
|---|---:|---:|---:|---:|---:|
| Seagrass CPU, 1 thread | 245 | 251 | 252 | 109 | 101 |
| Seagrass CPU, 24 threads | 19 | 100 | 362 | 237 | 199 |
| Seagrass SIMD, 1 thread | 310 | 416 | 446 | 151 | 133 |
| Seagrass SIMD, 24 threads | 19 | 136 | 383 | 222 | 197 |

### STA (1,3) Lorentz transformation of a vector

| Implementation | 1k | 10k | 100k | 1M | 4M |
|---|---:|---:|---:|---:|---:|
| Seagrass CPU, 1 thread | 87 | 87 | 90 | 81 | 78 |
| Seagrass CPU, 24 threads | 19 | 224 | 532 | 524 | 384 |
| Seagrass SIMD, 1 thread | 145 | 153 | 194 | 182 | 176 |
| Seagrass SIMD, 24 threads | 19 | 282 | 919 | 443 | 388 |

### 6D bivector × vector

| Implementation | 1k | 10k | 100k | 1M | 4M |
|---|---:|---:|---:|---:|---:|
| Seagrass CPU, 1 thread | 101 | 101 | 88 | 46 | 45 |
| Seagrass CPU, 24 threads | 18 | 61 | 239 | 110 | 95 |
| Seagrass SIMD, 1 thread | 145 | 165 | 158 | 71 | 68 |
| Seagrass SIMD, 24 threads | 18 | 78 | 249 | 107 | 96 |

`SimdBackend` sums each output coefficient in registers for a group of vectors of elements at once
(8 vectors for products, 4 for sandwiches), with each product output's negative entries summed
before its positive ones so that no entry branches on its sign. Because the plan is interpreted,
every multiply-add reads both of its operands from the planes, where code generated for one
algebra keeps them in registers. A one-off version of the motor × motor loop written with raw
AVX2 intrinsics and compiled for the native CPU ran at 13 to 17 cycles per element, against 14 for
`SimdBackend` at 10k, so the portable SIMD code costs nothing measurable over intrinsics.
