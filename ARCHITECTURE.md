# Architecture

How Seagrass works inside. The [readme](readme.md) covers using the library; this document follows a product from the bit patterns that name blades to the cached plan that computes it, then covers batches and their backends. [BENCHMARKS.md](BENCHMARKS.md) has the measurements.

* [Blades are bit patterns](#blades-are-bit-patterns)
* [A multivector is one integer and a short vector](#a-multivector-is-one-integer-and-a-short-vector)
* [The trick: multiply a whole set of blades at once](#the-trick-multiply-a-whole-set-of-blades-at-once)
* [Finding every pair with one AND](#finding-every-pair-with-one-and)
* [Outer products and contractions take one shift](#outer-products-and-contractions-take-one-shift)
* [Plans: the structure once, the numbers every time](#plans-the-structure-once-the-numbers-every-time)
* [One table of masks for every dimension](#one-table-of-masks-for-every-dimension)
* [Unary operations are word-sized](#unary-operations-are-word-sized)
* [Sandwiches and matrices](#sandwiches-and-matrices)
* [Batches and backends](#batches-and-backends)

## Blades are bit patterns

Each basis blade is indexed by a bitmask of its basis vectors, so the blade produced by multiplying two blades is the XOR of their indices.

```text
blade   1     e1    e2    e12   e3    e13   e23   e123
index   0     1     2     3     4     5     6     7
bits    000   001   010   011   100   101   110   111
```

## A multivector is one integer and a short vector

A second bitset with one bit per blade (8 bits in 3D, 64 in 6D) records which blades are present. Coefficients are packed in index order, and a blade's slot is the number of set bits below it, which is one `count_ones` call. Up to 16 coefficients live inline in the value, so no 3D or 4D value touches the heap.

```text
3 + 2e1 - e23                   3D, so the bitset is a u8

bitset         0b0100_0011      blades 0 (1), 1 (e1), and 6 (e23)
coefficients   [3.0, 2.0, -1.0]
```

## The trick: multiply a whole set of blades at once

Multiplying a blade by the basis vector eₖ flips bit k − 1 of its index. Applied to every blade in a bitset, that flip always reorders the bitset's bits the same way: neighbouring groups of 2ᵏ⁻¹ bits trade places. One mask, two shifts, and an OR do it inside a single machine word, however many blades the set holds.

```text
{1, e12, e3} × e2: neighbouring groups of 2 bits trade places

blade      1   e1     e2  e12      e3  e13    e23 e123
before    [x   .]    [.   x]      [x   .]    [.   .]      {1, e12, e3}
after     [.   x]    [x   .]      [.   .]    [x   .]      {e1, e2, e23}
```

```rust
// CONTAINS[k] marks the 3D blades that hold e_k, such as 0b1100_1100 for e2.
const CONTAINS: [u8; 4] = [0, 0b1010_1010, 0b1100_1100, 0b1111_0000];

// Multiply every blade in `set` by e_k. Blades that hold e_k move down and lose it;
// the rest move up and gain it.
fn times_basis_vector(set: u8, k: usize) -> u8 {
    let shift = 1 << (k - 1);
    ((set & !CONTAINS[k]) << shift) | ((set & CONTAINS[k]) >> shift)
}

// {1, e12, e3} × e2 = {e1, e2, e23}, as in the picture above.
assert_eq!(times_basis_vector(0b0001_1001, 2), 0b0100_0110);
```

The bitset records only which blades are present, so signs are tracked separately (below). Multiplying by a blade with several basis vectors applies one swap per vector, so any blade takes at most n swaps (6 in 6D). The real code is `BladeSet::times_vector` in [`src/blade/set.rs`](src/blade/set.rs).

## Finding every pair with one AND

In a product A × B, a stored blade x of A lands on output blade o exactly when its partner x ⊕ o is stored in B. Multiplying B's bitset by o with the swaps above gives every such x at once, and ANDing the result with A's bitset leaves exactly the pairs that land on o.

```text
Rotor × vector in Euclidean 3D: the pairs that land on e1

blade          1   e1  e2  e12 e3  e13 e23 e123
vector B       .   x   x   .   x   .   .   .
B × e1         x   .   .   x   .   x   .   .     groups of 1 bit trade places
rotor A        x   .   .   x   .   x   x   .
AND            x   .   .   x   .   x   .   .     pairs (1, e1), (e12, e2), (e13, e3)
```

Two more masks complete the pair list:

* Null vectors. A pair whose blades share a null basis vector multiplies to zero. The pairs that survive are those whose left blade x uses only vectors in o or non-null vectors, so one more AND with the sub-blades of (o ∪ non-null vectors) removes the rest before any coefficient is read, which matters in PGA. That set of sub-blades is built by doubling: start from {1} and, for each vector, OR in a copy of the set shifted by that vector's index. Basis vectors are ordered positive, then negative, then null, so the non-null vectors are always the first p + q, and testing whether eₖ is null is the single comparison k > p + q.
* Signs. With o fixed, a pair is identified by its left blade x, so every sign fits in one mask laid out like the bitset. The reordering sign of x × (x ⊕ o) works out to the reversion sign of x, flipped once more for every pair of a vector in o and a higher vector in x. Both parts are precomputed masks, so the reordering signs for o are one constant XORed with one constant per basis vector in o. Each basis vector that squares to −1 and is not in o then flips the sign of every left blade that contains it, because both blades of such a pair contain it.

The outer product and the contractions use a cheaper variant, described next, and the regressive product is an outer product between duals.

## Outer products and contractions take one shift

When every vector of x is also in o, subtracting x from o never borrows, so o − x and o ⊕ x are the same number. In the outer product and the contractions, one blade of each pair always sits inside another, so every pair is the same numeric distance apart, and a single shift of B's bitset lines all of them up. No swaps are needed.

* Left contraction. x ⌋ y lands on o when y = x + o and x shares no vectors with o. Shifting B right by o moves every such y onto its partner x, and a mask keeps only blades without o's vectors.
* Right contraction. The same with a left shift, keeping only blades that contain all of o.
* Outer product. x ∧ y lands on o when y = o − x, which reflects positions instead of shifting them. Reversing B's bits first sends y to position MAX − y, where MAX = 2ⁿ − 1 is the pseudoscalar's index, and that turns the reflection into a right shift by MAX − o. A mask then keeps the sub-blades of o.

```text
Left contraction onto o = e2: shift B right by 2, keep blades without e2

blade          1   e1  e2  e12 e3  e13 e23 e123
B              .   .   .   x   .   x   x   x
B >> 2         .   x   .   x   x   x   .   .
without e2     x   x   .   .   x   x   .   .
result         .   x   .   .   x   x   .   .     pairs e1 ⌋ e12, e3 ⌋ e23, e13 ⌋ e123
```

e13 shifts onto e12 and is masked out, because e13 has no e2 to contribute. ANDing the result with A's bitset gives the pair mask, as before.

## Plans: the structure once, the numbers every time

A product runs in two phases.

1. Structure. The output bitset is the union of B's bitset multiplied by each stored blade of A, skipping pairs that share a null vector. Each output blade then gets its pair mask and sign mask as above.
2. Numbers. For each output blade, multiply the coefficients of its pairs, negate where the sign mask says so, and add.

The structure depends only on the operation, the signature, and which blades the two operands hold, not on their coefficients. So the first phase runs once per such combination and is recorded as a plan: the output bitset and, for each output blade, a list of (left slot, right slot, sign). The plan is cached per thread, and later products with the same structures skip straight to the second phase: one loop of multiply-adds over coefficient slots, with no popcounts or masks. A loop that repeats one product finds its plan with a four-word comparison. Each output coefficient is a fixed sum that no other output touches.

The sandwich M x M̃ gets its own plan. It is derived from the two products M x and (M x) M̃, with M̃'s signs folded in as described [below](#sandwiches-and-matrices), and then expanded into terms m_i x_j m_k. Terms (i, j, k) and (k, j, i) are equal because coefficients commute, so they merge, and terms that cancel disappear before the plan is stored. A PGA3D motor applied to a point needs 40 such terms where the two products have 68 pairs.

## One table of masks for every dimension

Every grade, sign, and "contains eₖ" mask is a single `u64` written for 6D. A blade's index means the same thing in every dimension, so the 3D mask is the low 8 bits of the 6D one, and `as u8` produces it. The odd-grade mask is the Thue–Morse sequence (bit i is the parity of i's bit count), and the even-grade mask, `0x9669_6996_6996_9669`, is its complement. The reordering-sign masks are XORs of the "contains eₖ" masks, and the dual and complement signs, which depend on the dimension's pseudoscalar, come from the reordering sign of each blade and its complement. Both are computed at compile time, and tests check every mask against its definition.

## Unary operations are word-sized

* Reversion, grade involution, and Clifford conjugation flip signs under one constant mask.
* The dual is a bit reversal. A blade's complement has index MAX − i, with every bit flipped, and because the bitset holds exactly 2ⁿ bits, `reverse_bits()` sends position i to MAX − i. Coefficients are stored in index order, so reversing them puts each one in its complement's slot, and a sign mask finishes the job. There is no multiplication by the pseudoscalar, which is why the dual also works in degenerate metrics such as PGA, where the pseudoscalar has no inverse.
* `exp()` uses a closed form (cosine and sine, hyperbolic cosine and sine, or 1 + x) when its argument squares to a scalar. Otherwise it halves the argument until the sum of its absolute coefficients is at most ½, which bounds every power of the argument, runs a 14-term Taylor series that is then exact to `f64` precision, and squares the result back. `log()` on PGA3D motors is closed form.

## Sandwiches and matrices

* The sandwich M x M̃ never builds M̃. In the second product, each blade z of M x pairs with z ⊕ o in M. Multiplying the reversion sign mask by o with the same swap moves the reverse sign of z ⊕ o onto position z, so one XOR folds M̃'s signs into the pair's sign mask.
* A motor's matrix and its inverse come from one expression. Reversing a motor negates its bivector coefficients and keeps its scalar and e1234 coefficients. Each matrix entry is a sum of terms that keep their sign under reversal (constants and products of two bivector coefficients) and terms that flip it (a scalar or e1234 coefficient times a bivector coefficient). Grouping the first kind as A and the second as B gives M = A + B, and because a unit motor's inverse is its reverse, M⁻¹ = A − B.

## Batches and backends

A `Batch` holds many values that store the same blades, kept as one array per blade ("structure of arrays"): the coefficient of blade slot k of element i sits at index `k * len + i`. An elementwise product looks up its plan once for the whole batch and runs it over the arrays, so there is no per-element plan lookup. Values that store fewer blades are padded with zeros.

* `CpuBackend` runs the plan in blocks of 256 elements, small enough that a block of every array stays in the CPU cache. With the `rayon` feature and `CpuBackend::parallel()`, the elements are split into a few chunks per core, each a whole number of blocks. Every element sums its pairs in plan order, as the single-value product does, so the results agree with single-value products bit for bit.
* `SimdBackend` runs the same plans with explicit SIMD instructions for `f32` and `f64`, on the widest instruction set that fearless_simd finds at run time. It takes elements a group of SIMD vectors at a time (8 vectors for products, 4 for sandwiches), keeps each output coefficient's sums for the group in registers, and sums each output's negative pairs before its positive ones so that no pair branches on its sign. Fused multiply-adds and that summation order change the last bits, so results agree with `CpuBackend` to within rounding. Other coefficient types run the `CpuBackend` loops.
* `GpuBackend` keeps batches in device memory and queues products into one command encoder. Only `download` and `finish` submit and wait, so a chain of products never leaves the GPU. Every element is one GPU thread. By default (`ShaderMode::Specialized`) each plan becomes its own WGSL shader with every multiply-add written out, compiled the first time a pair of structures is seen and then cached. `ShaderMode::Interpreted` uses one shader that reads the plan from a buffer, which never compiles anything new but runs several times slower. Batches smaller than `GpuOptions::min_batch` (20,000 elements by default) stay on the CPU, where launching a dispatch would cost more than the work.
