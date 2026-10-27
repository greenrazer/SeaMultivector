# Rust Code Standards

Rules for writing, reviewing and maintaining Rust code in this project. They are
written so that a new contributor or an AI coding agent can follow them without further
context. Rules are imperative. *Must* rules are enforced by CI or review; *prefer* rules
are the default unless a comment gives a reason to deviate.

Where this document is silent, follow, in order: rustfmt defaults, Clippy defaults, the
[Rust API Guidelines](https://rust-lang.github.io/api-guidelines/), and the conventions
of the standard library.

Contents: 1 Scope · 2 Formatting and tooling · 3 Layout · 4 Naming · 5 API design ·
6 Errors and panics · 7 Documentation · 8 Testing · 9 Dependencies · 10 Performance ·
11 Checklist

---

## 1. Purpose and scope

- Applies to every Rust file in the project: library code, unit tests, integration
  tests, doctests, examples, benchmarks and build scripts.
- When adding to an existing module, mirror the structure of the nearest sibling that does
  the same kind of thing. If that sibling violates these standards, do not copy the
  violation; fix it as a separate change or leave it untouched.
- Prefer the simplest construct that is correct; optimise only with a measurement (§10).

## 2. Formatting and tooling

### Formatting

- Format with rustfmt using the default configuration. Do not add `rustfmt.toml` options
  other than `edition`. Run `cargo fmt` before handing code over; CI runs
  `cargo fmt --check`.
- Indent with 4 spaces, never tabs. Maximum line width is 100 columns (rustfmt default).
- Source files use LF line endings; configure editors accordingly.
- Every file ends with exactly one newline. No trailing whitespace. At most one blank
  line in a row.
- Imports: three groups separated by one blank line, each sorted by rustfmt:
  1. `std`, `core`, `alloc`
  2. external crates
  3. `crate::`, `super::`, `self::`
- Let rustfmt lay out long signatures, bound lists and multi-line calls. Do not
  hand-align anything. `#[rustfmt::skip]` requires a comment saying why.
- Keep formatting-only changes separate from behavioural changes, so each can be reviewed
  on its own.

### Lints

- `cargo clippy --all-targets --all-features -- -D warnings` must pass in CI.
- Configure lints once in `Cargo.toml` under `[lints.rust]` / `[lints.clippy]`, not with
  scattered crate attributes. Enable at least: `missing_docs = "warn"`,
  `rust_2018_idioms = "warn"`, `unreachable_pub = "warn"`, `unsafe_op_in_unsafe_fn = "deny"`.
  Crates that need no `unsafe` add `unsafe_code = "forbid"`.
- Never silence a whole crate with `#![allow(dead_code)]`, `#![allow(unused_imports)]` or
  similar. Remove the dead code or the unused import.
- A local `#[allow(...)]` goes on the smallest possible item and carries a comment stating
  the reason:
  ```rust
  // The RHS is used only for its position; see the docs on `insert_at`.
  #[allow(clippy::suspicious_op_assign_impl)]
  fn add_assign(&mut self, rhs: &Entry<T>) { /* ... */ }
  ```

### Toolchain and CI

- Use the newest stable edition. Migrate editions as a dedicated change with
  `cargo fix --edition`.
- Declare `rust-version` (MSRV) in `Cargo.toml`. CI tests both stable and the MSRV.
  Raising the MSRV is a minor-version change and is noted in the changelog.
- CI runs `cargo fmt --check`, the Clippy command above,
  `cargo test --all-targets --all-features`, and
  `RUSTDOCFLAGS="-D warnings" cargo doc --no-deps`.

## 3. Project and module layout

- Standard Cargo layout: `src/lib.rs` (or `main.rs`), `tests/` for integration tests,
  `benches/` for benchmarks, `examples/` for runnable examples.
- Module files use the `name.rs` + `name/` convention. Use `mod.rs` only for a directory
  that is a pure grouping of submodules: it contains `mod` declarations (and at most a
  `use super::*;`), nothing else. Never create both `name.rs` and `name/mod.rs`.
- Every source file is reachable through exactly one `mod` declaration. No orphan files.
- One primary type (or one trait family) per file. When a type has a large operation
  surface, split it: **one operation per file**, each file holding that operation's
  inherent methods, trait impls, operator impls and unit tests. Group operation files into
  directories by category (for example `linear_ops/`, `scalar_ops/`, `conversions/`,
  `iterators/`).
- A file that only continues a type's impl blocks may start with `use super::*;` to
  inherit the parent's imports. This is the only place a glob import of the parent module
  is allowed outside test modules. All other files import explicitly.
- Order inside a type's main file: file doc comment, imports, `mod` declarations, `pub use`
  re-exports, type aliases and constants, the type definition, constructors and core
  inherent impl, trait impls, `#[cfg(test)] mod tests`.
- The crate root holds the crate docs, lint configuration, `mod` declarations and a
  `pub mod prelude` that re-exports the public API. Keep internal module paths private and
  flatten public names with `pub use` so files can move without breaking users.
- Error types live in one dedicated module (for example `error.rs`).

## 4. Naming

- Casing follows the API Guidelines (C-CASE): `UpperCamelCase` for types, traits and
  enum variants (acronyms are words: `HttpClient`, `Vec3d`, not `HTTPClient`, `Vec3D`);
  `snake_case` for functions, methods, modules, variables and files;
  `SCREAMING_SNAKE_CASE` for constants and statics; single capital letters (`T`, `I`) or
  short `UpperCamelCase` words for generic parameters; short lowercase lifetimes (`'a`).
- Do not prefix type names with the crate or project name; modules provide namespacing.
  If the public API already carries a prefix, apply it to every public type or remove it
  everywhere in one release. Never half-apply it.
- Spell identifiers correctly. Check new public names before they ship: renaming a public
  item later is a breaking change.
- Use full words in public APIs. Short names (`idx`, `len`, `lhs`, `rhs`, `tmp`) are fine
  for locals in short functions.
- Verb forms distinguish in-place from returning operations:

  | Kind | Form | Example |
  |---|---|---|
  | mutate in place | imperative verb, `&mut self`, returns `()` | `normalize()`, `negate()`, `sort()` |
  | return a transformed copy | past participle, `&self` | `normalized()`, `negated()`, `sorted()` |
  | predicate | `is_`, `has_`, `contains_` | `is_empty()`, `has_key()` |
  | accessor | field name, no `get_` | `len()`, `inner()` |
  | mutable accessor or variant | `_mut` suffix | `inner_mut()`, `iter_mut()` |
  | iteration | `iter()`, `iter_mut()`, `into_iter()`, `<subset>_iter()` | `row_iter()` |
  | conversion (C-CONV) | `as_` free borrow, `to_` costly borrow, `into_` consuming | `as_slice()`, `to_string()`, `into_inner()` |
  | constructor | `new(...)` with all required state; named alternatives | `new`, `identity()`, `zero()`, `with_capacity(n)`, `from_parts(a, b)` |

- Prefer `From` / `TryFrom` / `FromStr` impls over `from_x` / `into_x` methods when the
  source is a single value. Use `from_<description>(a, b, c)` only for constructors that
  take several parts.
- A trait of in-place variants takes the base trait's name plus `Mut`
  (`Normalize` / `NormalizeMut`).
- Error types are named `<Operation>Error` (`ParseError`); variants name the cause
  (`MissingSeparator`, `IndexOutOfRange`), not the operation.
- Test functions describe the behaviour: `parses_negative_numbers`,
  `rejects_mismatched_lengths`, `round_trips_through_string`. No `test_` prefix, no
  `works` / `_works` names.

## 5. API design

### Ownership and borrowing

- Non-consuming methods take `&self`; mutating methods take `&mut self`; take `self` only
  when ownership moves into the result.
- Parameters: `&T` for read-only access to non-`Copy` values, `T` by value for `Copy`
  types, `&[T]` not `&Vec<T>`, `&str` not `&String`.
- Return owned values from computations. Return `&T` or `Option<&T>` for lookups into
  internal storage, and `&[T]` rather than `&Vec<T>`.
- Derives, in this order when present:
  `Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Default`. Every public type
  derives `Debug`. Derive `Copy` for small plain data, `Eq` whenever the type contains no
  floats, `Default` only when a meaningful default exists.

### Operator overloading

- Overload an operator only when its meaning for the type is conventional in the domain,
  and document the full operator table once in the crate docs.
- Every operator impl forwards to a named method (an inherent method or an operation
  trait's method). The named method is the single source of truth; operator impls contain
  no logic.
- For a binary operator between `A` and `B`, implement all four combinations:
  `A op B`, `&A op B`, `A op &B`, `&A op &B`. Exactly one is the primitive, normally
  `&A op &B`, or `op_assign(&mut self, rhs: &B)` when updating in place avoids an
  allocation. The others delegate to it:
  ```rust
  impl<T: Scalar> Add<&Self> for &Vector<T> {
      type Output = Vector<T>;
      fn add(self, rhs: &Self) -> Vector<T> { /* the one real implementation */ }
  }
  impl<T: Scalar> Add<Self> for Vector<T> {
      type Output = Vector<T>;
      fn add(self, rhs: Self) -> Vector<T> { &self + &rhs }
  }
  ```
- Implement the matching `*Assign` trait whenever the operation can update in place.
- Write `type Output` explicitly. Name the concrete type when it differs from `Self`.
- Never implement a trait method whose body is `unimplemented!()`, `todo!()` or
  `panic!()` for a particular type pair. If the operation is undefined for that pair, do
  not implement the trait; let the compiler reject the call.

### Conversions and constructors

- Implement `From<X> for Y`, never `Into<Y> for X` directly (the blanket impl supplies
  `Into`). Provide `From<&X>` alongside `From<X>` when the borrowed form avoids a clone.
- Fallible conversions use `TryFrom` or `FromStr` with the crate's error type.
- `new` validates invariants. A type whose invariants matter keeps its fields private and
  is constructed only through `new`, named constructors or `From`.

### In-place and returning variants

- Offer both forms when the operation is cheap to perform in place. Implement the
  returning form as clone-then-mutate unless a faster direct path exists:
  ```rust
  pub fn negated(&self) -> Self {
      let mut out = self.clone();
      out.negate();
      out
  }
  ```
- Keep in-place variants in the sibling `…Mut` trait so immutable types can omit them.

### Generics and trait bounds

- Name a recurring bound set once with an alias trait and a blanket impl, then use only
  that trait as the bound:
  ```rust
  pub trait Scalar: Copy + Debug + Add<Output = Self> + Mul<Output = Self> {}
  impl<T: Copy + Debug + Add<Output = Self> + Mul<Output = Self>> Scalar for T {}
  ```
- Short bounds go inline (`impl<T: Scalar>`). Bounds that span more than one line, or
  involve higher-ranked lifetimes or associated types, go in a `where` clause.
- Add `'static` or other lifetime bounds only when an external API requires them.
- Model operations as small single-purpose traits with one method and an associated
  `Output` type. Compose them into an umbrella trait for the full capability, listing
  supertraits one per line. Factor a long repeated `where` block into a helper trait
  rather than copying it across several traits.

### Visibility and re-exports

- Default to private. Use `pub(crate)` for items shared between modules of the crate and
  `pub` only for the supported public API. Enable the `unreachable_pub` lint.
- Struct fields are private unless the type is a plain data record with no invariants.
  Inner fields of newtypes are private; cross-module construction uses `pub(crate)`.
- Do not implement `Deref` / `DerefMut` on a newtype to emulate inheritance (C-DEREF).
  Provide `inner()`, `inner_mut()`, `into_inner()` or `AsRef`, and implement the traits
  the wrapper actually needs. Reserve `Deref` for smart-pointer-like types.
- Expose the public API through `pub mod prelude` with `pub use` of the types and traits a
  typical user needs, including external types that appear in public signatures. Keep
  everything else reachable by explicit path.

## 6. Error handling and panics

- Fallible operations return `Result<T, E>` with a crate-defined error. Use one `enum` per
  operation family; wrap an inner cause as a variant (`InvalidHeader(HeaderError)`).
- Error enums derive `Debug, Clone, PartialEq, Eq` and implement `Display` and
  `std::error::Error` (with `source()` for wrapped causes). Choose hand-written impls or
  `thiserror` once per crate and stay consistent.
- Propagate with `?`. Do not write `match x { Ok(v) => v, Err(e) => return Err(e) }`.
  Implement `From<Inner> for Outer` so `?` converts between error types.
- No `unwrap()` or `expect()` on values that depend on input or external data.
  `expect("reason")` is allowed only when the invariant is local and provable (for
  example converting a small integer literal to a generic numeric type); the message states
  the invariant.
- Library code must not panic on valid input. Every remaining panic condition is listed in a
  `# Panics` doc section and its message names the violated precondition. Prefer
  returning `Result`, or types that make the invalid state unrepresentable, to asserting.
- `assert!` / `assert_eq!` in non-test code only guard programmer errors. Use
  `debug_assert!` for expensive invariant checks in hot paths.
- `todo!()` and `unimplemented!()` never reach a release or `main`.
- `unsafe` is allowed only when no safe form exists or the safe form was measured and is
  materially slower (§10). Every `unsafe` block has a `// SAFETY:` comment stating the
  invariants that make it sound, and is as small as possible. Never create uninitialised
  values of ordinary types (use `MaybeUninit` or initialise with a default), and never
  extend a lifetime through a raw-pointer cast to satisfy the borrow checker; restructure
  the code instead. Run Miri in CI for crates that contain `unsafe`.

## 7. Documentation

- Enable `missing_docs`. Every public item has a `///` comment: one summary sentence,
  then details. Add `# Examples` (as doctests) to every public type and principal
  function, and `# Errors`, `# Panics`, `# Safety` sections where applicable.
- The crate root has `//!` docs with the purpose, a minimal complete example, a map of the
  main types and the operator table.
- Document what is not obvious from the signature: invariants, units, orientation and
  ordering conventions, numeric tolerances, complexity. Do not write `/// Returns x` for
  `fn x()`.
- Comments explain *why*. Do not leave commented-out code or timing logs in source;
  benchmark numbers belong in `benches/`. A one-line `// Name` comment (with a space
  after `//`) may head a group of impl blocks inside a file.
- Keep the README user-facing and in sync with the API. Compile its examples, either by
  including it as crate docs (`#![doc = include_str!("../README.md")]`) or by keeping the
  snippets in `examples/`.
- Maintain `CHANGELOG.md` (Keep a Changelog format). Record every user-visible change
  under its `Unreleased` section as part of the change that introduces it.

## 8. Testing

- Unit tests live in the file they test, in `#[cfg(test)] mod tests { use super::*; … }`
  at the end of the file. Each operation file tests its own operation. Integration tests
  of the public API go in `tests/`; documented examples are doctests.
- One behaviour per test, named for the behaviour. Several assertions are fine when they
  exercise the same behaviour on several inputs; prefer a table of `(input, expected)`
  pairs over copy-pasted blocks.
- Shared helpers go in `#[cfg(test)] pub(crate) mod test_util` at the crate root and are
  imported as `use crate::test_util::…`. Do not chain `pub use super::super::…` to share
  helpers between sibling test modules.
- Build fixtures with public constructors, or with `FromStr` when a parser exists and is
  itself tested. Keep fixture strings short and readable.
- Never compare computed floating-point results with `assert_eq!`. Use one shared
  approximate-equality helper that takes an explicit tolerance, and keep the tolerance
  tight (guideline: `1e-6` for `f32`, `1e-12` for `f64`). A test that needs a looser
  tolerance documents why in a comment, rather than widening silently.
- Test panics with `#[should_panic(expected = "...")]`. Do not write custom
  `catch_unwind` wrappers.
- Test errors by asserting the exact variant:
  `assert_eq!(parse("x").unwrap_err(), ParseError::InvalidDigit)`.
- Tests are warning-free: no unused imports, no `println!` or `dbg!`, no timing code, no
  `#[ignore]` without a reason string.
- Every bug fix adds a regression test. Every public function has at least one unit test
  or doctest. Consider property-based tests for algebraic laws (identities, inverses,
  round trips).
- `cargo test` passes locally before code is handed over for review. CI runs
  `cargo test --all-targets --all-features` and, when features exist,
  `--no-default-features`.

## 9. Dependencies

- Keep dependencies minimal and justify each new one in writing when it is introduced.
  Prefer the standard library: `std::sync::LazyLock` / `OnceLock` instead of
  `lazy_static` or `once_cell`.
- Specify versions as `"1.4.0"` (caret semantics with explicit minor and patch). No
  wildcards; released versions depend only on published crates. Disable default features
  you do not use (`default-features = false`).
- Keep `Cargo.lock` in the project, for libraries as well as binaries. Update
  dependencies (`cargo update`) as a separate change from code changes.
- Put heavy optional integrations (interop with a large external crate) behind a Cargo
  feature named after that crate.
- Test-only and benchmark-only crates go in `[dev-dependencies]`.
- `Cargo.toml` carries `description`, `license`, `readme`, `keywords`, `categories`,
  `edition` and `rust-version`. The project has a `LICENSE` file.
- Version the crate with semver. Before 1.0, breaking changes bump the minor version.
- Run `cargo audit` (or `cargo deny`) in CI once the crate is published.

## 10. Performance-sensitive code

- Measure first. Benchmarks live in `benches/` (criterion); report results alongside the
  change they justify, never in source comments.
- In hot loops: no allocation (preallocate with `with_capacity`, reuse buffers), no
  `clone()` of `Copy` values, no intermediate `collect()`, slices instead of `Vec`
  parameters, constants converted once outside the loop.
- Keep small plain-data types `Copy` and operate on them by value.
- Bit-level tricks are welcome when the invariants are documented next to the code with a
  worked example.
- Add `#[inline]` only to small, non-generic public functions where a benchmark showed
  cross-crate inlining matters. Generic functions are already inlined across crates.
- `unsafe` for speed requires a benchmark showing a material win (rule of thumb: at least
  10 percent on a realistic workload) plus the `// SAFETY:` comment from §6. Otherwise use
  the safe form (`vec![T::default(); n]`).
- Document numerical algorithms with their accuracy bounds and the input range over which
  they were validated.

## 11. Checklist for AI agents and reviewers

Before writing code:
- Read this file and the nearest sibling file that does the same kind of thing; copy its
  structure exactly.
- Decide where the code lives (§3): one operation per file, declared with `mod` in the
  parent, tests at the bottom of the same file.

While writing:
- Public item: doc comment with example, at least one test, correct verb form (§4).
- Fallible path: `Result` with an error variant; no new `unwrap()` outside tests.
- Operator: four ownership combinations delegating to one primitive, matching `*Assign`,
  no `unimplemented!()` stubs.
- `From`, not `Into`. No `Deref` inheritance. Fields private; `pub(crate)` for internals.
- Floats in tests use the approximate-equality helper.
- No commented-out code, timing notes, `println!`, or crate-level `#![allow]`.
- Every `unsafe` block has a `// SAFETY:` comment and a measured justification.

Before handing code over for review:
- `cargo fmt`
- `cargo clippy --all-targets --all-features -- -D warnings`
- `cargo test --all-targets --all-features`
- `RUSTDOCFLAGS="-D warnings" cargo doc --no-deps`
- Formatting-only changes kept separate from behavioural changes.

Reviewers additionally check: naming and spelling of new public items, growth of the
public API surface, documented panics, SAFETY comments, justification for new
dependencies, and the changelog entry.
