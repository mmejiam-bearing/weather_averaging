# AGENTS.md — Weather Daily Averaging Pipeline

## Project Overview

This repository implements a **physically rigorous daily energy-averaging pipeline** for global
ocean/atmosphere gridded weather data stored as `.npy` files on S3 (GFS-Wave / WAVEWATCH III,
0.25° resolution, available at 00/03/06/09/12/15/18/21 UTC). The primary consumer of most of
these averaged fields is a maritime routing engine that uses them as heuristic inputs for ship
speed-loss estimation and sailability graph edge weights; `wad`/`pwd` specifically are produced
for a separate downstream consumer in another repository (see README.md Section 4.7), not for
this repo's own routing heuristic.

The archive uses two different grids (confirmed against real downloaded data, see README.md
Section 2): wave fields (`sigwh`/`wsh`/`wsp`/`wsd`/`pwh`/`pwd`) are `621×1440`, wind fields
(`was`/`wad`) are `721×1440`. On disk, continuous fields are `float16` and direction fields
(`wsd`/`wad`/`pwd`) are `int8` 16-point-compass codes — **not** `float64` degrees. `src/io.cpp`
upcasts every field to `double` immediately on load and `src/averaging.cpp`'s
`compass_code_to_degrees` decodes the compass codes before any averaging; nothing downstream of
the I/O layer should ever see anything but `double` degree values.

The pipeline is implemented in **C++23** using **Eigen 3** for matrix operations and **OpenMP**
for shared-memory parallelism. cnpy is used for `.npy` I/O.

---

## Mandatory Conventions

### Language and Standard
- All source files: **C++23**. Compiler must be invoked with `-std=c++23`.
- No exceptions to this rule. No C++17 fallbacks, no `#ifdef __cplusplus` guards for older standards.

### Memory and Resource Management

#### Rule of Zero — Gold Standard

The **Rule of Zero** is the primary target for every class and struct in this codebase: if a
type does not directly manage a resource (raw memory, file handle, socket, lock), it should
declare **none** of the five special member functions (destructor, copy constructor, copy
assignment, move constructor, move assignment). The compiler-generated defaults are correct and
optimal when all members are themselves resource-safe types.

Concretely:

- Prefer standard library containers (`std::vector`, `std::string`, `std::array`) and smart
  pointers over raw owning pointers. These types already implement the Rule of Five correctly;
  composing them means your type needs no special members.
- Eigen matrix types (`Eigen::MatrixXd`) own their heap allocation internally and are fully
  RAII-compliant — use them as plain value members, no manual management required.
- cnpy `NpyArray` objects manage their own buffer; treat them as value types, load at scope
  entry, and let them destruct at scope exit. Do not cache the raw data pointer past the
  `NpyArray` lifetime.
- Never use `new`, `delete`, `malloc`, `free`, `realloc`, or `calloc` directly. If a third-party
  API returns a raw owning pointer, wrap it immediately in a `std::unique_ptr` with an
  appropriate deleter before doing anything else with it.

#### Smart Pointer Preference: `unique_ptr` First

`std::unique_ptr` is the **default smart pointer** for all heap-allocated resources in this
codebase. `std::shared_ptr` is a last resort, used only when shared ownership is genuinely
required and cannot be redesigned away.

**Use `std::unique_ptr` when:**
- A single owner is responsible for the lifetime of the resource — which covers the vast
  majority of cases in this pipeline (loaded grids, intermediate buffers, I/O handles).
- Transferring ownership between scopes or into a container: `std::move` makes the transfer
  explicit and auditable.
- Wrapping a raw owning pointer returned by a C API:
  ```cpp
  // Correct: wrap immediately, never touch the raw pointer again
  auto handle = std::unique_ptr<CHandle, decltype(&c_api_free)>(
      c_api_alloc(), &c_api_free
  );
  ```

**Use `std::shared_ptr` only when:**
- The resource has genuinely shared, non-deterministic lifetime — i.e., multiple independent
  owners exist and it is architecturally impossible to designate one as primary owner.
- Even then, construct exclusively via `std::make_shared<T>(...)`, never via
  `std::shared_ptr<T>(new T(...))`. The latter performs two separate heap allocations (object +
  control block); `make_shared` fuses them into one, which is both faster and safer (no
  exception-safety gap between the `new` and the constructor of `shared_ptr`).
- When exposing a shared resource to observers that must not extend its lifetime, give them
  `std::weak_ptr`, not a copy of the `std::shared_ptr`. Observers must call `.lock()` before
  use and handle the null case.

**`shared_ptr` is explicitly banned when:**
- A `unique_ptr` would suffice. If you find yourself reaching for `shared_ptr` for convenience
  (e.g., to avoid passing by reference or to sidestep a design question about ownership), that
  is a design smell — resolve the ownership question instead.
- As a function parameter type. Pass the underlying object by reference (`const T&` or `T&`)
  or by raw non-owning pointer (`T*`) to callees that do not participate in ownership. Passing
  `shared_ptr` by value to a function that does not store it is an unnecessary reference-count
  increment/decrement pair on every call.

**Summary:**

```
Does this object need heap allocation?
│
├─ No  →  Stack / value member. No smart pointer needed.
│
└─ Yes →  Does it have a single, clearly identifiable owner?
          │
          ├─ Yes →  std::unique_ptr<T>  (default choice)
          │
          └─ No  →  Is shared ownership genuinely unavoidable?
                    │
                    ├─ No  →  Redesign to establish single ownership.
                    │
                    └─ Yes →  std::make_shared<T>(...)
                               Observers get std::weak_ptr<T>.
```

#### Rule of Five — Fallback When Rule of Zero is Impossible

When a class **must** directly manage a resource (e.g., a custom RAII wrapper around a C API
that has no suitable standard-library analogue), apply the **Rule of Five**: if any one of the
five special member functions needs a non-default implementation, **all five must be explicitly
defined or explicitly defaulted**. No exceptions.

The five members are:

```cpp
class ResourceHolder {
public:
    ~ResourceHolder();                                        // destructor
    ResourceHolder(const ResourceHolder&);                    // copy constructor
    ResourceHolder& operator=(const ResourceHolder&);         // copy assignment
    ResourceHolder(ResourceHolder&&) noexcept;                // move constructor
    ResourceHolder& operator=(ResourceHolder&&) noexcept;     // move assignment
};
```

Rules for Rule-of-Five types:

- If the resource is not copyable (e.g., an exclusive file handle), **delete** the copy
  constructor and copy assignment operator explicitly (`= delete`). Do not leave them undefined
  (undefined triggers a linker error, which is worse than a clear compile-time diagnostic).
- Move constructor and move assignment must be marked `noexcept` whenever the move cannot throw.
  This is required for standard library containers to use the move path instead of the copy path
  (via `std::move_if_noexcept`), which is critical for performance when storing these objects in
  `std::vector` or passing them across thread boundaries.
- If a Rule-of-Five type is needed, add a comment above the class declaration explaining *why*
  Rule of Zero could not be applied, naming the specific resource and the reason no standard
  wrapper exists for it.

#### Summary Decision Tree

```
Does the type directly own a raw resource?
│
├─ No  →  Rule of Zero: declare no special members. Done.
│
└─ Yes →  Can the ownership be delegated to a standard-library type
          (unique_ptr, vector, string, etc.)?
          │
          ├─ Yes →  Refactor to use that type. Apply Rule of Zero. Done.
          │
          └─ No  →  Rule of Five: define or default/delete all five.
                     Add a comment explaining why Rule of Zero was impossible.
```

### Naming Conventions
- Types: `PascalCase` (e.g., `WindSeaResult`, `EnergyStack`)
- Functions: `snake_case` (e.g., `energy_average`, `swell_residual`)
- Constants: `SCREAMING_SNAKE_CASE` (e.g., `GRID_ROWS`, `GRID_COLS`)
- Member variables: `snake_case` with no prefix (no `m_`, no `_` suffix)
- Files: `snake_case.cpp` / `snake_case.hpp`

### OpenMP
- Parallelism strategy: **outer loop over days** using `#pragma omp parallel for schedule(dynamic)`
  at the day-dispatch level; inner cell loops (`collapse(2)` over rows×cols) are used only for
  computations that are too coarse-grained at the day level.
- Never use `omp_set_nested(1)` — nested parallelism is disabled. One level of parallelism only.
- All accumulator variables inside parallel regions must be declared with `reduction` clauses or
  thread-local copies. No shared mutable state without explicit synchronization.

### Eigen Usage
- Use `Eigen::MatrixXd` (dynamic, double precision) for every field in memory, regardless of its
  on-disk dtype or grid (wave-grid fields are `621×1440`, wind-grid fields are `721×1440` — see
  `config::WAVE_GRID_ROWS` / `WIND_GRID_ROWS`). Never hard-code a single fixed row count.
- Use `Eigen::RowMajor` storage when mapping from cnpy (cnpy loads row-major C-order arrays).
- Never use fixed-size Eigen matrices for grid-sized objects (triggers stack overflow).
- Vectorization is handled by Eigen's internal SIMD; do not add manual SIMD intrinsics.

### NaN Handling
- Land cells and missing data are represented as `NaN` (`std::numeric_limits<double>::quiet_NaN()`).
- All averaging functions must use NaN-aware accumulation: skip cells where `std::isnan(v)` is true.
- Output cells where all timesteps are NaN must propagate NaN — never output zero for missing data.
- The swell residual must be clipped to zero (not NaN) when the arithmetic produces a small
  negative value due to floating-point noise: `std::max(sigwh² - wsh², 0.0)`.

### Mathematical Derivations
See the README.md for full derivations. This section summarizes constraints that must be
respected in code:
- Energy-averaged height: `H_eff = sqrt(mean(H²))` over valid (non-NaN) timesteps only.
- Energy-weighted period: `T_eff = sum(H² * T) / sum(H²)` — denominator is the sum of weights,
  not the count of timesteps.
- Energy-weighted direction: circular mean weighted by `H²`, computed via unit-vector components
  `u = sin(θ)`, `v = cos(θ)`, result via `atan2(sum(H²·u), sum(H²·v)) mod 360`.
- `atan2` usage: `std::atan2(sum_u, sum_v)` — note argument order is `(sin, cos)` not `(y, x)`;
  this implements the meteorological clockwise-from-North convention correctly. Do not swap.

### Testing
- Every averaging function must have a corresponding unit test in `tests/`.
- Tests are implemented with **GoogleTest** (fetched via CMake `FetchContent`, pinned at `v1.14.0`).
- Mandatory test cases:
  - Scalar energy average: known input → known output, verified analytically.
  - Circular mean: the (359°, 1°) → ~0° case must be a test. The (90°, 270°) → undefined/NaN
    case (cancellation) must be a test.
  - NaN propagation: all-NaN input → NaN output; one-valid-timestep input → valid output.
  - Swell residual clamp: `wsh > sigwh` by epsilon → result is 0.0, not NaN.
  - Period weighting: zero-height timestep → contributes zero weight to T_eff (not equal weight).

### No Implicit Conversions
- All casts between numeric types must be explicit (`static_cast<double>(...)`).
- Do not mix `float` and `double`. All fields are `double`.

---

## Directory Structure

```
weather-daily-avg/
├── AGENTS.md               ← this file
├── README.md               ← full mathematical documentation
├── CMakeLists.txt          ← top-level build
├── cmake/
│   └── FindCnpy.cmake      ← find-module for cnpy if not using add_subdirectory
├── src/
│   ├── main.cpp            ← CLI entry point, date-range iteration
│   ├── averaging.hpp       ← declarations: energy_average, daily_energy_average, swell_residual
│   ├── averaging.cpp       ← implementations
│   ├── io.hpp              ← declarations: load_npy, save_npy, path_for
│   ├── io.cpp              ← implementations
│   ├── types.hpp           ← shared types: WindSeaResult, EnergyStack, GridShape
│   └── config.hpp          ← compile-time constants: GRID_ROWS, GRID_COLS, HOURS
├── tests/
│   ├── CMakeLists.txt
│   ├── test_averaging.cpp
│   └── test_io.cpp
├── cnpy/                   ← cnpy as git submodule
│   ├── cnpy.h
│   └── cnpy.cpp
└── scripts/
    └── sync_s3.sh          ← S3 sync script (see README for usage)
```

---

## What Agents / Contributors Must NOT Do

- Do not change the energy-averaging math to arithmetic mean without explicit sign-off and a
  documented physical justification. The H² weighting is load-bearing for the routing heuristic.
- Do not add `float` versions of any function to "save memory" — even though the S3 archive
  stores fields as `float16`/`int8` on disk, `src/io.cpp` upcasts everything to `double`
  immediately on load, and precision loss in the in-memory averaging math propagates into
  routing decisions.
- Do not use `std::thread` directly — all parallelism goes through OpenMP to respect the
  compile-time thread-count configuration.
- Do not add external dependencies beyond Eigen, OpenMP, cnpy, and zlib without updating both
  CMakeLists.txt and this file.
- Do not remove NaN guards. Replacing `std::isnan` with `v == v` (NaN-check trick) is acceptable
  but must be documented inline.
- Do not write output files in place of input files. Input paths and output paths are always
  distinct; the output directory is set via `--output-dir` CLI argument.

---

## CLI Contract

The compiled binary must accept the following arguments:

```
weather_daily_avg \
  --data-dir   <path>          # root of local S3 mirror, e.g. ./data
  --output-dir <path>          # root for daily-averaged output, e.g. ./averaged
  --years      <y1> [y2 ...]   # one or more years, e.g. 2024 2025
  --months     <m1> [m2 ...]   # one or more months (1–12), e.g. 1 2 3
  [--threads   <n>]            # optional: override OMP_NUM_THREADS
```

Output files are written to:
```
<output-dir>/<year>/<month>/<day>/
    sigwh.npy
    wsh.npy
    wsp.npy
    wsd.npy
    was.npy
    wad.npy
    pwd.npy
    swell_residual.npy
```

`wad`/`pwd` are produced for a downstream consumer in another repository, not for this
pipeline's own routing heuristic (see README.md Section 4.7). `pwh` is read as an input (to
weight `pwd`'s energy-weighted circular mean) but is never itself written as an output.

Days with fewer than 2 valid hourly timesteps (out of 8 possible) are skipped and a warning is
emitted to `stderr`. Days with 0 valid timesteps produce no output file (not an error).
