# DeepNAV

DeepNAV is an orbit-determination and flight-dynamics sandbox for DSN/Voyager
tracking work.
The current implementation is centered on C++/Eigen numerical utilities,
CSPICE-backed station geometry, synthetic observations, perturbation models,
batch weighted least-squares filtering, and an initial optical-navigation
geometry and camera-model layer.

## Current Scope

Implemented modules:

- CSPICE-backed DSN station catalog.
- Adaptive RKF45 and DP853 numerical integration with reusable ephemeris
  interpolation.
- Batch Weighted Least Squares filtering with a priori information and an
  iterative convergence driver.
- Synthetic one-way range, range-rate, and VLBI differential range observations.
- Interchangeable propagated-ephemeris and direct-CSPICE target-state sources.
- Solar Shapiro delay applied to synthetic radiometric observables.
- Modular perturbation models for third-body gravity and cannonball SRP.
- A reusable Cartesian state type and simultaneous-epoch geometric
  camera-to-target line of sight, independent of CSPICE.
- Pinhole camera projection with optional OpenCV-compatible radial and
  tangential distortion.
- Two-state residual clock dynamics, seeded truth simulation, Allan-data
  configuration, covariance propagation, and conservative calibration budgets.
- Focused tests for numerical integration, ephemeris interpolation, synthetic
  observations, optical geometry and projection, station-kernel sanity checks,
  and a Voyager 1 position OD case that exercises WLS and the perturbation
  models.

Planned or partial areas:

- Promotion of the estimator's computed-observation and dynamics code from the
  Voyager test into reusable library modules.
- Full dynamics builder that sums central gravity, third bodies, SRP, and future force models.
- Sequential filters and richer OD diagnostics.
- Optical-navigation attitude transformations, light-time/apparent-direction
  corrections, image measurements, and visualization work.

Longer-term directions include nonlinear optimization for trajectory and
maneuver design, plus multi-agent spaceborne/ground antenna tracking. CppAD,
IFOPT, and IPOPT are not yet wired into the current CMake build.

## Repository Layout

The active library and test modules wired into CMake are:

```text
include/
  Clocks/
    Types.hpp
    Clocks.hpp
    ClockHistory.hpp
    ClockModel.hpp
    LocalOscillator.hpp
    DSAC.hpp
    ClockTruthSimulator.hpp
    Allan.hpp
    Calibration.hpp
  DP853Integrator.hpp
  RKF45Integrator.hpp
  dynamics/
    CartesianState.hpp
    EphemerisInterpolator.hpp
    SpiceInterpolator.hpp
  opnav/
    CameraModel.hpp
    DistortionModel.hpp
    Types.hpp
    core/
      GeometricLineOfSight.hpp
  stations/
    ElevationMask.hpp
    StationCatalog.hpp
    Stations.hpp
  filters/
    BatchLeastSquaresDriver.hpp
    filter.hpp
    WLS.hpp
  observations/synth/
    obs_synth.hpp
    RangeSynth.hpp
    RangeRateSynth.hpp
    TargetStateProvider.hpp
    VLBISynth.hpp
  perturbations/
    Gravitational.hpp
    SRP.hpp
    Shapiro.hpp
  utils/
    ClockCsvWriter.hpp
    CSPICE/
      SpiceError.hpp
      SpiceErrorModeGuard.hpp

src/
  Clocks/
    Clocks.cpp
    ClockHistory.cpp
    ClockModel.cpp
    LocalOscillator.cpp
    DSAC.cpp
    ClockTruthSimulator.cpp
    Allan.cpp
    Calibration.cpp
    Validation.hpp
  utils/
    ClockCsvWriter.cpp
  RKF45Integrator.cpp
  dynamics/
    EphemerisInterpolator.cpp
    SpiceInterpolator.cpp
  opnav/
    CameraModel.cpp
    DistortionModel.cpp
    core/
      GeometricLineOfSight.cpp
  stations/
    ElevationMask.cpp
    StationCatalog.cpp
  filters/
    BatchLeastSquaresDriver.cpp
    filter.cpp
    WLS.cpp
  observations/synth/
    obs_synth.cpp
    RangeSynth.cpp
    RangeRateSynth.cpp
    TargetStateProvider.cpp
    VLBISynth.cpp
  perturbations/
    Gravitational.cpp
    SRP.cpp

tests/
  test_clocks.cpp
  test_clock_history.cpp
  test_plot_clocks.py
  test_ephemeris_interpolator.cpp
  test_geometric_line_of_sight.cpp
  test_optical.cpp
  test_synth_observations.cpp
  test_synth_source_comparison.cpp
  test_voyager_position_od.cpp

test_rkf45.cpp
test_stations.cpp
station_catalog_demo.cpp
kernels.tm
plot_clocks.py
```

## Dependencies

The CMake project is named `DeepNAV`, uses ISO C++23, and links these dependencies:

- Eigen3
- CSPICE
- GLFW, GLEW, OpenGL
- local `third_party/imgui`, `third_party/implot`, and `third_party/glm`

On macOS, the checked-in preset uses `/usr/bin/clang` and `/usr/bin/clang++`.
This keeps the compiler and macOS SDK matched and avoids the missing
`INFINITY`/`NAN` definitions observed with Homebrew LLVM 21 in strict C++23
mode. ImGui and ImPlot include paths are marked as third-party `SYSTEM` paths;
their intentional raw-memory implementation warning is disabled only for the
`imgui_lib` target, not for DeepNAV sources.

The optional diagnostic scripts under `tests/` use Python 3, NumPy, Matplotlib,
and, for `plot_observability.py`, pandas.

`CSPICE_HOME` must point to the CSPICE installation:

```sh
export CSPICE_HOME=/path/to/cspice
```

The project expects the Voyager/DSN kernels under `Kernels/` and the meta-kernel
at `kernels.tm`. The Voyager Jupiter-encounter source-comparison and OD tests
require `jup310.bsp` for 1979 coverage of Jupiter body `599` and Galilean moons
`501-504`.

Important path convention: `kernels.tm` uses `PATH_VALUES = ( '../Kernels' )`.
For tests or demos that load `../kernels.tm`, run from `build-clang`.

## Build

On Apple Silicon, the included Apple Clang preset uses Homebrew dependency
paths and creates the conventional `build-clang` directory:

```sh
cmake --preset apple-clang
```

When changing compiler presets in an existing build directory, refresh the
CMake cache before rebuilding:

```sh
cmake --fresh --preset apple-clang
```

For other toolchains, configure the same directory directly:

```sh
cmake -S . -B build-clang -DCMAKE_BUILD_TYPE=Release
```

```sh
cmake --build build-clang --target test_rkf45 -j4
cmake --build build-clang --target test_ephemeris_interpolator -j4
cmake --build build-clang --target test_synth_observations -j4
cmake --build build-clang --target test_synth_source_comparison -j4
cmake --build build-clang --target test_voyager_position_od -j4
cmake --build build-clang --target test_optical -j4
cmake --build build-clang --target test_geometric_line_of_sight -j4
cmake --build build-clang --target station_catalog_demo -j4
```

The current `CMakeLists.txt` filters source candidates with `EXISTS` so missing
legacy examples do not break generation.

## Tests

Run focused tests:

```sh
ctest --test-dir build-clang -R test_rkf45 --output-on-failure
ctest --test-dir build-clang -R test_ephemeris_interpolator --output-on-failure
ctest --test-dir build-clang -R test_synth_observations --output-on-failure
ctest --test-dir build-clang -R test_synth_source_comparison --output-on-failure
ctest --test-dir build-clang -R test_voyager_position_od --output-on-failure
ctest --test-dir build-clang -R test_optical --output-on-failure
ctest --test-dir build-clang -R test_geometric_line_of_sight --output-on-failure
ctest --test-dir build-clang -R '^test_clocks$' --output-on-failure
ctest --test-dir build-clang -R '^test_clock_history$' --output-on-failure
python3 test_clocks.py # Independent reference checks; requires NumPy.
python3 -m unittest discover -s tests -p test_plot_clocks.py # Requires Pandas/Matplotlib.
```

SPICE-backed tests are registered with `build-clang` as their working directory
so `../kernels.tm` resolves correctly.

Run `test_synth_observations` before `test_voyager_position_od` when the
synthetic report needs to be generated or refreshed; the OD test consumes that
report.

`test_stations.cpp` is an older kernel smoke test and currently hard-codes
`../Kernels.tm`. On case-sensitive filesystems, either update that literal to
`../kernels.tm` or provide a matching compatibility file before relying on it in
an all-test run.

## Implemented Components

### Clocks

The independent `include/Clocks` module uses the `fd::clocks` namespace.
`Clocks/Clocks.hpp` defines the abstract interface; `LocalOscillator`, `DSAC`,
`ClockModel`, and `ClockTruthSimulator` each have their own header and source
file. Common types, Allan analysis, and calibration helpers are separated into
`Types.hpp`, `Allan.hpp`, and `Calibration.hpp`. Bias is clock time minus reference time in seconds;
fractional frequency is dimensionless. Use elapsed time since calibration.
Known deterministic frequency drift has units s^-1. The model contains white
frequency diffusion (`q_bias_s`, units s) and random walk frequency diffusion
(`q_frequency_per_s`, units s^-1).

Both clock classes expose `fromParameters` and `fromAllanData`. Allan fitting
uses variance, requires explicit drift removal (or zero drift), a declared
valid tau interval and a noise assumption, and rejects negative coefficients,
ill-conditioned data and excessive residuals. One Allan point is accepted only
with a single dominant noise assumption. The default maximum relative variance
residual is 5%, configurable by the caller. All DSAC configurations are custom;
no published hardware preset is supplied. Allan data do not determine initial
bias, fractional frequency, their covariance, or the removed drift.

```cpp
#include "Clocks/Calibration.hpp"
#include "Clocks/ClockTruthSimulator.hpp"
#include "Clocks/LocalOscillator.hpp"
using namespace fd::clocks;

// Synthetic coefficients for illustration, not hardware specifications.
const auto clock = LocalOscillator::fromParameters({2e-15, 1e-24, 3e-26});
ClockTruthSimulator truth(clock, 42); // clock must outlive truth.
const ClockState initial{2e-9, 1e-11};
const auto predicted = clock.propagate(initial, 10.0);
const auto simulated = truth.step(initial, 10.0);
const auto covariance = propagateClockCovariance(
    clock, ClockCovariance::Zero(), 10.0);

const auto interval = clockCalibrationInterval(1.0, 30.0 * 86400.0,
    clock.parameters(), ClockBudget{.initial_state = initial});
// nullopt means > the specified horizon; zero means initial threshold contact.
```

`ClockTruthSimulator` uses exact correlated increments, including zero and
single-noise cases. Seeds reproduce samples within the same standard-library
implementation; they do not reproduce Python samples. Deterministic propagation
and covariance calculations draw no random numbers. Covariance propagation
supports a general symmetric positive semidefinite initial covariance.

`overlappingAllanDeviation` consumes uniformly sampled bias and an integer
averaging factor; `theoreticalAllanDeviation` optionally includes unremoved
deterministic drift. The calibration helper instead uses a monotonic envelope
with absolute deterministic terms and independent initial uncertainties. Its
default k=3 is a pointwise budget, not a probability guarantee for an entire
trajectory. Thresholds bound the equivalent one-way range contribution in
metres, not 3D position error or operational contact schedules. This two-noise
model excludes environmental effects, flicker noise and uncertain drift.

Build with `cmake --build build-clang --target test_clocks -j4`, then run the
focused CTest command above. The test target links Eigen only and exercises
both implementations through the base interface, Allan conversion/rejection,
empirical Q including cross-covariance, pure and combined Allan noise, drift,
general covariance propagation, and analytic calibration limits. The module is
also included in `big_functions`; it is not yet connected to an EKF or to the
radiometric measurement generators.

Clock history is optional and owned by the base `Clocks` class through
`ClockHistory`. Clocks are movable and noncopyable, so a large trajectory is
never copied implicitly. `enableHistory(capacity, includeCovariance)` creates
a fresh empty collection and reserves one contiguous sample vector; recording
cannot exceed the explicit sample limit or grow the allocation. Five doubles
per sample use approximately 40 bytes of payload (about 40 MB for one million
samples). `clearHistory()` retains the allocation; `disableHistory()` releases
it. Manual `recordSample()` calls are ignored when history is disabled.

`propagate()` and `step()` do not record. `ClockTruthSimulator::run()` explicitly
records accepted truth states, starting at elapsed time zero, and requires an
empty history. Pass a mutable clock reference to the simulator for recording.
The stride records the initial sample, every `recordEvery` steps, and the final
sample once. `ClockHistory::requiredCapacity(N, r)` computes the required limit
with overflow checks. The variable-step overload accepts a span of time steps
and uses the actual accumulated elapsed time. Each Monte Carlo realization
must own its clock/history; there is no shared-writer synchronization.

```cpp
#include "Clocks/ClockTruthSimulator.hpp"
#include "Clocks/LocalOscillator.hpp"
#include "utils/ClockCsvWriter.hpp"
using namespace fd::clocks;

auto clock = LocalOscillator::fromParameters({2e-15, 1e-24, 3e-26});
constexpr std::size_t steps = 360, stride = 6;
clock.enableHistory(ClockHistory::requiredCapacity(steps, stride), true);
ClockTruthSimulator truth(clock, 42);
const ClockCovariance p0 = ClockCovariance::Zero(); // Explicit known initial uncertainty.
const auto finalState = truth.run(steps, 240.0, {2e-9, 1e-11}, stride, p0);
fd::utils::writeClockCsv("local.csv", *clock.history()); // Refuses existing files.
// To replace a file, explicitly pass fd::utils::ClockCsvWriteMode::Overwrite.
```

Covariance history requires an explicit P0 and stores only `sqrt(P_bb)` and
`sqrt(P_yy)` alongside time, bias and fractional frequency. It does not infer
uncertainty from truth samples. CSV export occurs after simulation, preserves
double precision and units, and omits both sigma columns when covariance is
absent. An empty history exports only the header. `writeClockCsv` reports file
opening/writing errors and never appends implicitly.

`plot_clocks.py` is standalone and reads one or more CSV files into separate
Pandas DataFrames. It rejects invalid schemas, duplicate headers, missing or
nonfinite numbers, negative sigma and nonincreasing/negative timestamps. It
creates three separate figures/windows, with an independent vertical scale
for every clock so the smaller DSAC curve remains visible:

- Main figure: signed one-way range error and total pointwise budget
  `c*(abs(mu_bias) + k*sigma_bias_s)` (default k=3).
- Uncertainty figure: stochastic contribution `k*c*sigma_bias_s`, excluding
  the deterministic mean.
- Support figure: bias in ns and fractional frequency.

The total budget requires explicit deterministic mean histories supplied via
`--mean-csv`, one per input clock, at identical recorded timestamps. Those CSVs
use the same primitive schema and store the expected bias obtained by noiseless
propagation; the noisy realization is never treated as its own mean. Without
mean histories or covariance, the corresponding panel explains the missing
inputs instead of inventing a budget. This pointwise budget is not Allan
deviation or a guarantee for an entire sample path.

`--threshold-m` draws signed thresholds for range and positive thresholds for
budget/uncertainty, annotating the first recorded contact from all input rows
before display reduction. It reports the preceding sample time too. These are
sampled contacts; an earlier unobserved crossing between samples can be missed.

```sh
python3 -m pip install pandas matplotlib
python3 plot_clocks.py local.csv dsac.csv --mean-csv local_mean.csv dsac_mean.csv --labels Local DSAC --threshold-m 1 --output clocks.png
python3 plot_clocks.py local.csv --time-unit days --threshold-m 1 --output clocks.pdf
python3 plot_clocks.py local.csv --max-points 0 --show --output clocks.svg
```

The plotter uses a noninteractive backend unless `--show` is supplied. PNG,
PDF and SVG are supported. `--output clocks.png` saves the main figure there
and the other figures to `clocks_uncertainty.png` and `clocks_diagnostics.png`.
With `--show`, all three windows open together. Its default limit of 20000 points per curve preserves
the endpoints and only affects display; narrow peaks can be missed. Original
CSV data are untouched. Reduced histories (especially an extra final sample)
may be nonuniform and must not automatically be used for Allan estimates or
first-crossing analysis. Keep simulation parameters and seeds in the run's
configuration; filenames do not establish hardware performance.

Build the recording checks with
`cmake --build build-clang --target test_clock_history -j4`. Running
`build-clang/test_clock_history` automatically creates `Output clocks` in the
project root and saves `local.csv`, `dsac.csv`, `local_mean.csv`, and
`dsac_mean.csv` there, replacing previous synthetic outputs. Both realizations
cover the same 24-hour interval and include covariance. Mean files are exported
separately from deterministic propagation. This destination is independent of the working directory;
change `defaultOutputDirectory` in `tests/test_clock_history.cpp` to customize it,
or pass an output directory as a command-line argument. For example:

```sh
./build-clang/test_clock_history
python3 plot_clocks.py "Output clocks/local.csv" "Output clocks/dsac.csv" \
  --mean-csv "Output clocks/local_mean.csv" "Output clocks/dsac_mean.csv" \
  --labels "Local (synthetic)" "DSAC (custom synthetic)" \
  --threshold-m 1 --output "Output clocks/clocks.png" --show
```

CSV validation artifacts remain temporary; only simulation CSVs are retained.
Generated files in `Output clocks` are ignored by Git. Their
parameters and seeds are defined in `tests/test_clock_history.cpp` and printed
at execution, including initial conditions, drift, noise intensities, P0 and
sampling. The coefficients remain synthetic validation values; the mathematical
Allan checks do not establish OCXO/DSAC hardware performance or ESA suitability.
Python tests
also exercise C++ exports when that executable is available.

### Optical Navigation

The optical-navigation code is split into a CSPICE-independent geometry layer
and a camera projection layer.

`fd::dynamics::CartesianState` stores position in kilometres and velocity in
kilometres per second. Its `allFinite()` validation covers both state
components. Velocity is retained for future apparent-direction and light-time
work, although the current simultaneous-epoch geometry uses position only.

`fd::opnav::core::computeGeometricLineOfSight(...)` calculates:

```text
vectorKm      = target.positionKm - camera.positionKm
rangeKm       = |vectorKm|
unitDirection = vectorKm / rangeKm
```

The function rejects non-finite states and coincident camera/target positions.
It does not perform CSPICE retrieval, light-time iteration, stellar aberration,
attitude transformations, or camera projection. `test_geometric_line_of_sight`
covers a hand-computable geometry, unit direction and range, translation
invariance, coincident positions, and non-finite input.

`fd::opnav::CameraModel` projects a direction already expressed in the camera
frame onto pixel coordinates. It supports an identity distortion model and an
OpenCV-compatible rational radial plus tangential distortion model. Directions
with non-finite components or a non-positive camera-frame Z component are
rejected. `test_optical` covers ideal projection, zero and non-zero distortion,
boresight projection, and targets behind the camera.

The current processing boundary is therefore:

```text
Cartesian states -> geometric LOS -> camera-frame direction -> pixel
```

The attitude transformation needed to produce the camera-frame direction is
not implemented yet.

### Station Catalog

`include/stations/StationCatalog.hpp` and `src/stations/StationCatalog.cpp` provide:

- `defaultDsnStationCatalog()`
- `stationNaifIdFromName(...)`
- `buildStationFromKernel(...)`
- `buildDefaultDsnCatalogFromKernel(...)`

Kernel loading remains the caller's responsibility. The catalog code queries
CSPICE and returns `od::Station` objects without mutating the base station model.

`include/stations/ElevationMask.hpp` and `src/stations/ElevationMask.cpp`
compute station-target elevation from CSPICE ITRF93 geometry. The current kernel
set does not include station topo-frame definitions such as `DSS-43_TOPO`, so
the helper derives the local up vector from station geodetic coordinates.

### Numerical Integration and Ephemeris Interpolation

`RKF45Integrator` uses:

- `Eigen::VectorXd` state representation.
- `Eigen::Ref` derivative callbacks.
- Fehlberg 4(5) embedded error estimate.
- Adaptive step acceptance/rejection.
- Preallocated trial, fourth-order, fifth-order, and stage matrices.
- Accepted-step history for downstream Hermite interpolation.

`od::EphemerisInterpolator` stores RKF45 state and derivative nodes and evaluates
cubic Hermite states with logarithmic interval lookup.

`od::DP853Integrator` is a header-only adaptive eighth-order integrator for
six-element states. It can retain accepted-step history in the same node format
used by `od::EphemerisInterpolator` and is used by the Voyager OD test.

`od::SpiceInterpolator` samples CSPICE states over a configured interval and
uses cubic Hermite interpolation for repeated body and station state queries.

`test_rkf45` propagates a simple circular Kepler orbit for one period and checks
position closure, velocity closure, energy, angular momentum, and runtime.
`test_ephemeris_interpolator` covers interpolation nodes, interior values, and
out-of-range error handling.

### WLS Filter

`fd::filters::WLS` implements a batch Weighted Least Squares update with a priori
information:

```text
Lambda = H^T R^-1 H + P0^-1
N      = H^T R^-1 r + P0^-1 (x_prior - x_nominal)
dx     = Lambda^-1 N
```

The implementation uses Eigen solves rather than explicitly forming `R^-1`.
LDLT is attempted first, with QR fallback for robustness.

`fd::filters::BatchLeastSquaresDriver` repeatedly linearizes and applies WLS
updates. It supports iteration callbacks and convergence checks based on
weighted RMS and whole-state, component, or state-block correction tolerances.

### Synthetic Observations

Synthetic observation classes live under `fd::observations::synth`.

Shared features:

- `GeometryConfig`: target, station, frame, aberration correction.
- `NoiseConfig`: enabled flag, sigma, random seed.
- Configurable time grid via `makeEpochGrid(start, end, step)`.
- CSPICE station names resolved through the DSN station catalog.
- A `TargetStateProvider` interface supporting either propagated/interpolated
  target states or direct geometric states from CSPICE.

`RangeSynth`:

- Gets light-time corrected station-to-target geometry through CSPICE.
- Adds solar Shapiro range delay.
- Adds optional Gaussian range noise.

`RangeRateSynth`:

- Computes one-way range-rate as centered differenced range over a configurable
  count time, defaulting to 60 s.
- Includes solar Shapiro delay in both endpoint ranges, so the differenced
  observable naturally includes the Shapiro rate.
- Adds optional Gaussian range-rate noise with sigma scaled by
  `1 / sqrt(count_time_s)`.

`VLBISynth`:

- Generates differential one-way range (`station2 - station1`) for a station pair.
- Keeps samples only when the target is visible from both stations.
- Applies the same one-way light-time and solar Shapiro treatment used by the range observable.
- Adds optional Gaussian VLBI delay noise in km.

`test_synth_observations` seeds the Voyager 1 state from CSPICE once, propagates
that state with RKF45 and the same dynamics used by the OD test, applies manual
one-way light-time and Shapiro delay against the propagated ephemeris, and
generates elevation-filtered station and VLBI observations:

- Stations: DSS-43, DSS-63
- Target: Voyager 1 (`-31`), initial state seeded at arc start
- Start: `1979-01-10T00:00:00`
- Duration: 4 days
- Elevation mask: 10 degrees
- Cadence: 3 minutes
- VLBI cadence: 20 seconds
- Samples after mask: 1655 total; DSS-43 659, DSS-63 996
- VLBI samples after dual-station visibility mask: 2902 total; DSS-43/DSS-63
  82, DSS-14/DSS-43 2820
- Range sigma: 0.010 km
- Range-rate base sigma: `1.0e-6 km/s` at 1 s
- Range-rate count time: 60 s
- Range-rate row sigma: `1.290994e-7 km/s`
- VLBI sigma: 0.001 km

The output report keeps the historical filename but now includes a `station`
column:

```text
synthetic_observations_dss43_voyager1.txt
```

`test_synth_source_comparison` evaluates range, range-rate, and VLBI with a
propagated target-state provider and reports differences from the direct
CSPICE source. This is a diagnostic comparison: it verifies that every
comparison set is populated but does not impose accuracy thresholds.

### Perturbations

`fd::perturbations::Shapiro` provides inline reusable functions:

- `computeShapiroRangeDelay(...)`
- `computeShapiroTimeDelay(...)`

`fd::perturbations::ThirdBodyGravity` computes:

```text
a_3rd = mu_i * (r_sc_to_i / |r_sc_to_i|^3 - r_central_to_i / |r_central_to_i|^3)
```

`fd::perturbations::SolarRadiationPressure` computes cannonball SRP:

- anti-sunward direction
- configurable `C_R`, area, and mass
- returns `km/s^2`

## Units

Unless otherwise noted:

- distance: km
- time: s
- velocity: km/s
- acceleration: km/s^2
- SRP pressure: N/m^2 internally, converted to km/s^2 at output

## SPICE Notes

- Reusable modules do not call `furnsh_c`.
- Tests and demos load meta-kernels explicitly.
- The canonical repository meta-kernel is lowercase `kernels.tm`; some legacy
  files may still mention `Kernels.tm`.
- SPICE error mode is temporarily set to `RETURN` inside reusable modules that
  query CSPICE, then restored with RAII guards.
- Pure dynamics perturbation geometry uses no aberration correction (`"NONE"`).
- Synthetic observations currently default to light-time correction (`"LT"`) and
  then add Shapiro delay explicitly.
- The Voyager synthetic report is an exception: it uses one CSPICE initial state
  and then RKF45-propagated truth with manual light-time iteration.

## Generated Files

`test_synth_observations` writes:

```text
synthetic_observations_dss43_voyager1.txt
```

This file is useful for inspecting synthetic observation values and for checking
whether light-time, Shapiro, noise, or elevation-mask settings changed output.

`test_voyager_position_od` consumes that report and estimates the Voyager 1
initial state using DP853 propagation with cached SPICE body/station states:
Sun gravity, Jupiter body `599`, Galilean moons `501-504`, Saturn barycenter
`6`, SRP, light-time, Shapiro delay, centered count-time range-rate, and VLBI
differential range. It writes:

```text
tests/voyager_position_estimation_report.txt
tests/voyager_od_postfit_diagnostics_VLBI.csv
tests/voyager_od_trajectory_error_VLBI.csv
tests/voyager_station_observability_windows_VLBI.csv
tests/voyager_wls_iteration_residuals_VLBI.csv
tests/voyager_wls_iteration_summary_VLBI.csv
tests/posterior_covariance.csv
```

The CSV exports are intended for Python plotting. The post-fit diagnostics CSV
contains final residuals and state errors at observation epochs. The trajectory
CSV contains truth-estimate state error at the arc cadence. The observability
CSV contains contiguous station visibility windows from the elevation-masked
observation schedule. The WLS CSVs contain residuals and state-error summaries
for each solver iteration, while `posterior_covariance.csv` stores the final
6-by-6 state covariance. Historical non-VLBI CSVs without the `_VLBI` suffix
may be kept for comparison plots.

Run the companion diagnostics from the repository root after generating the
CSVs:

```sh
python3 tests/plot_wls_iteration_residuals.py
python3 tests/analyze_posterior_covariance.py
python3 tests/plot_vlbi_comparison.py
```

Latest verified `jup310.bsp` run:

```text
Active third bodies: 599 501 502 503 504 6
Prior position error:      70.710678 km
Posterior position error:  5.766723 km
Solver status:             CONVERGED in 5 iterations
Solver weighted RMS:       0.694371
Range RMS:                 150770.580 m -> 9.671 m
Range-rate RMS:            528.089 mm/s -> 0.179 mm/s
VLBI RMS:                  1.598 m -> 1.005 m
```

This is a diagnostic OD test, not a final high-precision solution. The longer
visibility-filtered arc intentionally keeps state truth error as a report
diagnostic rather than a hard test assertion; use residuals and future
conditioning diagnostics before interpreting state accuracy from residual
reduction alone.
