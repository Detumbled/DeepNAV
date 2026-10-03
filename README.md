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
- Iterative light-time and first-order stellar aberration for body targets,
  validated inertial-to-camera attitude rotations, and pinhole camera projection
  with optional OpenCV-compatible radial and tangential distortion.
- Pixel-integrated circular Gaussian image rendering and weighted five-parameter
  centroid fitting with covariance, plus explicit unresolved-body photocenter
  corrections from projected brightness models.
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
- Optical star catalogs, limb/terminator and landmark geometry,
  reflectance models, real-image calibration, optical OD partials/filter
  integration, and visualization work.

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
    CameraAttitude.hpp
    DistortionModel.hpp
    Types.hpp
    core/
      GeometricLineOfSight.hpp
      LightTimeSolver.hpp
      ApparentDirection.hpp
    image/
      CircularGaussian.hpp
      CentroidEstimator.hpp
      PhotocenterCorrection.hpp
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
    CameraAttitude.cpp
    DistortionModel.cpp
    core/
      GeometricLineOfSight.cpp
      LightTimeSolver.cpp
      ApparentDirection.cpp
    image/
      CircularGaussian.cpp
      CentroidEstimator.cpp
      PhotocenterCorrection.cpp
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
plot_centroid_demo.py
  test_ephemeris_interpolator.cpp
  test_geometric_line_of_sight.cpp
  test_optical.cpp
  test_apparent_direction.cpp
  test_centroid_estimator.cpp
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
cmake --build build-clang --target test_apparent_direction -j4
cmake --build build-clang --target test_centroid_estimator -j4
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
ctest --test-dir build-clang -R '^test_apparent_direction$' --output-on-failure
ctest --test-dir build-clang -R '^test_centroid_estimator$' --output-on-failure
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
residual is 5%, configurable by the caller. `fromParameters` and `fromAllanData`
produce custom configurations. `DSAC::shortTermWhiteFmBaseline()` adds the
named DSAC-inspired baseline described below; it is not a fitted hardware model.
Allan data do not determine initial
bias, fractional frequency, their covariance, or the removed drift.

The diffusion convention is
`db = y*dt + sqrt(q_b)*dW_b`, `dy = D*dt + sqrt(q_y)*dW_y`, with independent
Wiener processes. Diffusion intensities are not one-sided PSD coefficients or
per-step standard deviations. The existing exact two-state covariance and
correlated-increment sampler are retained, including singular zero-q_y cases.

`DSAC::shortTermWhiteFmBaseline()` identifies its configuration as
`DSAC_inspired_short_term_white_FM`, with deterministic drift
`D=3e-16/86400 s^-1`, white-FM diffusion `q_b=2.25e-26 s`, and `q_y=0`.
[Burt et al. (2021)](https://doi.org/10.1038/s41586-021-03571-7), abstract,
report the ground-test short-term law
`sigma_A(tau)=1.5e-13/sqrt(tau/(1 s))`. With time in seconds,
`q_b=tau*sigma_A(tau)^2=(1.5e-13)^2*1 s=2.25e-26 s`.
This identifies a white-FM regime rather than inferring it from a single one-day
value. The drift coefficient comes separately from the reported flight drift.
Combining this ground-test noise law and flight drift defines an illustrative
DSAC-inspired preset, not a fit to a complete device operating condition.
The 20-day white-FM continuation is an illustrative extrapolation and does not
reproduce measured long-term behaviour or a stability floor; flicker is omitted.
The old `DSAC::dayMatchedWhiteFmBaseline()` remains available for explicit legacy
comparisons (`q_b=7.776e-25 s`), but is not used by the demo or plotter SETUP.

Initial state and covariance remain separate caller inputs. The comparison
defaults to ideal initial calibration (`b0=y0=0`, `P0=0`) and an ideal ground
clock, with no measurement noise. `localInitial`/`dsacInitial` and
`localP0`/`dsacP0` in the demo can be changed without changing process noise.
At two days this preset predicts approximately 0.062354 ns stochastic bias sigma,
0.056080 m of 3-sigma clock-only range uncertainty, and 0.015541 m of drift-only
range error. These are analytical model quantities, not required recalibration
intervals. Tests verify timestep-independent covariance, sqrt(dt) increments,
noiseless drift evolution, and empirical short-tau ADEV versus the ground-test law.

Published local-oscillator comparators from
[Ely et al. (2025), Table 2](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2025RS008244):

| Reference | Fractional-frequency ADEV at 1 s | ADEV at 1000 s | Fractional-frequency aging/day |
|---|---:|---:|---:|
| Representative OCXO | 5e-13 | 8e-12 | 7e-10 |
| Representative deep-space USO | 5e-13 | 6e-13 | 1e-10 |

Separately, the [AccuBeat manufacturer page](https://www.accubeat.com/uso)
claims ADEV below 5e-13 over averaging times 1–1000 s for its USO and describes
JUICE flight-model delivery. That manufacturer claim is not substituted for
the review table's representative entries.

`LocalOscillator::representativeUsoWhiteFmWithAging()` is the demo preset:
`D=1e-10/86400 s^-1`, `q_b=2.5e-25 s`, `q_y=0`, configuration name
`USO_aging_simplified_white_FM`. It is labelled **representative USO-like preset**:
5e-13 ADEV at 1 s is an illustrative assumption, informed by the review rather
than measured for a particular oscillator. It does not fit a full device spectrum: at 1000 s,
retained drift alone gives about 8.18e-13 ADEV, exceeding the review's 6e-13.
`representativeUsoAgingOnly()` remains available for deterministic comparisons.
For zero initial covariance, DSAC-inspired stochastic range uncertainty is
`sqrt(2.25e-26/2.5e-25)=0.30` times the USO-like uncertainty at every positive
elapsed time. With equal nonzero initial covariance, this ratio applies to the
white-FM process-noise contributions, not necessarily the total uncertainty.
These illustrative presets do not rank real hardware. Aging is
fractional-frequency change per day.

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
creates four separate figures/windows. Range, budget and diagnostics use independent
vertical scales; uncertainty and Allan panels share vertical limits:

- Main figure: signed one-way range error and total pointwise budget
  `c*(abs(mu_bias) + k*sigma_bias_s)` (default k=3).
- Uncertainty figure: stochastic contribution `k*c*sigma_bias_s`, excluding
  the deterministic mean, on a common linear scale.
- Support figure: bias in ns and fractional frequency state y. White FM enters
  bias increments and is excluded from the displayed y state.
- Allan figure: overlapping fractional-frequency Allan deviation, with one
  log-log subplot per clock side by side and common vertical limits for comparison.
  The axes are averaging time tau in seconds and dimensionless sigma_y(tau).
  SETUP supplies explicit q_b and D for dashed theory curves
  `sqrt(q_b/tau + 0.5*(D*tau)^2)`; drift is retained.

The total budget requires explicit deterministic mean histories supplied via
`--mean-csv`, one per input clock, at identical recorded timestamps. Those CSVs
use the same primitive schema and store the expected bias obtained by noiseless
propagation; the noisy realization is never treated as its own mean. Without
mean histories or covariance, the corresponding panel explains the missing
inputs instead of inventing a budget. This pointwise budget is not Allan
deviation or a guarantee for an entire sample path.

`--threshold-m` draws signed thresholds for range and positive thresholds for
budget/uncertainty, annotating the first sampled threshold crossing from all input rows
before display reduction. It reports the preceding sample time with sufficient precision. These are
realization-dependent sampled crossings; an earlier unobserved crossing between samples can be missed.

```sh
python3 -m pip install pandas matplotlib
python3 plot_clocks.py local.csv dsac.csv --mean-csv local_mean.csv dsac_mean.csv --labels Local DSAC --threshold-m 1 --output clocks.png
python3 plot_clocks.py local.csv --time-unit days --threshold-m 1 --output clocks.pdf
python3 plot_clocks.py local.csv --max-points 0 --show --output clocks.svg
```

The plotter uses a noninteractive backend unless `--show` is supplied. PNG,
PDF and SVG are supported. `--output clocks.png` saves the main figure there
and the other figures to `clocks_uncertainty.png` and `clocks_diagnostics.png`.
The additional Allan figure is saved as `clocks_allan.png` when full uniformly
sampled histories are available. The plotter automatically reads sibling
`local_allan.csv`/`dsac_allan.csv` files produced by the demo, or accepts explicit
`--allan-csv` paths in the same clock order. It does not silently use the reduced
display histories for Allan analysis, interpolate nonuniform data, or apply
`--max-points` to the estimator. The estimator uses second differences of bias
and integer averaging factors. Drift is retained and labelled; this is not a
drift-removed noise-only stability curve. Tau is limited to approximately a
tenth of the record duration; long-tau estimates have few independent averages.
Zero ADEV cannot be represented on a logarithmic axis and is not replaced by an
artificial epsilon.

With `--show`, all available windows open together. By default curves display
at most 10000 representative points. `--max-points 0` plots every sample.
Simulation/CSV resolution remains one second and image output remains at 180 DPI.
Display selection preserves endpoints and crossing markers; it can miss narrow peaks.
ADEV and crossing calculations always precede display thinning. Original
CSV data are untouched. Reduced histories (especially an extra final sample)
may be nonuniform and must not automatically be used for Allan estimates or
first-crossing analysis. Keep simulation parameters and seeds in the run's
configuration; filenames do not establish hardware performance.

Build the recording checks with
`cmake --build build-clang --target test_clock_history -j4`. Running
`build-clang/test_clock_history` automatically creates `Output clocks` in the
project root and saves `local.csv`, `dsac.csv`, `local_mean.csv`, and
`dsac_mean.csv`, plus full-step `local_allan.csv` and `dsac_allan.csv`, replacing
previous reference-model outputs. Both realizations cover the same configured interval
and include covariance. The current demo uses 1728000 steps of 1 second (20 days)
and a recording stride of 1: 1728001 samples per history. Edit the clearly labelled
DEMO SETUP constants in `tests/test_clock_history.cpp` to configure duration, step,
threshold, recording stride and seeds; keep the plotter threshold/model parameters in sync.
Display and Allan files contain the same realization; no additional random run
is generated for Allan. Mean files are exported
separately from deterministic propagation. This destination is independent of the working directory;
change `defaultOutputDirectory` in `tests/test_clock_history.cpp` to customize it,
or pass an output directory as a command-line argument. For example:

```sh
./build-clang/test_clock_history
python3 plot_clocks.py
```

With no arguments, the plotter uses its clearly labelled `SETUP` dictionary near
the top of `plot_clocks.py`: CSV/mean/Allan filenames, labels, days on the time
axis, k=3, a 1 m threshold, a 10000-point display limit, explicit theoretical
model parameters, a two-hour USO range inset, image output and window display.
The default output directory is `Output clocks` beside the script, independent
of the current working directory. Edit these settings in one place; set `show`
to False for saving without opening windows. Explicit command-line arguments
remain supported and use the normal CLI defaults instead of SETUP.

CSV validation artifacts remain temporary; only simulation CSVs are retained.
Generated files in `Output clocks` are ignored by Git. Their
parameters and seeds are defined in `tests/test_clock_history.cpp` and printed
at execution, including initial conditions, drift, noise intensities, P0 and
sampling, and the ideal ground-clock assumption. The default comparison uses
the DSAC-inspired short-term white-FM baseline and representative USO-like preset;
the other numerical tests retain synthetic validation coefficients. Plot labels
distinguish these approximations. Mathematical Allan checks do not establish
complete OCXO/DSAC hardware performance or ESA suitability.
Python tests
also exercise C++ exports when that executable is available.

The result is a simplified clock holdover comparison. Figure footnotes state ideal
initial bias/frequency calibration (`b0=y0=0`, `P0=0`), ideal ground reference,
constant uncompensated aging, simplified white FM, and omitted flicker/other
long-term noise. The illustrative 1 m clock-only allocation is not a universal
navigation requirement. Clocks are never reset at crossings; these do not impose
ground-contact intervals or represent first-passage probabilities.

| Days | DSAC deterministic range [m] | DSAC stochastic 3 sigma [m] | Total budget [m] |
|---|---:|---:|---:|
| 2 | 0.01554 | 0.05608 | 0.07162 |
| 10 | 0.38853 | 0.12540 | 0.51393 |
| 20 | 1.55412 | 0.17734 | 1.73146 |

For ideal initial conditions, the analytical DSAC budget reaches 1 m at
14.77008 days; its stochastic-only envelope stays below 1 m over the 20-day run. The plot labels
analytical budget crossings separately from sampled realization/envelope crossings.
At 20 days USO deterministic range is about 518041 m and stochastic 3 sigma
about 0.591 m. Its deterministic-only 1 m crossing is approximately 2400.83 s;
the actual noisy sampled crossing can differ. Extending the run does not validate
hardware over 20 days; plotted Allan tau is at most two days for this run.
Slide caption:

> Illustrative one-way ranging clock budget: USO-like and DSAC-inspired models
> with white frequency noise, initial calibration uncertainty, and deterministic
> aging. Flicker noise is omitted; long-term stochastic behaviour is extrapolated.

The default demo has zero initial calibration covariance; nonzero initial
uncertainty is propagated separately when supplied. This qualifier is shown
in figure footnotes so the caption does not imply a nonzero calibration error.

Allan analysis follows the bias-second-difference definition in
[NIST SP 1065](https://nvlpubs.nist.gov/nistpubs/Legacy/SP/nistspecialpublication1065.pdf).

### Optical Navigation

The optical-navigation layer predicts body-target pixel coordinates following
Owen, *Spacecraft Optical Navigation*, Sections 4.2 and 4.5-4.8. Geometry and
attitude calculations depend only on Eigen; ephemerides enter through
`fd::dynamics::StateProvider` and may be analytic or supplied by the caller.

The complete body-target prediction path is:

```text
Geometric barycentric states -> light-time iteration -> stellar aberration
                            -> inertial-to-camera rotation -> distorted pixel
```

Use `fd::dynamics::CartesianState`: position is in kilometres and velocity in
kilometres per second. The reception epoch is the exposure midpoint, in TDB
seconds past J2000. Observer and target states must be geometric states in the
same barycentric inertial frame (for example J2000/ICRF with SSB origin).
Observer velocity must be barycentric, not relative to the target. The raw
observer state has no frame metadata, so the caller is responsible for matching
it to the target provider. Do not supply states already corrected for light time
or aberration. The provider must cover the retarded emission epochs.

`fd::opnav::LightTimeSolver` holds the observer at reception and iterates the
emission epoch of the target:

```text
tau = |S(t - tau) - R(t)| / c
T   = S(t - tau) - R(t)
```

Its defaults are a 1e-9 s convergence tolerance and eight iterations. The result
contains reception/emission epochs, target state at emission, the retarded LOS
vector, iteration count, and the absolute fixed-point residual in seconds.
`lightTimeSeconds` is the tau used to evaluate the returned emission state;
`|T|/c` agrees within the reported residual. Accuracy is also limited by the
provider and floating-point epoch resolution. Invalid epochs/options/states,
zero range, and failure to converge throw exceptions. This is a geometric
light-time calculation and does not include gravitational delay.

`fd::opnav::core::computeApparentDirection(...)` wraps this solver and adds the
paper's first-order stellar aberration (Eq. 4.3):

```text
A = T + tau * observer.velocityKmPerSec
```

It returns the light-time solution, apparent vector in kilometres, and its unit
direction. The vector magnitude is not an optical range measurement; camera
projection depends on its direction. This follows the paper's Newtonian
approximation for ordinary spacecraft velocities, without relativistic
aberration or gravitational light bending.

`fd::opnav::CameraAttitude` accepts a proper 3x3 rotation matrix `C` with the
explicit convention `directionCamera = C * directionInertial`. Its rows are
camera axes expressed in inertial coordinates. Camera +Z is the viewing axis.
A camera-to-inertial matrix must be transposed before construction. The matrix
is checked once for finite entries, orthogonality, and determinant +1 (absolute
tolerance 1e-12), then reused for every target in that exposure. Attitude
retrieval, estimation, interpolation, and spacecraft-to-camera mounting
calibration remain caller responsibilities; compose those rotations into `C`.

`fd::opnav::CameraModel::project(directionInertial, attitude)` rotates and uses
the existing camera-frame projection. The original `project(directionCamera)`
overload remains available. Both apparent vectors and unit directions can be
projected. The camera uses gnomonic projection, identity or OpenCV-compatible
rational radial/tangential distortion, and a 2x2 focal-plane-to-pixel matrix
plus principal point. It rejects non-finite directions and camera-frame Z <= 0.
Projection does not clip to image bounds. The book's optional detector `xy`
terms and additional misalignment models are not implemented.

```cpp
#include "dynamics/LinearStateProvider.hpp"
#include "opnav/core/ApparentDirection.hpp"
#include "opnav/CameraModel.hpp"

using namespace fd::opnav;
const TdbEpoch exposure(0.0);
const fd::dynamics::CartesianState observer{{0, 0, 0}, {0, 30, 0}};
const fd::dynamics::LinearStateProvider target(
    exposure, {{0, 0, 299792.458}, {0, 0, 0}}, "J2000", "SSB");
const CameraAttitude attitude(Eigen::Matrix3d::Identity());
const CameraModel camera(CameraIntrinsics{
    50.0, {1000.0, 500.0}, 20.0 * Eigen::Matrix2d::Identity()});
const auto apparent = core::computeApparentDirection(exposure, observer, target);
const auto pixel = camera.project(apparent.vectorKm, attitude);
// About (1000, 500.100069) pixels: aberration displaces the image toward +Y.
```

`core::computeGeometricLineOfSight(...)` remains a simultaneous-epoch helper
returning vector, range, and unit direction. It applies no light time or
aberration. Both geometric and camera projection tests remain in place.
`test_apparent_direction` links Eigen only and checks analytic stationary and
moving-target light time, convergence failure, invalid inputs, aberration sign,
rotation direction and inverse, projection scale invariance, and the complete
state-to-pixel path. Fixed-size vectors/matrices avoid per-target heap storage;
no numerical differentiation or extra ephemeris passes are needed.

This completes the basic body-center prediction path. The image estimator below
provides local centroid measurements and covariance. Star catalog propagation,
limb/terminator geometry, body-fixed landmarks, general reflectance models,
optical OD derivatives, and filter integration remain future work. Existing distortion support is a selectable practical model rather
than an implementation of every camera calibration model in the book.

### Image Centerfinding: Sections 5.2 and 5.4

`opnav/image` implements a local image estimator for one isolated star or
unresolved body. Inputs are calibrated DN grids with detector bias, flat-field,
and readout-smear corrections already applied, plus known independent pixel
variances in DN squared. Calibration and variance estimation are upstream
responsibilities. This is a per-image fit, not a temporal navigation filter.

`CircularGaussian` uses the standard deviation convention

```text
I(s,l) = h * exp(-((s-sc)^2 + (l-lc)^2)/(2*sigma^2))
DN(pixel) = b + integral_over_pixel(I)
parameters = [sc, lc, h, sigma, b]
```

Here `h` is continuous peak intensity, `sigma` is in pixels (FWHM is about
2.355*sigma), and `b` is constant background DN per pixel. Total source signal
is `2*pi*h*sigma^2`. This explicit convention differs from the book's printed
Gaussian width/amplitude normalization; the pixel-integrated fitting approach
is the same. Pixel centers have integer sample/line coordinates, boundaries
are +/-0.5, and `Image(line, sample)` stores the grid in row-major order.
An optional origin preserves full-image coordinates when fitting a cutout.

`renderCircularGaussian` integrates each pixel with error functions, rather than
sampling the intensity at its center. `simulateCircularGaussian` adds Poisson
shot noise and Gaussian read noise with a seed. Gain is electrons per DN, read
noise is in DN, and the returned expected variance is
`meanDN/gain + readNoiseDN^2`. These simulated variances use the injected model;
a real image needs its own calibrated noise estimate. Seeds reproduce draws
within the same C++ standard-library implementation.

`fitCircularGaussian` jointly estimates center, height, width, and background
using weighted nonlinear least squares with analytic pixel derivatives and
Levenberg-Marquardt damping. Logarithmic height/width keep them positive;
scaled 5x5 normal equations limit numerical conditioning problems. One-dimensional
pixel integrals are cached per axis and their storage is reused each iteration.
Without an initial model, border-median background and positive signal moments
provide a starting point. Supply a local window containing one source; this is
not full-image source catalog detection or deblending.

A nonzero `PixelMask` entry excludes a bad or saturated pixel. Unmasked DN and
variances must be finite, variances strictly positive, and more than five pixels
must remain in a window of at least 3x3. The fitted center must stay inside the
window; trial widths are limited to 0.05 pixels through twice its largest
dimension. Windows should include source wings and background. The default
iteration limit is 80, step tolerance is 1e-7, and minimum fitted height SNR is 5.

On success, `GaussianFitResult` contains the model, full five-parameter
covariance, and a `PixelMeasurement` with sample/line coordinates and the full
2x2 covariance in pixels squared. This is the inverse local WLS information,
including nuisance-parameter coupling. Supplied pixel variances are treated as
known, so covariance is not multiplied by the reduced chi-squared. The result
also reports chi-squared, degrees of freedom, and iterations. Inspect reduced
chi-squared for model mismatch; fit convergence alone does not establish image
quality. PSF errors, illumination errors, correlated detector noise, and attitude
uncertainty are not included in this covariance.

`FitStatus` distinguishes `Converged`, `NoSignal`, `NonConverged`, and `Singular`.
Only successful, sufficiently significant, observable fits return a measurement
and covariance. Invalid grids/options throw exceptions. NoSignal may be reported
for a converged fit below the configured height SNR. Strongly undersampled
sources can remain unobservable: a literal one-pixel impulse contains no
recoverable subpixel information.

Section 5.4 uses the same fit to find an unresolved body's **photocenter**.
`brightnessPhotocenterOffset` implements the discrete first moments of
Eqs. 5.16-5.17 using a caller-supplied nonnegative, background-free projected
body-brightness grid and its known geometric center. Coordinates and spacing
must use the same pixel units as the measurement. Its offset is light center
minus geometric center. This model grid is intrinsic body brightness, not a
noisy observed image; its shape/illumination comes from the caller. Numerical
first moments avoid reliance on the printed spherical shortcut in Eq. 5.18.
General reflectance laws and a physical moon renderer are not yet implemented.

`correctPhotocenter` subtracts a supplied `PhotocenterOffset` and returns a
separate geometric-center measurement. It preserves the raw photocenter and
adds offset-model covariance under an explicit independence assumption. The
brightness-grid helper treats its model as exact and defaults that covariance
to zero; supply nonzero uncertainty when the offset is uncertain. Apply the
correction once, and compare the corrected observation with the predicted body
center. Alternatively retain the raw photocenter and offset the prediction,
but do not do both. A geometric center represents center of mass only when the
body model justifies that assumption; a Gaussian fit cannot establish it.

```cpp
#include "opnav/image/PhotocenterCorrection.hpp"
using namespace fd::opnav::image;

const CircularGaussian source{{5.25, 5.4}, 2000.0, 0.7, 20.0};
const auto synthetic = simulateCircularGaussian({11, 11}, source, {2.0, 3.0}, 42);
const auto fit = fitCircularGaussian(synthetic.dn, synthetic.varianceDn2);
if (fit.measurement) {
    const PhotocenterOffset offset{Eigen::Vector2d{0.25, 0.0},
                                  0.0001 * Eigen::Matrix2d::Identity()};
    const auto geometric = correctPhotocenter(*fit.measurement, offset);
    // Attach exposure midpoint, target/camera IDs, attitude, and quality metadata
    // before passing geometric.center and geometric.covariance to an OD filter.
}
```

The estimator does not attach target identity or time itself; it supplies the
local measurement component. OD measurement construction, attitude uncertainty,
optical Jacobians, and sequential-filter ingestion remain future work.

`test_centroid_estimator` links Eigen only. It checks pixel integration against
independent numerical quadrature, analytic derivatives against finite differences,
flux conservation, five-parameter noiseless recovery over subpixel phases and
PSF widths, masked pixels, and failure cases. A 200-image seeded Monte Carlo
compares measured bias/scatter with reported covariance and detector-noise
residuals. An end-to-end test uses the apparent-direction and camera models to
place a synthetic unresolved moon, then independently fits its noisy image.
The synthetic moon is a symmetric point source with zero physical photocenter
offset; a separate asymmetric brightness-model test checks correction sign and
propagation of offset uncertainty. These tests validate the assumed Gaussian
and noise models, not general real-camera or moon performance.

To export and inspect that demonstration:

```sh
./build-clang/test_centroid_estimator /tmp/centroid_demo.csv
python3 plot_centroid_demo.py /tmp/centroid_demo.csv /tmp/centroid_demo.png
```

The figure shows pixel intensities, injected/fitted centers, and a joint 95%
local Gaussian uncertainty ellipse. Normal CTest runs generate no output files.

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
