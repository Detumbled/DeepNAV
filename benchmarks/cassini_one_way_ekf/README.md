# Cassini one-way onboard EKF

This benchmark generates Earth-to-spacecraft pseudorange and 60-second counted
Doppler in C++, and runs the existing forward EKF against the retained two-week
Cassini simulated trajectory (9–23 October 2004). Python validates exported CSV
files and plots results. It does not generate observations or estimate states.

## Reused modules

- `CassiniForceModel`, `CassiniAttitude`, `Sat427Gravity`: unchanged retained multibody, direct J2, uniform ring, attitude HGA SRP and thermal-recoil model, with mass 3,220 kg, area 12 m², Cr 1.42 and fixed nominal force scales. Indirect J2 and relativity remain disabled.
- `SavedSimulation`: the existing WLS truth reader, extracted into a shared header; only `multibody_j2_srp_thermal` rows are used. No spacecraft SPK state is used for estimation or the accuracy reference.
- `OneWayLink`, `OneWayMeasurement`, `EndpointClock`, library light-time solver, solar Shapiro and troposphere corrections.
- `EKF`, `CartesianPropagator`, RKF45 orbit/STM integration and Joseph covariance updates.
- `DSAC`, `ClockTruthSimulator`, `ClockModel`: existing clock noise, known deterministic aging, clock STM and exact discrete process covariance.
- The retained realistic pass CSV supplies target windows and station assignments. New uplink visibility checks determine the actual observation counts.

`ReceptionTrajectory` is a short local orbit/STM adapter for the receiver during
an end-tagged count. It integrates backward and forward around the current EKF
epoch. Its affine trajectory factory avoids repeated nonlinear integrations for
measurement finite differences. Large symmetric stencils resolve GHz count
arithmetic; half-sized stencils are checked every 500 updates. A unit test checks
the local trajectory and STM against independent nonlinear propagation.

## State and assumptions

The eight states are `[rx, ry, rz, vx, vy, vz, b, y]`, in km, km/s, seconds and
fractional frequency. Clock bias is clock minus the nominal reference. SRP and
thermal scale factors are fixed at 1. There is no orbital process noise or VLBI.

Two cases are run:

- `ideal`: exact clock, zero clock bias/frequency covariance and zero clock process noise. Orbit states are still estimated by the eight-state EKF.
- `dsac`: estimated bias/frequency, using `DSAC::shortTermWhiteFmBaseline()`. Truth starts at 2 µs bias and 2×10⁻¹² fractional frequency; the estimate starts at zero with 10 µs and 10⁻¹⁰ prior standard deviations. Independent clock noise is generated at 60 s intervals, including tracking gaps.

The DSAC-inspired white-FM model has bias diffusion 2.25×10⁻²⁶ s, no random-walk
frequency diffusion, and known frequency aging 3×10⁻¹⁶/day. White-FM bias increments
are linearly interpolated within a 60 s interval; deterministic frequency aging
is evaluated exactly. A fourteen-day white-FM continuation is illustrative, not
a flight-calibrated long-term oscillator model. Bias diffusion also perturbs the
counted observable; this first filter uses the local deterministic clock model
within a count. At the default carrier/noise, the white-FM count contribution is
approximately 0.00014 Hz, compared with 0.01 Hz receiver noise. The small induced
state/measurement noise correlation is not represented in this baseline.

Both clocks use the existing smooth TDB−TT nominal offset, with spacecraft clock
errors defined relative to this TT-equivalent reference. This is a controlled
clock/navigation baseline; complete spacecraft proper-time frequency transfer
and stochastic ground clock errors are not modeled. Background ephemerides,
attitude, ground time references and uplink carrier schedule are assumed available
to the onboard estimator.

The orbit estimate starts at [1, −0.8, 0.5] km and
[0.01, −0.008, 0.006] m/s offsets, with the same 2 km and 0.1 m/s prior standard
deviations used in the two-way WLS. Unlike a batch optimizer starting point,
an EKF initial estimate is its prior mean. `--initial-guess-scale` changes that
mean for an explicit stress test; it does not change the prior standard deviations.

The filter receives the raw spacecraft clock tag. It infers a coordinate epoch
from its current clock estimate, propagates to that epoch and evaluates the
one-way measurement at the raw tag. True coordinate times and true clock states
are used only by the simulator and diagnostics, never to pre-correct filter tags.
No WLS estimate initializes the EKF. Updates are not iterated and NIS gating is
disabled so numerical/model problems remain visible.

## Schedule and measurement errors

The original fourteen realistic pass windows target DSS-14/43/63. Each candidate
uplink is checked at both count endpoints with the station at its actual retarded
transmission epoch and a 10° elevation mask. Some target-window observations may
be rejected; reports use actual counts. Doppler is sampled every 60 s with 60 s
receiver-clock counts; range is activated every 20 minutes relative to each pass
start. The X-band carrier is 7.175 GHz, with no coherent two-way turnaround ratio.

One-way range noise is independent 10 m Gaussian noise, and Doppler noise is
independent 0.01 Hz Gaussian noise. These are declared synthetic one-way assumptions,
not a Cassini instrument calibration. The previous two-way pass biases are not
copied, and this experiment does not reproduce their correlated covariance.
Named independent RNG seeds keep instrument noise identical across clock cases;
clock noise uses a separate stream. Range and Doppler covariance is diagonal.

## Build and run

From the repository root, with the existing kernels, retained truth and pass CSV present:

```sh
cmake --build build-clang --target cassini_one_way_ekf -j 4
build-clang/cassini_one_way_ekf --output results/cassini_one_way_ekf
MPLCONFIGDIR=/tmp/deepnav-matplotlib .venv-cassini/bin/python benchmarks/cassini_one_way_ekf/analyze.py
```

Options: `--days` (whole days, 1–14), `--clock both|ideal|dsac`, `--max-step-s`
(up to 60), `--range-sigma-m`, `--doppler-sigma-hz`, `--initial-guess-scale`,
`--clock-frequency-sigma` (initial fractional-frequency prior standard deviation,
default 10⁻¹⁰), `--seed`, `--no-noise`, and `--output`. `--no-noise` disables instrument noise;
the DSAC case still has its stochastic clock. Short checks should use `/tmp`:

```sh
build-clang/cassini_one_way_ekf --days 1 --no-noise --output /tmp/cassini-one-way-check
ctest --test-dir build-clang -R '^(test_reception_trajectory|test_one_way_link|test_ekf|test_clocks|test_cassini_one_way_ekf_smoke)$' --output-on-failure
```

## Outputs and interpretation

Each case exports `observations.csv`, `trajectory.csv`, `innovations.csv`,
`covariance.csv` (full 8×8 covariance), `initial_prior.csv`, `clock_truth.csv`, and
`summary.json`. Analysis exports source/input/output hashes, validation statistics,
RMS analysis, PNG/PDF plots, a comparison CSV and the shared `REPORT.md`.

The one-day 60 s versus 30 s step check changed the solutions by at most
16.6 m (ideal) and 40.1 m (DSAC-inspired), below 0.32% of formal position
uncertainty. This is adequate for this initial kilometre-scale benchmark, but
metre-level conclusions require tighter numerical validation. See
`results/cassini_one_way_ekf/numerical_checks.json`.

Orbit diagnostics use the nominal 600 s grid, including gaps and the initial
prior, plus the last 60 s epoch; the final point is one minute before the saved
truth boundary to leave room for clock-dependent epoch assignment. Error is
compared with simulation truth at the filter's coordinate epoch. RMS panels have
symmetric −3σ to +3σ cumulative covariance reference shading without a separate
formal uncertainty curve. This envelope is not a confidence interval on RMS.
Orbit-error and pre-update innovation panels are unshaded. Clock panels have
component ±3σ bands. Equivalent range-rate innovations use `−c/f` times frequency
innovations; the CSV preserves the original Hz observable.

NIS uses one observable at Doppler-only epochs and two at range/Doppler epochs.
NEES uses six active states for the ideal-clock case and eight for the stochastic
case. These are single-run consistency diagnostics, not independent ensemble
samples. An EKF has only past data; full-arc WLS has future observations and
different observables, so their RMS values are context rather than a controlled
one-way/two-way performance comparison.

Sources: [DSAC mission results, NASA](https://ntrs.nasa.gov/citations/20230006982),
[Burt et al. (2021), clock stability measurements](https://doi.org/10.1038/s41586-021-03571-7).
The existing DSAC module documents which measurements inform its illustrative
short-term white-FM configuration.

## Clock-frequency prior sensitivity

The completed [frequency-prior sweep](../../results/cassini_one_way_ekf/frequency_prior_sweep/REPORT.md)
compares initial frequency standard deviations 10⁻¹⁰, 10⁻¹¹ and 10⁻¹². Only this
prior uncertainty changes; its mean remains zero, while truth starts at 2×10⁻¹².
Clock-bias uncertainty remains 10 µs. The original 10⁻¹⁰ run is copied into the
sweep with its provenance preserved; the original baseline stays intact.

For this one realization, position RMS is 13.560, 8.025 and 8.274 km respectively.
Endpoint errors are 17.534, 5.326 and 3.379 km. The tightest prior reduces the
along-track/frequency correlation to −0.307, but leaves a frequency error of
−2.02×10⁻¹², approximately two posterior standard deviations. A smaller
covariance does not establish accurate calibration. No new operational prior
is selected from this synthetic sensitivity experiment; the default is unchanged.

Reproduce the two tighter-prior C++ fits and analyze the comparison:

```sh
build-clang/cassini_one_way_ekf --clock dsac --clock-frequency-sigma 1e-11 --output results/cassini_one_way_ekf/frequency_prior_sweep/1e-11
build-clang/cassini_one_way_ekf --clock dsac --clock-frequency-sigma 1e-12 --output results/cassini_one_way_ekf/frequency_prior_sweep/1e-12
MPLCONFIGDIR=/tmp/deepnav-matplotlib .venv-cassini/bin/python benchmarks/cassini_one_way_ekf/analyze_frequency_sweep.py
```

For a fresh sweep without the copied baseline, generate the 10⁻¹⁰ case with the
same C++ command and `--clock-frequency-sigma 1e-10 --output
results/cassini_one_way_ekf/frequency_prior_sweep/1e-10`. Individual fit validation
and plots use `analyze.py --output <sweep-case-directory>`. The sweep analyzer
verifies exact observation and clock-truth identity, unchanged remaining priors
and settings, numerical derivative checks, and positive covariance. It exports
errors, cumulative RMS with symmetric ±3σ reference shading, clock errors,
comparison CSV and the shared report.
