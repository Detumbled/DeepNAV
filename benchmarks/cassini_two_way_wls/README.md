# Cassini two-week realistic-pass WLS

The retained two-way experiment is `results/cassini_two_way_wls/pass_schedule`,
covering 9–23 October 2004. See its [pass report](../../results/cassini_two_way_wls/pass_schedule/PASS_EXPERIMENT.md),
[paired comparison](../../results/cassini_two_way_wls/pass_schedule/COMPARISON.md)
and [VLBI residual/RMS report](../../results/cassini_two_way_wls/pass_schedule/with_vlbi/VLBI_UNCERTAINTY.md).
Earlier result folders and Monte Carlo experiment outputs were removed.

C++ generates observations, propagates dynamics and performs WLS. Python only
validates, analyzes and plots the C++ exports. Truth is the saved
`multibody_j2_srp_thermal` trajectory in `pass_schedule/truth/trajectory.csv`,
interpolated in C++ with cubic Hermite position/velocity histories. The accompanying
`metadata.json` and `summary.json` preserve its force settings and provenance.
These three files were relocated without changing their contents. Spacecraft SPK
is not the observation truth or error reference; planetary, station and attitude
kernels remain model inputs.

## Model and estimation

The retained model uses multibody gravity, direct Saturn J2, uniform ring gravity,
attitude-dependent HGA SRP and thermal recoil. Mass is 3,220 kg, nominal area
12 m² and Cr 1.42. Indirect J2 and relativistic terms remain disabled.

Eight fitted parameters are the six initial Saturn-relative Cartesian components,
SRP scale `C_s` and thermal scale `lambda`. Both scale truth values and initial
values are 1. Gaussian scale priors have standard deviations 0.2 and 0.5.
Numerical starting offsets are [10, −8, 5] km and [0.10, −0.08, 0.06] m/s.
Prior means retain the original offsets [1, −0.8, 0.5] km and
[0.01, −0.008, 0.006] m/s, with standard deviations 2 km and 0.1 m/s.
Mass, area and other force coefficients remain fixed.

The nonlinear batch integrates the augmented state and variational equations,
repropagates accepted corrections, and recomputes undamped covariance. Convergence
requires a joint undamped correction below 0.1 posterior standard deviations,
with at least two iterations. Both retained fits converged in two iterations.
Runtime checks include force/STM derivatives, radio half-stencil consistency and
half-step integration. No Cartesian truth samples are fitted as observations.

## Tracking assumptions

Fourteen non-overlapping passes use DSS-14/43/63 (5/5/4 passes), targeting up to
six hours each. Every pass midpoint lies in a different elapsed TDB day. Both
light-time legs and both Doppler count endpoints satisfy a 10° elevation mask.

- 4,889 Doppler observations: contiguous 60 s counts at 60 s cadence, 0.01 Hz independent noise.
- 242 range observations: 20-minute cadence, 1 m independent noise plus a shared 10 m Gaussian offset per pass.
- 60 VLBI measurements in the paired fit: jointly visible station pairs, 0.01 m independent path noise (approximately 33 ps).

`PassRangeWhitening` applies the exact covariance block
`R = sigma_ind² I + sigma_shared² 11ᵀ` to range residuals and design rows.
It marginalizes the pass bias without supplying its injected value to WLS or
adding fitted states. Independent VLBI noise leaves the paired radio data identical.
The common-emission VLBI observable is `c × (TDB reception two − reception one)`.
Both links include solar Shapiro and troposphere; station clocks are assumed known.
The two-way model uses an X/X turnaround ratio of 880/749, 7.175 GHz uplink and
known 1 µs transponder delay.

This is a synthetic schedule and noise budget. It does not reconstruct actual
DSN allocations, plasma fluctuations, clock errors or quasar-calibrated delta-DOR.
Force coefficients were previously calibrated against SPK, so the simulation
comparison is not independent physical validation. Formal covariance is conditional
on the assumed observation errors, priors and fixed forces.

| Result | Without VLBI | With 60 VLBI |
|---|---:|---:|
| Position RMS (m) | 1,216.74 | 437.48 |
| Endpoint position error (m) | 1,368.21 | 145.18 |
| Formal position RMS (m) | 1,925.16 | 445.33 |

## Reproduce the retained example

Run from the repository root with the required kernels and saved truth inputs present.
Full estimation commands regenerate the paired fits:

```sh
cmake --build build-clang --target cassini_two_way_wls -j 4
build-clang/cassini_two_way_wls --pass-schedule --range-sigma-m 1 --range-bias-sigma-m 10 --initial-guess-scale 10 --vlbi-count 0 --output results/cassini_two_way_wls/pass_schedule/no_vlbi
build-clang/cassini_two_way_wls --pass-schedule --range-sigma-m 1 --range-bias-sigma-m 10 --initial-guess-scale 10 --vlbi-count 60 --output results/cassini_two_way_wls/pass_schedule/with_vlbi
MPLCONFIGDIR=/tmp/deepnav-matplotlib .venv-cassini/bin/python benchmarks/cassini_two_way_wls/analyze_baseline.py --output results/cassini_two_way_wls/pass_schedule/no_vlbi
MPLCONFIGDIR=/tmp/deepnav-matplotlib .venv-cassini/bin/python benchmarks/cassini_two_way_wls/analyze_baseline.py --output results/cassini_two_way_wls/pass_schedule/with_vlbi
MPLCONFIGDIR=/tmp/deepnav-matplotlib .venv-cassini/bin/python benchmarks/cassini_two_way_wls/plot_vlbi_comparison.py
MPLCONFIGDIR=/tmp/deepnav-matplotlib .venv-cassini/bin/python benchmarks/cassini_two_way_wls/plot_pass_schedule.py
```

Export covariance from the saved VLBI solution and plot residuals and cumulative
position/velocity RMS:

```sh
build-clang/cassini_two_way_wls --export-covariance --pass-schedule --output results/cassini_two_way_wls/pass_schedule/with_vlbi
MPLCONFIGDIR=/tmp/deepnav-matplotlib .venv-cassini/bin/python benchmarks/cassini_two_way_wls/plot_vlbi_uncertainty.py
```

Covariance export repropagates the saved eight-parameter estimate without refitting.
RMS panels have symmetric −3σ to +3σ shading and no separate formal uncertainty
curve. The shading is a covariance reference envelope, not a confidence interval
on the nonnegative RMS. Range, equivalent counted range-rate and VLBI residuals
have no shading.

Each fit retains observation, residual, trajectory, initial estimate/covariance,
iteration and pass CSV exports, `summary.json`, `BASELINE.md`, and PNG/PDF/interactive
3D plots. The VLBI fit additionally retains propagated covariance and RMS analysis.
The shared `truth/` files must be preserved when cleaning generated outputs.

The reusable Monte Carlo core and Cassini adapter remain available; see
[the module documentation](../../docs/monte_carlo.md). No Monte Carlo results are retained.

For historical sampling context, the
[JPL Cassini navigation assessment, §1.3.1](https://descanso.jpl.nasa.gov/DPSummary/DESCANSO17_Cassini_RevA.pdf)
describes Doppler compression and pass-level range weighting. Our cadence and
shared-error covariance are synthetic assumptions.
