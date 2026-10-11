# Cassini two-week simulated-truth two-way WLS

C++ generates coherent two-way observations from the retained Cassini simulation,
then estimates eight parameters: six initial Cartesian components, an SRP scale
`C_s`, and a thermal-recoil scale `lambda`. Python only reads C++ CSV outputs for
validation, analysis and plotting.

```sh
cmake --build build-clang --target cassini_two_way_wls -j 4
build-clang/cassini_two_way_wls --truth simulation --output results/cassini_two_way_wls
MPLCONFIGDIR=/tmp/deepnav-matplotlib .venv-cassini/bin/python benchmarks/cassini_two_way_wls/analyze_baseline.py
```

Truth comes only from the `multibody_j2_srp_thermal` rows of
`results/cassini_trajectory/trajectory.csv`, interpolated in C++ with cubic Hermite
position/velocity histories. The saved initial state initializes the estimator.
Spacecraft SPK is not used as observation truth or as the error reference.
Planetary, station and attitude kernels remain model inputs. The original retained
simulation is unchanged.

Defaults are the 14-day 9–23 October 2004 arc, 3,220 kg, 12 m², Cr 1.42, and the
retained calibrated XYZ recoil vector. The shared force model includes multibody
tides, direct Saturn J2, uniform rings, attitude HGA SRP and thermal recoil.
Indirect J2 and relativity remain disabled, as in the retained simulation.

The two new coefficients are dimensionless constant scales:

- `a_SRP = C_s * a_SRP_nominal`, including nominal attitude and shadow geometry.
- `a_thermal = lambda * a_thermal_nominal`, preserving recoil direction,
  attitude dependence and radioactive decay.

Both truth coefficients and initial guesses are 1. Gaussian scale priors are
centered at 1 with assumed standard deviations 0.2 and 0.5 respectively. These
are illustrative regularization assumptions, not Cassini calibration errors.
Cartesian initial offsets are [1, −0.8, 0.5] km and [0.01, −0.008, 0.006] m/s;
Cartesian priors are centered on that guess with 2 km and 0.1 m/s uncertainties.
Remaining force coefficients, mass and area are fixed.

The existing C++ `WLS` and `BatchLeastSquaresDriver` solve the nonlinear batch.
RKF45 integrates eight states and the 8×8 variational matrix. The SRP/recoil
sensitivity columns are driven by their nominal accelerations. The complete
radio model is differentiated using affine sensitivity trajectories; accepted
corrections are evaluated with a fresh nonlinear orbit. Large affine coefficient
stencils resolve their weak radio signatures above numerical Doppler precision;
they do not propagate actual orbits with those large coefficients. Covariance
is recomputed without damping. No Cartesian truth samples are fitted by WLS.

Observation assumptions:

- DSS-14/43/63; one station per epoch, choosing the highest minimum elevation
  over both legs and both count endpoints; cutoff 10°.
- 600 s cadence; end-tagged 60 s counted Doppler; reception begins 6000 s after
  the initial spacecraft epoch so retarded spacecraft events stay inside the arc.
- Independent Gaussian noise: 10 m range and 0.01 Hz Doppler; seed 2026.
- Coherent X/X ratio 880/749, 7.175 GHz constant uplink, known 1 µs delay.
- Matching solar Shapiro, 2.3 m zenith troposphere and smooth TDB−TT ground clock.

This is a synthetic visibility schedule, not flight tracking data. The original
force coefficients were previously calibrated against SPK; that prior history
does not make the generated two-way observations independent physical validation.
A good radio fit does not guarantee accurate recovery of all eight parameters.
Posterior uncertainties and correlations show how much the data constrain the
force scales beyond their priors.

Runtime verification checks force derivatives, the augmented STM against a
one-day nonlinear DP853 propagation, radio half-stencil consistency and a final
half-step integration. The half-stencil discrepancy over posterior parameter
uncertainty must stay below 1% of observation noise. Final convergence requires an undamped correction below
0.1 posterior standard deviations. For a trial below that threshold, the line
search permits a 0.1% objective-resolution allowance for counted Doppler;
larger steps must decrease the objective. Python verifies CSV dimensions,
RMS values, covariance positivity, coefficient exports and the saved truth source.

Outputs in `results/cassini_two_way_wls`:

- `observations.csv`: simulation truth, noise and generated two-way data.
- `initial_estimate.csv`: all eight estimates, truth, seed and formal uncertainty.
- `trajectory.csv`: estimated position/velocity, saved truth positions, RTN
  errors and propagated formal position uncertainty.
- `residuals.csv`, `iterations.csv`, `initial_covariance.csv` (8×8), `summary.json`.
- `BASELINE.md`: numerical results, assumptions and validation.
- Overview, residual, convergence, parameter/correlation and 3D PNG/PDF plots;
  standalone `trajectory_3d.html` when Plotly is available.

Short diagnostic:

```sh
build-clang/cassini_two_way_wls --days 0.25 --no-noise --output /tmp/cassini-wls-smoke
```

The short diagnostic checks numerical execution; its weak geometry and finite
prior prevent treating it as an exact parameter-recovery test.
