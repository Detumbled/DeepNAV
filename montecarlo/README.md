# EKF Monte Carlo

A standalone runner reuses the EKF, RKF45, CSPICE Earth scenario, geometric
radiometric model and DSAC clock simulator. Results go to `Output_montecarlo/`.

```sh
cmake --build build-clang --target ekf_montecarlo test_montecarlo -j4
ctest --test-dir build-clang -R '^test_(ekf|montecarlo(_runner)?)$' --output-on-failure
./build-clang/ekf_montecarlo --runs 100 --seed 2026
./build-clang/ekf_montecarlo --runs 100 --seed 2026 --acceleration-noise-density 0 \
  --output Output_montecarlo/no_process_noise
python3 montecarlo/plot_montecarlo.py --input Output_montecarlo/no_process_noise --no-show
python3 montecarlo/plot_montecarlo.py --baseline Output_montecarlo/no_process_noise --no-show
```

The plotter requires NumPy, Matplotlib and SciPy. For faster ensembles, use a
separate CMake Release build. CSPICE is process-global, so each runner is sequential.
`--output DIR` and `--kernels DIR` override their defaults; the plotter accepts
`--input DIR`. Each study records its actual configuration in `study.json`.

## Initial conditions and pairing

The default initial estimated state is the nominal truth plus independent,
zero-mean Gaussian errors: **100 m per position component, 1 m/s per velocity
component**, 1 microsecond clock bias and `1e-10` fractional frequency. Both
sampled errors and P0 use these sigmas. `--position-sigma-m X` and
`--velocity-sigma-m-s X` change the orbital priors independently. The prior is
currently diagonal; the sampling helper supports full correlations.

The same 8 h orbit, epoch, forces, 60 s measurements and tracking gaps as the
single-run demo are used. Two cases compare paired eight/nine-state filters:

- **matched:** true SRP scale is 1, known to both filters. The ninth state's
  variance is zero; full NEES uses eight stochastic degrees of freedom in both.
- **srp:** each run draws a constant truth scale from `N(1, 0.2^2)`. The eight-state
  filter fixes it at 1; the nine-state filter estimates it, starting at 1 with
  sigma 0.2. Nonpositive physical truth draws fail rather than being discarded.

Initial errors, measurement noise and clock noise use separate reproducible
streams. Each pair receives the same truth and measurements. Contact timing and
station selection are fixed from the nominal trajectory, including elevation
mask and scheduled gaps. CSPICE supplies station states. Forces are Earth point
mass + J2 + solar/lunar differential gravity + eclipsed SRP. Clock parameters
are matched; the clock's existing stochastic simulation and covariance remain.

## Simple orbital process noise

NASA *Navigation Filter Best Practices*, section 3.2.3.1, equations (3.44),
(3.47)-(3.49), describes State Noise Compensation (SNC). We use a single isotropic
white-acceleration diffusion intensity q. Isotropy makes RTN and inertial
intensities equivalent, avoiding a frame rotation in this first implementation.
For a short interval and free-motion STM, its discrete covariance is:

```text
S(dt) = q * [ (dt^3/3) I   (dt^2/2) I ]
            [ (dt^2/2) I       dt I   ]
```

The existing Cartesian propagator evaluates the more general equation (3.44)
by integrating `dS/dt = A*S + S*A^T + B*qI*B^T` alongside the orbit and STM.
It therefore retains gravity effects and position/velocity cross terms rather
than adding an arbitrary diagonal covariance per measurement.

`--acceleration-noise-density X` sets **sqrt(q) in m/s^(3/2)**, not an acceleration
sigma in m/s^2. The default is `1e-6`, hence q = `1e-12 m^2/s^3` =
`1e-18 km^2/s^3`. This is an illustrative tuning choice, not a NASA-prescribed
value. The free-motion per-axis contribution over a two-hour gap is about
0.353 m position sigma and 0.0000849 m/s velocity sigma. Actual orbital dynamics
change this growth. Zero explicitly disables orbital SNC.

Noise is added only to the estimator's orbital covariance. Truth has no added
stochastic acceleration: this is a model-compensation experiment, not a truth
ensemble driven by the same white noise. The SRP scale remains constant with no
parameter process noise. Q does not remove measurement nonlinearities or turn
a constant SRP error into white acceleration; compare consistency and accuracy
against the zero-Q baseline before drawing conclusions about its tuning.

## Outputs and interpretation

Four CSVs (`matched_8.csv`, `matched_9.csv`, `srp_8.csv`, `srp_9.csv`) retain
per-run/per-epoch errors, sigmas, position norm, SRP truth/estimate, observations,
full-state NEES, marginal six-state orbital NEES and pre-update NIS. Errors are
estimate minus truth in km, km/s, clock seconds, fractional frequency and SRP scale.
NEES uses the full covariance with unit-scaled Cholesky solves; no pseudoinverse,
eigenvalue clipping or failed-run removal is used. NIS is NaN during gaps.

The plotter verifies pairing and produces:

- `position_ensemble.png`: median and 5-95% position errors, first-gap zoom,
  and median `3 sqrt(lambda_max(Prr))` reference; this is not a calibrated bound
  on the 3D norm.
- `first_gap_growth.png`: per-run end/start error growth through the first gap.
- `consistency.png`: mean NEES/dof, orbital NEES/6 and NIS/2 with pointwise 95%
  chi-square intervals. Time samples are correlated; intervals are not simultaneous.
- `coverage.png`: marginal +/-3 sigma coverage over runs/components in each block;
  this is not joint state coverage.
- `srp_estimation.png`: SRP error distribution and empirical RMSE versus formal sigma.
- `ensemble.csv`, `summary.json`: aggregate curves and diagnostics.
- `process_noise_comparison.png` and `.json`: optional paired zero-Q/SNC comparison;
  the plotter checks priors, seeds, truth, observations and contact schedule match.

The single-run demo is separate and retains its own prior and zero orbital Q.
Old diagnostic subfolders reflect earlier configurations; read their metadata
before comparing results. Incomplete interrupted experiments are not valid studies.

## Seeded trial results

For 100 paired runs, seed 2026, with 100 m / 1 m/s initial sigmas:

| Case | Final mean NEES, orbital Q=0 | Final mean NEES, SNC sqrt(q)=1e-6 | Expected |
|---|---:|---:|---:|
| Matched, 8/9 states | 8.74 | 2.90 | 8 |
| Uncertain SRP, 8 states | 12504.76 | 10.32 | 8 |
| Uncertain SRP, 9 states | 9.30 | 3.45 | 9 |

The zero-Q matched/nine-state results are within the final pointwise 95%
intervals. This SNC trial is conservative in those cases and is not a calibrated
consistent configuration: it introduces stochastic covariance while orbital
truth is deterministic. For uncertain SRP with eight states it compensates much
of the missing uncertainty, but remains above the 95% upper limit 8.80. Final
position medians change from 5.19 to 6.20 m (matched), 9.35 to 10.25 m (SRP/8),
and 9.13 to 9.46 m (SRP/9). Mean NIS changes from approximately 2.000 to 1.980.
The lower NEES alone therefore does not establish an accuracy improvement.
