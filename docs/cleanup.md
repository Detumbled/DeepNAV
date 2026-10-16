# Repository cleanup

## Source layout and build dependencies

Removed the vendored graphics libraries and the Try3D viewer. The numerical
build now depends on Eigen and CSPICE without a graphics stack. Console demos
live in `demos/`, Python plotters in `scripts/`, and the two former root-level
tests in `tests/`. Plotter defaults, imports and CMake source paths follow the
new locations. The root README is a short build/layout guide.

`mission_analysis/` provides include/source/test folders for future module work;
it does not add an empty build target or implement mission-analysis algorithms.

The entries below describe earlier result cleanups, not this source-layout change.

Validation: the full build succeeds, all 13 Python clock-plotter tests pass, and
34 of 35 CTest cases pass. The remaining `test_rkf45` accuracy failure also
reproduces when independently compiling the pre-cleanup HEAD test and integrator:
position closure error is about 0.01177 km over its 60-period propagation. The
relocated test and integrator mathematics were not changed during this cleanup.

Removed prior generated results from Output EKF, Output clocks, Output optical,
Output_montecarlo, tests, montecarlo/data, build benchmark result folders, and
the GMAT result reports. Removed the Cassini one-way
C++ demo, its plotting script, study runner, documentation and CMake targets.
Reusable one-way radio models and tests remain; the batch helper now lives in
include/filters/OneWayBatchBaseline.hpp. New trajectory results go under results/.

A recoverable copy of removed files is temporarily stored outside the repository at:

/var/folders/ny/9gsmq_q12v31spm3nfzgjzz80000gn/T/deepnav-previous-results-oeec0joo

Retained only the nominal SAT427/PDS-3220-kg Cassini trajectory. Removed the
alternate-mass/recoil outputs and study driver from the repository; numerical
validation of the retained physical model remains. Recoverable archive:

/tmp/deepnav-sensitivity-archive-8wnawsoh

Replaced the superseded single-fit two-way WLS output with the paired
`results/cassini_two_way_wls/no_vlbi` and `with_vlbi` fits plus comparison figures.
The original simulated truth in `results/cassini_trajectory` remains unchanged.
The old single-fit outputs are recoverable at:

/tmp/deepnav-pre-vlbi-baseline-onx98r6m

## Retained realistic-pass results

Only `results/cassini_two_way_wls/pass_schedule` is retained, including both paired
fits and their existing plots. The three unchanged simulation truth files were
relocated into its `truth/` subfolder; consumers and input provenance paths now
point there. Earlier trajectory, continuous-tracking, initial-guess and Monte Carlo
result folders were removed. Original numerical fit results and run source hashes
were preserved; no WLS refit was performed during this cleanup.
