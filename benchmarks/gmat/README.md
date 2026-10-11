# Cassini two-way observation benchmark against GMAT R2026a

Run from the DeepNAV directory on this Mac:

```sh
python3 benchmarks/gmat/run_benchmark.py
```

The runner builds `compare_gmat_two_way` in the existing `build-clang`, generates
three **native GMAT `.script` files**, runs the installed `GmatConsole`, and compares
GMAT observations with the actual DeepNAV `TwoWayLink` implementation. No MATLAB
or Python GMAT bindings are needed. GMAT's installed configuration is not edited.
Use `--prepare-only` to generate scripts without executing them. Use `--gmat`
and `--build` to select another installation or configured build directory.

Outputs are under `build-clang/gmat-benchmark/{zero_delay,delay_1us,delay_1ms}/`:

- `cassini_two_way.script`: open and run directly in GMAT's GUI, or run with the console.
- `measurements.gmd`: unmodified native GMAT measurements.
- `states.csv`: GMAT spacecraft and station states at one-second intervals.
- `comparison.csv`: converted GMAT value, DeepNAV prediction, and signed residual.
- `console.log` and `comparison.log`: execution log and maximum residuals.

The parent output directory also contains `delay_comparison.csv`, comparing each
delay case's range increase against zero delay in both tools.

A workspace-local `startup.txt` loads just the relevant plugins and directs logs
into the workspace. The equivalent standalone console command for the zero-delay case is:

```sh
'/Applications/GMAT R2026a/bin/GmatConsole' \
  -s '/Users/adb/Documents/DeepNAV/build-clang/gmat-benchmark/startup.txt' \
  -r '/Users/adb/Documents/DeepNAV/build-clang/gmat-benchmark/zero_delay/cassini_two_way.script'
```

## Scope and assumptions

The link is DSS-14 → Cassini (NAIF -82) → DSS-14. It uses the existing DE442 and
Cassini SPK kernels, with reception tags 2004-10-10 18:13:06.046, 18:23:06.046,
and 18:33:06.046 UTC. Spacecraft propagation reads the SPK; no orbit integration
or force-model comparison is involved. The SPK propagator explicitly starts
`FromSpacecraft`; otherwise GMAT starts at the beginning of the kernel.

The uplink is a constant 7175 MHz, turnaround ratio 880/749, with end-tagged
60-second Doppler counts. The coherent transponder has no independently modelled
spacecraft oscillator. Cases use effective hardware delays of 0, 1 µs, and 1 ms.
Measurement noise, bias, ionosphere, troposphere, relativistic path corrections,
and GMAT's ET-minus-TAI correction are disabled. The DeepNAV comparison uses
geometric light time and coordinate-TDB frequency, without a clock callback.
This does **not** validate Allan noise, frequency ramps, media, Shapiro delay,
or realistic clock corrections.

Station coordinates are approximate DSS-14 Earth-fixed coordinates in km:
(-2353.621251, -4641.341542, 3677.052370). GMAT supplies its own Earth rotation
and station states. Both station and spacecraft trajectories are exported in a
common barycentric MJ2000Eq frame, then cubic-Hermite interpolated in DeepNAV.
This avoids comparing GMAT's station geometry to a separate SPICE station model.
GMAT and DeepNAV still share ephemeris inputs, but **their observation solvers
and Doppler computations are independent**. CSPICE in the comparison executable
only loads the LSK and converts measurement time tags; it does not compute range
or light time.

## Conversions and residuals

GMAT's generic `Range` for the three-node strand is **full round-trip km**.
DeepNAV uses c·RTLT/2, so `gmat_converted = GMAT_Range / 2`. Range residuals in
`comparison.csv` are `(DeepNAV - gmat_converted) * 1000`, in metres.
This deliberately avoids DSN sequential-range units, modulo ambiguities, and
uplink-frequency conversion.

GMAT `DSN_TCP` is negative mean received carrier frequency in Hz. Convert to
DeepNAV's residual using `-GMAT_TCP - (880/749)*7.175e9`. Doppler errors are
`DeepNAV - gmat_converted`, in Hz. GMD tags are TAI modified Julian dates with
JD offset 2430000.0; state tags are TDB modified Julian dates with the same offset.
They must not be interpreted as ordinary MJD (JD offset 2400000.5).

The runner checks execution and requires three measurements of each type.
It deliberately reports residuals **without declaring a numerical validation
pass**. Set scientific acceptance tolerances only after resolving the remaining
Doppler discrepancy and time conventions. Dense exported states also introduce
interpolation and text/epoch precision limits.

For an independent sampling refinement run:

```sh
python3 benchmarks/gmat/run_benchmark.py --state-step 0.5 \
  --output build-clang/gmat-benchmark-refined
```

GMAT's current [Doppler adapter source](https://github.com/nasa/GMAT/blob/main/plugins/EstimationPlugin/src/base/adapter/DopplerAdapter.cpp)
uses a stable path-difference computation, including a Chebyshev-difference option.
Therefore its native TCP need not numerically equal Doppler reconstructed by
naively subtracting two rounded GMD ranges. The exact cause of the observed
Doppler residual has not been established; changing DeepNAV to force agreement
would be premature.

For the generic range convention, see the locally installed GMAT API example
`/Applications/GMAT R2026a/api/Ex_R2020a_RangeMeasurement.py`, which explicitly
compares the two-way `Range` with `2*rNorm`. Measurement settings and TCP conventions
can also be inspected in the installed `docs/help/html/TrackingDataTypes.html`
and `TrackingFileSet.html`.

## Full 24-hour, 260-observation comparison

```sh
python3 benchmarks/gmat/run_benchmark.py --full-day
MPLCONFIGDIR=/tmp/deepnav-matplotlib .venv-cassini/bin/python \
  benchmarks/gmat/plot_full_day.py build-clang/gmat-benchmark-24h
```

`--full-day` reads the saved EKF schedule from
`Output EKF/Cassini two-way/diagnostics.csv`: 145 candidate epochs, 130 tracking
updates, and 15 gaps. It runs a 24-hour GMAT simulation for each station and
retains precisely the original station's range and Doppler at each visible update:
49 DSS-14 pairs, 28 DSS-43 pairs, and 53 DSS-63 pairs. It fails if any requested
observation is missing or duplicated. `--schedule` selects another saved file;
the current full-day mode requires this fixed 2004 arc and its 145/130 schedule.

The full-day case uses a 1 µs transponder delay. DSS-43 and DSS-63 coordinates are
from the DSS table in `Kernels/Cassini/fk/earth_topo_201023.tf`, converted from
metres to km; DSS-14 retains the coordinates used in the short benchmark.
GMAT exports approximately 27.2 hours of states, including the light-time
lookback before the first observation and a margin after the final one.

The original EKF increments its reception epochs in TDB seconds. GMAT's simulator
increments its internal epoch in uniform atomic-time seconds; consequently,
converted TDB tags can differ slightly from that grid. `comparison_260.csv`
records the actual TDB elapsed time and its offset from the original grid.
The runner requires agreement within 1 ms and `summary.json` reports the maximum
actual offset. DeepNAV evaluates each observation at the **actual GMAT reception
tag**, rather than rounding its evaluation epoch to the schedule grid.

Results are stored under `build-clang/gmat-benchmark-24h/`:

- `comparison_260.csv`: all 130 range and 130 Doppler comparisons with station and epoch.
- `summary.json`: counts, mean/RMS/maximum errors, per-station statistics, epoch alignment,
  and a SHA-256 fingerprint of the source schedule.
- `schedule.csv`: the original schedule preserved with the results.
- `DSS-14/`, `DSS-43/`, `DSS-63/`: standalone scripts and unfiltered native outputs.
- `residuals_24h.png` and `.pdf`: residual plots, produced by the plotting command.
- `residuals_24h_by_station.png` and `.pdf`: range and Doppler panels in separate
  station columns, with common axis limits and all 260 selected observations.
- `residuals_24h_DSS-14`, `_DSS-43`, and `_DSS-63` (`.png` and `.pdf`): individual
  station plots, retaining the original tracking gaps and common axis limits.

These are deterministic observations on the **same tracking schedule**, not the
original noisy EKF measurements. The EKF run used solar Shapiro delay, a simple
2.3 m zenith troposphere model, a TT ground-clock callback, and seeded noise;
those terms remain disabled in this benchmark. The comparison tests geometric
observation modelling over the full arc, not corrected measurement or filter
accuracy. To repeat the full-day interpolation refinement, add `--state-step 0.5`
and `--output build-clang/gmat-benchmark-24h-refined`.
