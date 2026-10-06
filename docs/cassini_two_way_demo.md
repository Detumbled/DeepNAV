# Cassini coherent two-way console demo

Build and run from any working directory:

```sh
cmake --build build-clang --target cassini_two_way_demo test_two_way_link -j 4
./build-clang/cassini_two_way_demo
ctest --test-dir build-clang -R 'test_two_way_link|test_cassini_two_way_benchmark' --output-on-failure
```

The demo loads absolute paths from `DEEPNAV_SOURCE_DIR/Kernels`, including the
local Cassini September–October 2004 reconstructed SPK. It loads the Cassini SPK
after DE442 so the reconstruction's associated planetary ephemerides take
precedence. It needs the local Cassini and DSN kernel set described in
[`Kernels/Cassini/README.md`](../Kernels/Cassini/README.md); it downloads nothing.
It unloads its own kernels and does not clear an existing CSPICE pool. CSPICE
uses process-global state; do not run its kernel operations concurrently.

The default end reception is `2004-10-10T18:13:06.046` UTC, with three counts
spaced by 600 TDB seconds. Each count spans 60 TDB seconds. The station is fixed
within a count; automatic selection chooses the highest minimum elevation
among DSS-14/43/63 at uplink/downlink and start/end epochs. Counts below the
10-degree mask are skipped; an entirely invisible run returns an error.

## Event chain and observables

States are geometric J2000 vectors relative to the solar system barycenter,
queried with `NONE`. An Earth-centered origin moving between emission and
reception is unsuitable for subtracting these states.

Given Earth reception `t3`, solve backwards:

1. Cassini downlink transmission `t2out = t3 - tauDown`.
2. Cassini uplink reception `t2in = t2out - transponderDelay`.
3. Earth uplink transmission `t1 = t2in - tauUp`.

Both light times are solved independently. The downlink receiver is the station
at `t3`; the uplink receiver is Cassini at `t2in`. The station is evaluated at
`t1` on uplink, including Earth's motion and station rotation during the roughly
2.5-hour round trip. The uplink is not inferred by doubling the downlink.

The printed range is the coordinate-time convention `c * RTLT / 2`, with
`RTLT = tauUp + transponderDelay + tauDown`. It includes the configured delays;
it is not a calibrated DSN range in range units or either instantaneous distance.

The ideal locked transponder uses the Cassini X-up/X-down carrier ratio
`k = 880/749`. See [JPL, Station Delay Calibration for Ranging Measurements,
2009](https://tda.jpl.nasa.gov/progress_report/42-177/177C.pdf).
The ratio rescales carrier phase, not range or travel time. The default
7.175 GHz uplink is an illustrative X-band frequency, not a recovered Cassini
frequency record.

Doppler is a finite phase count, end tagged at `t3`:

```text
receivedCycles = k * integral[fUp(u) du] over the corresponding Earth TX interval
meanReceivedHz = receivedCycles / Earth RX clock interval
residualHz = meanReceivedHz - constantReceiverReferenceHz
```

Positive residual means received frequency above the declared receiver reference.
This sign is explicitly a demo convention, not an ODF data convention. The
transmitter ramp is linear and phase continuous; its reference epoch and slope
are measured on the selected ground clock. The receiver reference stays constant
at `k * fUp(referenceEpoch)`, including when the transmitter ramp is enabled.
The default station clock follows TT via a `SpiceTdbClock` instance, which loads
its constants once after kernel loading and uses the LSK's smooth periodic
TDB-TT model. The executable checks that offset against CSPICE `unitim` to within
absolute-epoch rounding. The count length supplied by `--count-s` is in TDB;
the received carrier is divided by the corresponding TT interval.

The implementation differences short round-trip delays to obtain the transmit
interval. Subtracting two absolute emission epochs would amplify epoch rounding
into artificial Hz-level errors at X-band. It also computes the residual before
adding the GHz receiver reference. No onboard free-running oscillator is injected
into this ideal coherently locked link.

## Two light-time backends and corrections

The benchmark compares:

- CSPICE `spkezr(..., "CN", ...)` for each reception leg.
- The existing `opnav::LightTimeSolver`, via `libraryReceptionSolver`.

Both use the same geometric SPK states and the same transponder delay. Neither
benchmark includes propagation corrections. `CN+S` is unnecessary for radio
travel times. CSPICE's `CN` does not include gravitational delay:
[NAIF SPKEZR documentation](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/spkezr_c.html).

A third solution iterates the library light time with additive corrections:

- First-order point-mass solar Shapiro, gamma=1, using solar GM from the PCK
  and retarded endpoint states. This is not a superior-conjunction-grade model.
- Non-dispersive tropospheric path `zenithPath / sin(elevation)`, default 2.3 m.
  This is an illustrative dry-path mapping, not measured atmospheric calibration.

The transponder delay defaults to **1 microsecond**, explicitly illustrative,
constant and non-dispersive in TDB seconds. It is not a Cassini flight-calibrated
hardware delay. Set it to zero for an ideal instantaneous transponder.

Propagation delay callbacks receive both endpoint providers and solved endpoint
states. They can compose or replace these models. Group and carrier phase delay
must be supplied separately for dispersive media; the same non-dispersive
configuration is suitable for the present demo. Different transponder group and
phase delays can likewise use separate `TwoWayLink` instances.

## Reusable modules

| Module | Responsibility |
|---|---|
| `dynamics/SpiceStateProvider` | Geometric J2000/SSB SPK states; caller owns the kernel pool |
| `radiometric/TwoWayLink` | Event chain, transponder ratio/delay, linear ramp, finite phase count, ground clock callback |
| `radiometric/PropagationCorrections` | Point-mass Shapiro and zenith-path mapping |
| `radiometric/SpiceLightTime` | CSPICE CN backend, DSN topocentric elevation, smooth TDB-TT clock offset |
| `cassini_two_way_demo.cpp` | Local kernel loading, CLI, station selection and console diagnostics |

`deepnav_radiometric` links only Eigen and the existing numerical light-time
solver; it has no CSPICE or graphics dependency. `deepnav_radiometric_spice`
adds only the CSPICE adapter, without linking the application's graphics stack.
Mission-specific ratio and frequency are supplied explicitly by the caller;
there are no Cassini radio defaults in the reusable interfaces. Neither target
depends on the Voyager demo, which is intended to be rewritten.
See [reusable radiometric modules](radiometric.md) for a standalone build and
examples with or without CSPICE. Providers and delay callbacks
must outlive any operation referencing them. The existing EKF and its
simultaneous-epoch radiometric model are unchanged: this demo is a separate
observation-model baseline, not an EKF integration.

## Options and examples

```sh
./build-clang/cassini_two_way_demo --help
./build-clang/cassini_two_way_demo --samples 1 --ramp-hz-s 0.5
./build-clang/cassini_two_way_demo --samples 1 --station DSS-14 \
  --transponder-delay-us 0 --zenith-delay-m 0 --no-shapiro
./build-clang/cassini_two_way_demo --start-utc 2004-10-11T18:00:00 \
  --samples 6 --step-s 300 --count-s 60
```

No plots or data files are produced. Console output includes all four event
epochs, elevation, both leg light times, round-trip delay, coordinate range,
Doppler count, backend differences and the corrections' range/Doppler impact.

The three default counts produced maximum backend differences of approximately
0.91 picoseconds per leg and 0.00026 Hz in count Doppler on the development
machine. These are numerical agreement figures, not physical measurement
accuracy. The integration benchmark rejects differences above 10 ns or 0.005 Hz.

Analytic tests cover stationary range, transponder delay, zero stationary
Doppler, a ramp evaluated at transmission, constant clock bias cancellation,
receding two-way Doppler at small and large epochs, additive corrections on
both legs, radial Shapiro, zenith-path mapping and invalid inputs.

## Remaining fidelity inputs

These are synthetic predictions from reconstructed ephemerides, not real Cassini
tracking records. A calibrated DSN observation model additionally requires
actual transmitter ramp tables and count metadata, spacecraft and station group/
phase delay calibration, antenna phase centers and axis offsets, measured wet/dry
meteorology, ionosphere and solar plasma, Earth/Saturn/ring occultation checks,
clock errors and proper-time effects, and full relativistic frequency transfer.
Visibility here checks ground elevation at the count endpoints; it does not prove
visibility throughout a long count, spacecraft antenna pointing or availability
of a real tracking pass. Solar Shapiro alone does not account for Saturn delay.

## Longer arc with an EKF

```sh
./build-clang/cassini_two_way_ekf_demo
MPLCONFIGDIR=/tmp/deepnav-matplotlib .venv-cassini/bin/python tests/plot_cassini_two_way_ekf.py
```

The default run spans 24 hours at 600-second cadence, with 60-second counts.
Truth is Cassini's reconstructed SPK (`--truth spice`); synthetic Gaussian noise
is added with a declared seed. These are not recorded tracking data. The EKF
estimates six Saturn-relative Cartesian components and propagates a separate
trajectory using Saturn point mass plus J2, solar tides and eight major-moon
tides. GMs, Saturn radius and pole come from the kernels. Rounded J2=0.016298 is
from the [NASA Saturn fact sheet](https://nssdc.gsfc.nasa.gov/planetary/factsheet/saturnfact.html).
The pole is fixed at the arc start; higher harmonics, SRP and maneuvers are omitted.

Range sigma is 10 m; counted Doppler sigma is 0.01 Hz. Orbital SNC density defaults
to `1e-5 m/s^(3/2)` and is configured separately from measurement R with
`--acceleration-noise-density`. An optional `--truth matched` mode propagates
truth with the same force model, useful for isolating numerical/filter behavior.
The deterministic matched smoke test sets SNC to zero.

The driver's state epoch is 5000 TDB seconds before Earth reception, keeping both
spacecraft count endpoints after the reference state on this Cassini arc. The
trajectory factory propagates the current estimate forward and evaluates the
independent retarded legs. It never resets an estimate to the SPK truth. State
propagation, STM and Q use the same acceleration model. The receiver-state reuse
and LSK snapshot optimizations apply equally to truth and estimate calculations.
A fixed lookback is an arc-specific driver convention, not an assumption in the
reusable radio interfaces. Q in the short forecast from the state epoch to the
measurement epochs is not separately added to R; this remains a simplified
latency/noise treatment rather than a complete operational delayed-data filter.

`diagnostics.csv` records prefit innovations and NIS only at visible updates;
gaps contain NaNs for those fields. `summary.json` declares the arc and noise
settings. The plotting script creates only `ekf_rms_nis.png`, with two panels:

- Running temporal RMS of the **3D post-update position error**, including
  propagation-only epochs: `sqrt(mean(||estimated_position - truth_position||^2))`.
  This is one realization, not a Monte Carlo ensemble RMS or a covariance bound.
- Raw prefit NIS for each range/Doppler pair, with reference mean 2 and the
  chi-square(2) single-update 95th percentile, 5.991. No NIS gating is enabled.

The initial 24-hour SPK-truth run had 145 epochs and 130 visible updates (260
scalar observations), mean NIS about 1.36 and final running position RMS about
3.07 km. The 3D position error increased despite small innovations. NIS therefore
does not establish full-state recovery; weak geometry and force-model mismatch
need further investigation. This figure is a diagnostic result, not a validated
navigation-accuracy claim.
