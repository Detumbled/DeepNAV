# Reusable coherent two-way radiometric model

The numerical radio model is independent of mission, kernel source, orbit
propagator and estimator. Cassini is one application of it. The existing Voyager
example is scheduled for replacement; it is not a dependency or foundation of
these modules. Its replacement should consume the reusable interfaces rather
than contain another copy of the observation physics.

## Build without CSPICE

All targets are defined in the root `CMakeLists.txt`. Use the same `build-clang/`
directory for either configuration; no auxiliary CMake or build directories are
required. To build only the numerical model without configuring CSPICE or graphics:

```sh
cmake -S . -B build-clang -DCMAKE_BUILD_TYPE=Release \
  -DDEEPNAV_RADIOMETRIC_ONLY=ON -DDEEPNAV_RADIOMETRIC_WITH_CSPICE=OFF
cmake --build build-clang -j 4
ctest --test-dir build-clang --output-on-failure
```

`deepnav_radiometric` requires Eigen and C++23. An external project can use
DeepNAV through `add_subdirectory` with the numerical-only option enabled and
link `deepnav_radiometric`; there is no separate CMake module to include.

## Restore the full build with CSPICE

```sh
CSPICE_HOME=/path/to/cspice cmake -S . -B build-clang -DCMAKE_BUILD_TYPE=Release \
  -DDEEPNAV_RADIOMETRIC_ONLY=OFF -DDEEPNAV_RADIOMETRIC_WITH_CSPICE=ON
cmake --build build-clang --target cassini_two_way_demo cassini_two_way_ekf_demo -j 4
```

`deepnav_radiometric_spice` adds CSPICE without graphics linkage. Kernel loading
belongs to the application. `SpiceKernelSet` owns a caller-selected ordered list;
`SpiceTdbClock` snapshots LSK constants on construction and performs no pool reads
while evaluating. Recreate the clock after changing the LSK. There is no static
cache of the mutable global kernel pool.

The library light-time backend evaluates its fixed receiver once per leg and
reuses the existing solver's converged emitter state. Corrections iterate only
the retarded emitter. This reuse is local to one solve, so EKF state changes do
not require invalidating a persistent trajectory cache.

`TwoWayMeasurement` supplies range and counted Doppler plus a central-difference
Jacobian for an arbitrary trajectory factory. `IntegratedStateProvider` supplies
numerically propagated, Hermite-interpolated states and optionally adds a moving
origin to expose a fixed barycentric origin. It belongs to `deepnav_ekf`, retains
no SPICE dependency, and accepts any acceleration callback.

## Numerical provider example

```cpp
#include "dynamics/LinearStateProvider.hpp"
#include "observations/radiometric/TwoWayLink.hpp"

namespace radio = fd::observations::radiometric;
fd::dynamics::LinearStateProvider ground(
    fd::opnav::TdbEpoch{0}, {{0, 0, 0}, {0, 0, 0}}, "J2000", "SSB");
fd::dynamics::LinearStateProvider probe(
    fd::opnav::TdbEpoch{0}, {{3e6, 0, 0}, {10, 0, 0}}, "J2000", "SSB");

// Mission values are mandatory: no implicit Cassini ratio or uplink frequency.
radio::LinkConfig transponder{1.25, 0};
radio::FrequencyRamp ramp{0, 2e9, 0};
radio::TwoWayLink link(radio::libraryReceptionSolver(), transponder);
auto range = link.solve(100, ground, probe);
auto count = link.count(100, 60, ground, probe, ramp,
                        transponder.turnaroundRatio * ramp.frequencyHz);
```

The numerical tests compile `LinearStateProvider.cpp` explicitly. A consumer of
this example must similarly compile that implementation or provide its own
`StateProvider`. The reusable target does not choose a state propagation model.

A mission can supply `StateProvider` implementations backed by a numerical
integrator, interpolated ephemerides, analytic trajectories, an estimated EKF
state or an SPK reader. Endpoint states must have the same inertial frame and
fixed origin, and cover retarded epochs as well as count endpoints. For a real
Earth–spacecraft link, use barycentric states rather than subtracting positions
with an Earth-centered origin evaluated at different times.

## CSPICE provider example

```cpp
#include "dynamics/SpiceStateProvider.hpp"
#include "observations/radiometric/SpiceLightTime.hpp"

// The caller loads the mission and ground-station kernels first.
fd::dynamics::SpiceStateProvider ground("DSS-43");
fd::dynamics::SpiceStateProvider probe("-12345"); // Replace with the mission ID.
radio::TwoWayLink link(radio::spiceReceptionLeg, transponder);
auto range = link.solve(receptionTdb, ground, probe);
```

`spiceReceptionLeg` requires SPICE-backed endpoints because its CN query uses
their target IDs. `libraryReceptionSolver` accepts any `StateProvider`, including
SPICE-backed providers or a mixture of SPICE and numerical providers, as long as
both endpoints share the same frame and origin. This allows the same two-way
model to generate truth from CSPICE and predictions from a propagated estimate.

## Extension boundaries

- `StateProvider`: geometric state source. No radio parameters or kernel loading.
- `LegSolver`: reception light time with solved states and convergence diagnostics.
- `LinkConfig`: explicit transponder ratio and effective non-dispersive delay.
- `FrequencyRamp`: explicit carrier, slope and reference epoch.
- `DelayModel`: additive delay callback receiving the solved leg and both endpoint
  providers. Compose atmosphere, gravitational terms and calibrated path delays.
- `ClockOffset`: TDB minus ground clock reading. Omit for a coordinate-time model,
  or supply a deterministic/stochastic clock history. The solver itself does not
  advance a clock simulator or draw noise.

The Shapiro adapter uses the shared primitive in `perturbations/Shapiro.hpp`;
that module is independent of the Voyager file. Cassini's body choice and GM
come from the caller. The radio model does not assume that the Sun is the only
source of gravitational delay.

Group and carrier-phase observables may need different propagation/transponder
calibrations, especially for plasma. Use separate link configurations and delay
callbacks. The current link is coherent and returns to the same station; it does
not yet implement three-way links, a free-running one-way onboard oscillator,
piecewise ramp records or complete proper-time frequency transfer.
