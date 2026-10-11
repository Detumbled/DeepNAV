# Cassini SPK versus open-loop dynamics

The default arc is **9–23 October 2004**, ending before OTM-4. Both initial
position and velocity come directly from reconstructed SPK, with geometric
`NONE` states in Saturn-centered J2000. No initial-state fitting or later state resets occur.
The retained recoil vector is calibrated on the first seven days, with bounded
component changes; see [RECOIL_UPDATE.md](RECOIL_UPDATE.md).
The plots show the SPK and one retained model: multibody gravity, direct Saturn
J2, uniform rings, attitude-dependent HGA SRP and thermal recoil. Indirect J2 and
point-mass 1PN are disabled after a controlled eight-combination SPK comparison.
The comparison table and trajectory CSV retain only this model and the SPK
reference. A finer-step integration is stored separately in `validation.csv`.

```sh
cmake --preset apple-clang
cmake --build build-clang --target compare_cassini_spk test_perturbation_models -j 4
build-clang/test_perturbation_models
.venv-cassini/bin/python benchmarks/cassini_trajectory/audit_gravity_mass.py
.venv-cassini/bin/python benchmarks/cassini_trajectory/audit_attitude.py
build-clang/compare_cassini_spk
MPLCONFIGDIR=/tmp/deepnav-matplotlib .venv-cassini/bin/python benchmarks/cassini_trajectory/plot_comparison.py
```

Python requires NumPy, Matplotlib and SpiceyPy; Plotly produces the standalone
rotatable HTML. Results are in ignored `results/cassini_trajectory/`: CSV states,
comparison table, metadata/hashes, PNG/PDF plots and `trajectory_3d.html`.
The gallery in `COMPARISON.md` includes:

- `diagnostic_overview`: position/velocity error norms, RTN errors and force magnitudes.
- `trajectory_projections`: XY, XZ and YZ SPK/full-model overlays.
- `position_error_components`: signed J2000 XYZ errors and norm in metres.
- `velocity_error_components`: signed J2000 XYZ errors and norm in mm/s.
- `rtn_position_errors`: radial, along-track and cross-track position errors in metres.
- `perturbation_accelerations`: SRP/thermal magnitudes in nm/s² and sampled attitude sources.
- `gravity_additions`: indirect J2, ring gravity and relative 1PN magnitudes along SPK.
- `orbital_plane_residual`: radial versus along-track error, colored by elapsed day.

RTN axes use the instantaneous SPK position and angular momentum. The plotted
components are projections of the Cartesian position difference into that frame.
The force curves show separate magnitudes, not a vector sum; the attitude strip
represents output samples, while the manifest gives exact interval coverage.

## Forces and conventions

[Thermal recoil and SRP source: Di Benedetto, Iess & Roth (2009)](https://issfd.org/ISSFD_2009/OrbitDeterminationI/DiBENEDETTO.pdf).

`include/perturbations/ThermalRecoil.hpp` implements Eq. (4), with caller-supplied
body-to-integration rotation and mass, signed body acceleration in km/s²,
inverse mass scaling and Pu-238 half-life 87.74 Julian years. Epochs use seconds
in a common convention. Eclipse does not suppress internally generated recoil.
The generic header has no SPICE dependency. Prescribed attitude and mass imply
zero Cartesian state Jacobians.

The original recoil hypothesis uses body X/Y zero and the paper's **later-tour
mean -5.45e-12 km/s² at T7 mass 2781 kg**, referenced by convention to
13 August 2005. The retained reference vector is now
**[-2.725e-13, -2.725e-13, -5.232237415805135e-12] km/s²** in body X/Y/Z.
It is a bounded first-seven-day SPK calibration, not an October measurement;
X/Y reach the negative search limits. The second half is propagated without
state resets, but was used in choosing the lowest full-arc RMS candidate.
See [the sensitivity and prediction audit](RECOIL_UPDATE.md).
`--thermal-z-km-s2` changes the signed reference Z
coefficient; zero disables recoil. Current mass defaults to **3220 kg**, NASA PDS’s pre-burn OTM-004
mass on 23 October 2004; surrounding OTM-003/005 entries have the same rounded
value. Huygens was still attached. The table has about 10 kg reporting resolution
and maneuver snapshots, not continuous mass telemetry, so the benchmark uses a
constant estimate without interpolating through burns. Both SRP and thermal
recoil use this mass. `--mass-kg` supplies an alternative constant mass. This
default is justified for the selected October arc, not arbitrary mission dates.
The original NASA PDS CSV/README and source hashes are in `data/`;
`audit_gravity_mass.py` verifies the entries and records provenance.

`include/perturbations/FlatPlateSRP.hpp` implements a one-sided opaque Lambertian
plate with analytic position partials and optional existing Saturn eclipse
geometry. With q the Sun-to-spacecraft unit vector, n the outward illuminated
face normal, h=max(0,-n·q), its acceleration is

`a = P_1AU (AU/d)^2 A h / (1000 m) [q - 2 nu n]` (km/s²).

The paper's nu=(1-alpha)/3 is not the diffuse energy fraction. The default
nu=0.21 gives normal-incidence equivalent **CR=1.42**, from the paper's table
and Appendix Eqs. (6)–(8). The `--cr` option now parameterizes this equivalence
and must be between 1 and 5/3. The HGA face normal is body -Z; its physical area
defaults to 12 m² (paper approximation to a 4 m antenna). The projected illuminated
area is computed from attitude. No force is assigned to its back face. This
omits bus surfaces, antenna backside, specular reflection and the paraboloid's
local normals, so it is a reduced HGA model rather than the full flight SRP model.

The point-mass and third-body terms retain the Sun, eight other planetary-system
barycenters, and Mimas, Enceladus, Tethys, Dione, Rhea, Titan, Hyperion, Iapetus
and Phoebe, plus Helene (the tenth nonzero moon mass in the release label).
Indirect Saturn-origin acceleration is included. `Sat427Gravity.hpp` reads
Saturn/moon GMs, Jupiter/Uranus/Neptune system GMs, J2=0.01629054544344128,
gravity reference radius **60330 km**, and the three-term analytic pole series
from the [SAT427.23 release label paired with this SPK](https://naif.jpl.nasa.gov/pub/naif/CASSINI/kernels/spk/200128RU_SCPSE_04251_04303.bsp.lbl).
RA uses sine and declination cosine; time is TDB Julian centuries from POLTIM.
The formula convention is [Jacobson (2022), Eqs. 1–2](https://doi.org/10.3847/1538-3881/ac90c9);
that paper’s later SAT441 coefficients are not substituted. Download the linked
`.bsp.lbl` next to its SPK binary if it is missing; the benchmark requires it. The independent
audit compares the pole with CSPICE at six epochs, including ±100 years.
This uses the release’s analytic approximation rather than a numerical pole
integration. Saturn’s physical shadow radius remains **60268 km**.

SAT427’s GM10 lumps the inner planets into the Sun. Since this benchmark
resolves planetary barycenters 1–4, their generic-kernel GMs are subtracted
before using GM10. This is a resolved approximation to the lumped model.
Other planetary GMs retain the generic DE440 kernel; paired SPK ephemerides
have load priority. Higher harmonics remain omitted, as requested.
The J2 origin reaction, separate ring gravity and barycentric point-mass 1PN
relativity are implemented and tested. The retained default enables rings only
from these additions (mask 2), keeping direct J2. Moon shadowing and maneuver accelerations remain
absent. See [GRAVITY_UPDATE.md](GRAVITY_UPDATE.md) for equations, assumptions,
validation and the measured before/after comparison.

## Attitude kernels

Five reconstructed CKs and two predicted CKs generated from **as-flown sequences**
are downloaded from [NASA NAIF](https://naif.jpl.nasa.gov/pub/naif/CASSINI/kernels/ck/).
`attitude_manifest.json` records exact filenames, URLs, SHA-256 hashes, types,
coverage durations and gaps. The binaries remain in ignored `Kernels/Cassini/ck`.
For a missing binary, download the URL in the manifest to its listed local path.
The existing `cas_v40.tf`, `cas00172.tsc` and `naif0012.tls` resolve the
`CASSINI_SC_COORD` frame (-82000). The as-flown CK comments explicitly identify
`cas00172.tsc` as their clock correlation kernel.

As-flown files load first and reconstructed files last, so measured telemetry
has priority wherever available. Reconstructed coverage is about 93.96% of
this arc; the remaining approximately 20.3 hours use as-flown sequence attitudes.
Those are predicted orientations, not reconstructed telemetry. A 6.309023 s
boundary gap on 18 October is bridged by quaternion SLERP between valid endpoint
rotations. The code rejects gaps longer than 10 s and any extrapolation beyond
coverage. The audit uses `ckcov` INTERVAL windows and validates frame resolution
at every covered minute and at the bridge endpoints. CSV rows identify the
attitude source at each output epoch. Kernels are queried at every integration
stage, not just at output times.

## Numerical checks and interpretation

DP853 uses nominal maximum step 120 s and absolute tolerance 1e-12 for the
attitude-based models. The refined full-force run uses 30 s and 1e-13; relative
tolerances are 1e-12 and 1e-14. A separate refined integration checks the combined trajectory. The plotter
checks finite values, exact initial state, cadence and error norms, and requires
nominal/refined agreement within 0.1 m. Numerical agreement validates integration
of this specified model, not its optical, mass or thermal assumptions.

Use `--start-utc`, `--days`, `--cadence-s`, `--area-m2`, `--mass-kg`, `--cr` and
`--thermal-x-km-s2`, `--thermal-y-km-s2` and `--thermal-z-km-s2` to vary inputs.
Signed recoil coefficients may be zero. The gallery requires a refined numerical
validation run of the same physical model. Run the attitude audit with matching
start and duration before plotting a different arc; insufficient CK coverage
fails explicitly.

The mass table lists OTM-004 at 06:16 UTC on 23 October; the
[weekly NASA report gives 07:29:28 UTC](https://science.nasa.gov/missions/cassini/significant-event-report-for-week-ending-10292004/).
Both times follow the endpoint, by at least 6 h 16 min; the benchmark uses
the earlier time conservatively. Absence of every attitude-control thruster
firing has not been established.


## Requested ±100 kg checks

A 3320 kg sensitivity run, with the selected gravity model, original nominal recoil and all other inputs
fixed, gives 110.281 m position RMS and 258.976 m endpoint error. Both exceed
the retained 3220 kg result (109.133 m / 249.203 m), so 3220 kg remains the
default and the sole retained simulation. The mass increase scales SRP and
thermal recoil down by 3.012%; it is not a revised PDS measurement. The full
comparison, RTN changes and numerical checks are recorded in
[GRAVITY_UPDATE.md](GRAVITY_UPDATE.md). The alternative output is removed.

The subsequent 3120 kg test increases SRP and thermal recoil by 3.205% and gives
119.380 m RMS and 263.148 m endpoint error. Its endpoint along-track residual
improves from -241.190 to -227.698 m, but cross-track grows from -55.187 to
-131.542 m. Both overall metrics worsen. Of the three masses tested, 3220 kg
remained closest to SPK in that mass study; the lower-mass
output is also removed. Its finer-step check agrees within 6.130 mm.

## Retained configuration and next diagnostics

Only **SAT427 + PDS 3220 kg + bounded XYZ recoil calibration** is retained.
Its full-arc RMS is **78.174 m**, versus 109.133 m for the original recoil, and
its endpoint error is **196.992 m**, versus 249.203 m. Second-half RMS improves
from 151.673 to 109.688 m. The negative-X-only trial gives a smaller endpoint
error but higher RMS, so it is not retained. The nine-candidate audit and exact
coefficients are in [RECOIL_UPDATE.md](RECOIL_UPDATE.md).
The alternate simulations from the current recoil study and its temporary driver
are removed after validation. `validation.csv` is the same physical model with finer numerical
settings, not a second force-model scenario. Among the three new gravity terms,
only rings are enabled in the retained model. The eight temporary gravity
ablation runs are removed after validation; the audit table is retained in
[GRAVITY_UPDATE.md](GRAVITY_UPDATE.md).

Beyond gravity, the last update replaced the approximate 2700 kg mass with the
rounded NASA PDS October estimate of 3220 kg. This reduces both SRP and thermal
accelerations by 16.15% at unchanged geometry and recoil normalization. The
nominal thermal and SRP coefficients were retained without fitting at that stage. Enabling all
three new gravity terms subsequently worsened position RMS to 164.327 m.
One-at-a-time and all-combination tests identified indirect J2 as the dominant
driver of that increase: it raises RMS by 52–55 m in every matched pair. Rings
reduce RMS in every pair; 1PN raises it slightly. With the original nominal recoil, the rings-only addition
gave **109.133 m RMS** and **249.203 m endpoint error**, versus the fresh baseline's
113.276 m and 257.951 m. Endpoint RTN is -29.7 m radial, -241.2 m along-track and
-55.2 m cross-track. The finer-integration check agrees within **3.470 mm**.
These results identify the recent worsening, not the original drift's cause.
The selection favors this SPK arc's agreement and does not establish that the
disabled J2 origin correction is physically unnecessary.

A preliminary review of the two NASA PDS small-force files covering this
period found no predicted RCS delta-V record inside the current arc. The nearest
preceding record ends on 6 October at 10:10:15.500 UTC; the next starts on
23 October at 05:16:51.250 UTC. These files label their records as predicted,
so actual unplanned firings are not ruled out. A 7–23 October **16-day** arc
is a candidate extension: existing CK coverage (including the same 6.309 s
bridge) and SPK endpoint availability have been checked, without changing the
retained simulation. See [ARC_REVIEW.md](ARC_REVIEW.md) for sources and details.

Saturn's oblateness changes the acceleration of the moving origin through its
interaction with the moons. The correction in JPL §4.4.6 is implemented but
disabled in the selected default. The SAT427 ring GMs sum to 1.123 km³/s²; this exactly
accounts for the difference between the label's system GM and the sum of its
Saturn/moon masses. Their field is approximated by uniform concentric annuli.
The available relativity implementation uses Moyer Eq. 4-26 with barycentric
velocities and subtracts Saturn's 1PN correction from Cassini's; it is also
disabled in the selected default. No optical/recoil coefficients were changed.

Next, improve the physical SRP geometry for the Cassini/Huygens configuration:
illuminated bus, antenna backside and boom surfaces, with separate optical
coefficients. The current one-sided HGA plate cannot represent these forces
when attitude exposes other surfaces. The paper describes optical evolution
and non-gravitational force estimation, so lowering CR alone is not an
established correction. The transferred later-tour recoil coefficient is also
not an October measurement, and body X/Y may be nonzero or affected by other
thermal emissions. [Di Benedetto, Iess & Roth (2009)](https://issfd.org/ISSFD_2009/OrbitDeterminationI/DiBENEDETTO.pdf).

After checking impulses and surface geometry, a constrained estimation of the
three body-frame recoil components could test a remaining smooth force bias.
Keep the SPK initial state fixed, estimate on an initial subarc, and validate
prediction on the withheld remainder or another quiet arc. A coefficient fit
against this reconstructed SPK would be a surrogate calibration, not independent
validation from tracking data. Recoil/SRP parameter correlation must be checked.
The paper's fitted multi-arc/stochastic navigation solutions should not be
interpreted as an expectation of zero error in a 14-day open-loop propagation.
Gravity physics still has approximations, including uniform ring density,
neglected relativistic spin/multipoles and omitted higher Saturn harmonics;
alignment of constants alone does not rule these out. Further
refinement of numerical tolerances is low priority given the millimetre check.
