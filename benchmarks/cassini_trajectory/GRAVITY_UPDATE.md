# Gravity ablation and retained SPK comparison

The retained default is **multibody + direct Saturn J2 + rings + attitude SRP +
bounded XYZ thermal recoil calibration**, with indirect J2 and relativity disabled (gravity mask
2). The arc remains 9–23 October 2004, 14 TDB days, with the exact initial
geometric SPK position and velocity, SAT427 constants/pole, 3220 kg, the same
attitude kernels and unchanged SRP coefficients. The current recoil vector and
78.174 m RMS are documented in [RECOIL_UPDATE.md](RECOIL_UPDATE.md).
The gravity and mass tests below preceded that calibration and used the original
nominal Z-only recoil; no coefficients were fitted in those historical tests.

## Controlled one-at-a-time and combination tests

All eight combinations were reintegrated at the same nominal tolerances, 120 s
maximum step and 600 s export cadence. The choice minimizes the sampled
whole-arc position RMS; the same configuration also has the smallest endpoint
position error among the candidates. Direct J2 and point-mass origin terms are
always enabled. The mask switches only the three recent additions.

| Additions to baseline | Mask | Position RMS (m) | Endpoint error (m) | Endpoint velocity error (mm/s) |
|---|---:|---:|---:|---:|
| Baseline: direct J2, no new terms | 0 | 113.276 | 257.951 | 0.438660 |
| Indirect J2 only | 1 | 166.084 | 350.466 | 0.427073 |
| Rings only — retained | 2 | 109.133 | 249.203 | 0.413828 |
| Relativity only | 4 | 114.671 | 261.158 | 0.448754 |
| Indirect J2 + rings | 3 | 163.877 | 344.186 | 0.392654 |
| Indirect J2 + relativity | 5 | 166.891 | 352.821 | 0.439470 |
| Rings + relativity | 6 | 110.029 | 250.956 | 0.418175 |
| All three additions | 7 | 164.327 | 345.448 | 0.399475 |

**The indirect J2 addition is the dominant cause of the recent worsening.**
Enabling it alone increases RMS by 52.808 m and changes the endpoint along-track
residual from -251.854 m to -348.586 m, a further -96.733 m displacement.
Across all four matched pairs (holding the other two switches fixed), it increases
RMS by 52.221–54.744 m. Ring gravity reduces RMS in every matched pair by
2.208–4.642 m. Relativity increases it by 0.450–1.394 m. These comparisons isolate
the cause of the **increase after the gravity update**, not the cause of the
original nonzero SPK residual.

The retained rings-only addition gives **109.133 m RMS** and **249.203 m endpoint
error**, a 3.66% RMS improvement over the freshly reintegrated baseline. Endpoint
RTN is approximately -29.7 m radial, -241.2 m along-track and -55.2 m cross-track.
The new baseline differs from the previous 113.275 m value by about 1 mm, within
the integration-convergence scale after changing the force-summation order.

This is a configuration selected for agreement with this SPK arc. A required
origin correction does not become physically unnecessary because disabling it
reduces the residual. The smaller residual may reflect compensating model errors;
that possibility has not been established. The selection is not independent
validation on another arc or tracking data.

## Requested +100 kg mass sensitivity

The requested 3320 kg test retains gravity mask 2, the original arc and exact
SPK initial state, with unchanged attitude, area, optical coefficient and thermal
reference coefficient. Only the constant spacecraft mass changes from 3220 kg.
This scales SRP and thermal recoil by 3220/3320 = 0.96987952, a 3.012% reduction.
The gravity parameters are unchanged. Both runs include finer-step validation.

| Constant mass (kg) | Position RMS (m) | Endpoint error (m) | Endpoint velocity error (mm/s) | Nominal/refined maximum (mm) |
|---:|---:|---:|---:|---:|
| 3220 | 109.133 | 249.203 | 0.413828 | 3.470 |
| 3320 | 110.281 | 258.976 | 0.438013 | 1.493 |

**3220 kg remains retained.** The 3320 kg run has higher full-arc RMS and endpoint
error, despite improving the endpoint cross-track component. Endpoint RTN changes
from [-29.731, -241.190, -55.187] m to [-48.411, -253.871, +16.566] m.
The additional mass therefore does not resolve the dominant along-track drift.
Early improvement (at day 9.2) did not persist in the full-arc comparison.

3320 kg is a sensitivity scenario, not a revision to the PDS nominal mass. This
change scales both non-gravitational forces together and does not isolate an SRP
coefficient from thermal recoil. The prior eight-mask gravity ranking was done
at 3220 kg and was not repeated at 3320 kg. The alternative mass trajectory is
removed after the comparison; the table preserves the result.

## Requested -100 kg mass sensitivity

The 3120 kg test changes only mass, keeping mask 2, the arc, exact SPK initial
state and all force coefficients fixed. SRP and thermal recoil both increase
by 3220/3120 - 1 = **3.205%** relative to the nominal mass.

| Constant mass (kg) | Position RMS (m) | Endpoint error (m) | Endpoint velocity error (mm/s) | Nominal/refined maximum (mm) |
|---:|---:|---:|---:|---:|
| 3120 | 119.380 | 263.148 | 0.419077 | 6.130 |
| 3220 | 109.133 | 249.203 | 0.413828 | 3.470 |
| 3320 | 110.281 | 258.976 | 0.438013 | 1.493 |

**3220 kg remains the sole retained simulation.** The lower-mass endpoint RTN
is [-9.853, -227.698, -131.542] m, versus the nominal
[-29.731, -241.190, -55.187] m. Radial and along-track endpoint errors improve,
but cross-track error grows enough to worsen both total endpoint error and
whole-arc RMS. Decreasing mass therefore does not resolve the residual.
The numerical differences are millimetres, far below the observed changes.

Both mass alternatives are sensitivity scenarios, not revised PDS measurements.
They scale SRP and recoil together and cannot identify mass independently of
force-coefficient errors. Of these three tested masses, 3220 kg gives the
lowest RMS; this does not locate a continuous optimum. Gravity switches were
not re-ranked at either alternative mass. The alternative outputs are removed;
this table preserves the comparison.

## Retention and reproduction

Only the selected model and SPK reference are retained in
`results/cassini_trajectory`. `validation.csv` uses the same selected forces with
finer numerical settings. The historical nominal-recoil winner agreed within
**3.470 mm**; the current calibrated-vector convergence check is recorded in
[RECOIL_UPDATE.md](RECOIL_UPDATE.md). Force tests and the retained plot/export
checks pass. The eight trial
trajectories and temporary study driver have been removed. The compact
table above preserves the audit without retaining alternate simulations.

`compare_cassini_spk` defaults to mask 2. `--gravity-mask 0..7` allows reproducible
diagnostic overrides (bits 1=indirect J2, 2=rings, 4=relativity); `--study-only`
skips the finer integration for temporary experiments. The plotting script
requires numerical validation and derives labels/enabled-force statements from
metadata, rather than calling every combination a complete physical model.

## Origin correction

This implementation remains tested and available but is disabled in the retained default.

For moon position `r_i` relative to Saturn, the Saturn J2 field produces an equal
and opposite force on the central mass. Subtracting the acceleration of our
moving origin therefore **adds**

`sum_i (GM_i / GM_Saturn) * a_J2(r_i)`

to Cassini's Saturn-centered acceleration. This is the reaction/origin term in
[Moyer §4.4.6](https://descanso.jpl.nasa.gov/monograph/series2/Descanso2_all.pdf).
It uses all ten modeled nonzero-mass moons, their prescribed SPK positions, and
the same SAT427 pole, radius and J2 as the direct Cassini term. With concentric
co-accelerating rings, the denominator becomes the composite central GM. This
changes the J2 reaction by only about 3e-8 fractionally. Because moon ephemerides
are prescribed, the origin term has zero spacecraft position/velocity partials.

`j2OriginCorrection` is implemented in the perturbation library. This does not
change J2 or add Saturn J4/J6. Reactions from the distant Sun/planets through
Saturn's non-spherical field are still neglected.

## Rings

The [paired SAT427 label](https://naif.jpl.nasa.gov/pub/naif/CASSINI/kernels/spk/200128RU_SCPSE_04251_04303.bsp.lbl)
lists these separate masses and radii:

| Ring | GM (km³/s²) | Inner radius (km) | Outer radius (km) |
|---|---:|---:|---:|
| C | 0.039 | 74500 | 92000 |
| B | 0.784 | 92000 | 117500 |
| A | 0.300 | 122357 | 136780 |

Their total GM, 1.123 km³/s², agrees with the system GM minus the sum of the
label's Saturn/moon masses. We therefore include them separately without
replacing Saturn's GM with the entire system GM (which would duplicate moons).

Each ring is approximated as a thin, uniform-surface-density, concentric annulus
with Saturn's pole. Its field is evaluated by the exterior Legendre expansion of
the Newton integral through degree 80. The area-weighted radial moment for even
degree `l` is

`<s^l> = (Ro^(l+2) - Ri^(l+2)) / ((l/2+1) * (Ro²-Ri²))`.

The potential is `GM/r * sum_l P_l(0) <s^l>/r^l P_l(pole·r/r)` and the acceleration
is its spatial gradient. The monopole is included. This is a ring-mass expansion,
not extra Saturn gravity coefficients. The interface rejects points within
1.2 outer radii; all evaluated spacecraft and moon positions are exterior.
[Krogh, Ng & Snyder (1982)](https://ntrs.nasa.gov/citations/19820045807) describe
the underlying uniform-disk field and annuli as differences of disks.

Under the concentric co-acceleration assumption, the moving origin uses composite
mass `GM_Saturn + GM_rings`. The ring monopole reaction cancels the corresponding
mass-normalization change in the point-mass origin acceleration already present.
Only the ring shape field contributes an extra moon-driven origin correction.
Uniform density and a shared center are explicit approximations: the label does
not establish that this reproduces SAT427's full ring force implementation.

## Relativity

This implementation remains tested and available but is disabled in the retained default.

`RelativisticSystem` implements the **correction only** from Moyer Eq. 4-26,
with GR parameters beta=gamma=1 and c=299792.458 km/s. The 20 sources are Saturn
plus the existing Sun, planetary-system barycenters and ten moons. Source
Newtonian potentials and accelerations supply the first post-Newtonian cross
terms. Saturn is omitted from its own outer force sum, but its mass remains in
other sources' potentials and accelerations.

Velocities are solar-system barycentric. Positions can share an instantaneous
translation because the formula uses only their differences. The relative
correction is `a_1PN(Cassini) - a_1PN(Saturn)`. This preserves the Saturn-centered
TDB/J2000 formulation and adds no duplicate Newtonian gravity. The paper lists
Saturn, Jupiter and Sun relativistic perturbations; this implementation retains
all the modeled point masses. [Di Benedetto et al., §1.2](https://issfd.org/ISSFD_2009/OrbitDeterminationI/DiBENEDETTO.pdf).
Relativistic spin, ring and multipole terms remain omitted. The acceleration-only
ring and 1PN interfaces do not provide variational partials.

## Verification and force audit

The perturbation tests check the J2 reaction sign/mass weighting, zero limits,
zero spacecraft partials and rotation covariance; annulus axial and far-field
analytic limits, independent direct Newton quadrature and rotation covariance;
and the relativistic fixed-mass Schwarzschild and finite-mass two-body relative
1PN formulas, translation invariance and singular-input rejection.

`force_audit.csv` records signed vectors along the prescribed SPK path at 600 s
cadence. These are evaluated force magnitudes, not isolated trajectory changes:

| Addition | Minimum (nm/s²) | Maximum (nm/s²) |
|---|---:|---:|
| J2 origin correction | 0.025896 | 1.299947 |
| Rings, direct plus origin shape | 0.031591 | 0.158141 |
| Relative 1PN | 0.008387 | 0.038965 |

The generated comparison report includes the retained ring-force plot,
updated 3D trajectory, signed XYZ/RTN residuals and the original SRP/recoil plots.
The next useful physical checks remain October spacecraft geometry/optical
properties, attached Huygens, recoil transfer from the later tour, and small
attitude-control impulses. A smaller residual from parameter fitting would need
validation on withheld epochs or another quiet arc.
