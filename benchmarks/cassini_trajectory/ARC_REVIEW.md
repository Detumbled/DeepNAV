# Arc extension and model review

The retained simulation remains 9–23 October 2004, SAT427 + 3220 kg, now with
the bounded XYZ recoil calibration described in [RECOIL_UPDATE.md](RECOIL_UPDATE.md).
This arc-extension review does not generate another simulation.

## Events near the arc (UTC)

| Record | Start | Stop | Predicted net delta-V |
|---|---|---|---:|
| AACS small-force prediction | 2004-10-02 12:11:33.125 | 12:23:29.875 | 3.050 mm/s |
| AACS small-force prediction | 2004-10-06 10:00:33.125 | 10:10:15.500 | 2.184 mm/s |
| AACS small-force prediction | 2004-10-23 05:16:51.250 | 05:28:59.250 | 2.117 mm/s |
| AACS small-force prediction | 2004-10-23 06:13:10.000 | 06:23:36.375 | 1.730 mm/s |

These are norms of J2000 delta-V vectors in two NASA PDS files. Their headers
explicitly label records `P` (predicted), not measured firings. The second
23 October interval overlaps the maneuver-history OTM-004 window; event
identity must be reconciled before adding both as independent model inputs.

Sources:

- [04256_04292pd_S4_psiv2.DV](https://pds-atmospheres.nmsu.edu/data_and_services/atmospheres_data/Cassini/Cassini/DOCS%20for%20events%20%26%20Configuration%20page/Small%20Forces/04256_04292pd_S4_psiv2.DV)
- [04292_04320pe_fsiv2.DV](https://pds-atmospheres.nmsu.edu/data_and_services/atmospheres_data/Cassini/Cassini/DOCS%20for%20events%20%26%20Configuration%20page/Small%20Forces/04292_04320pe_fsiv2.DV)
- [PDS engineering overview](https://pds-atmospheres.nmsu.edu/data_and_services/atmospheres_data/Cassini/eng-overview.html)

The maneuver performance table lists OTM-004 on 23 October at 06:16 UTC,
464.375 s RCS burn, and OTM-005 on 29 October at 06:15 UTC, 4.600 s main-engine
burn. The weekly NASA report gives a different OTM-004 start, 07:29:28 UTC.
The earlier predicted small-force event at 05:16 is the conservative forward
boundary regardless of that discrepancy. The following predictor file also
lists three RCS intervals on 26 October (13:41–16:57 UTC), around Titan-A;
extending beyond 23 October requires explicit maneuver treatment.

[Burn-performance table](https://pds-atmospheres.nmsu.edu/data_and_services/atmospheres_data/Cassini/Cassini/DOCS%20for%20events%20%26%20Configuration%20page/cassini_mea_rcs_burns.pdf),
[table README](https://pds-atmospheres.nmsu.edu/data_and_services/atmospheres_data/Cassini/Cassini/DOCS%20for%20events%20%26%20Configuration%20page/README_MEA_RCS_BURNS.txt),
[NASA weekly report](https://science.nasa.gov/missions/cassini/significant-event-report-for-week-ending-10292004/).

## Candidate extension

**7 October 00:00 UTC to 23 October 00:00 UTC (16 days)** adds two days
backward and remains between the neighboring predicted events. The relevant
files contain no predicted event in this interval. The existing attitude audit
passes on this interval: reconstructed coverage 1298846.947 s, predicted
as-flown coverage 83546.744 s, same 6.309 s short bridge. Geometric SPK
availability was checked for Cassini and every modeled perturber at both
endpoints. No extended trajectory has been propagated yet. Any extended run
must initialize all six state components from SPK at its new start.

## Priorities before tuning parameters

Update: priorities 1–2 have now been implemented and reintegrated on the original
14-day arc. Subsequent eight-combination ablation retains only rings from those
additions, with direct J2 still enabled; indirect J2 and 1PN are disabled in the
SPK-match-selected default. See [GRAVITY_UPDATE.md](GRAVITY_UPDATE.md). The following rationale
records the original diagnosis; the proposed 16-day extension remains unrun.

1. **Complete Saturn-centered gravity bookkeeping.** The code adds direct
   Saturn J2 to Cassini and subtracts point-mass Saturn-origin acceleration.
   It previously omitted the acceleration of Saturn due to its own oblateness interacting
   with the moons. JPL §4.4.6 treats this explicitly for outer planet systems.
   For point-mass moons, the additional relative-acceleration term is
   `sum((GM_moon/GM_Saturn) * a_J2(r_moon))`, with moon positions relative to
   Saturn and the same inertial pole. This is a missing indirect term, not
   a change to J2 or an addition of higher harmonics. A sampled force audit
   on the current arc gives a magnitude of about 0.026–1.300 nm/s².
   A simple double time integral gives a 116 m displacement scale; it excludes
   gravitational feedback and is not a predicted reduction of the residual.
   A correct reintegration may improve or worsen the SPK discrepancy.
   The paired label's ring GMs total 1.123 km³/s²; their monopole and shape fields
   are now included using uniform concentric annuli.
   [JPL Moyer formulation, §4.4.6](https://descanso.jpl.nasa.gov/monograph/series2/Descanso2_all.pdf).
2. **Match the relativistic dynamics.** The Di Benedetto paper includes
   Saturn/Jupiter/Sun relativistic perturbations and propagates in a
   Saturn-system barycentric frame. Our TDB, Saturn-centered model must use
   appropriate relative equations, not append a heliocentric formula blindly.
   [Paper §1.2](https://issfd.org/ISSFD_2009/OrbitDeterminationI/DiBENEDETTO.pdf).
3. **Improve SRP geometry and October optical coefficients.** Include the bus,
   antenna backside, boom and attached Huygens where illuminated. The current
   one-sided HGA plate assigns no force to these surfaces. The paper documents
   changing HGA optical properties and uses temperature measurements to
   constrain them; its tracking-oriented approximation is less representative
   of arbitrary science attitudes. Keep optical data and frame conventions
   consistent rather than lower CR solely to reduce this arc's error.
4. **Then constrain all three body-frame thermal components.** Our coefficient
   is a transferred later-tour mean and X/Y are assumed zero. Estimate limited
   adjustments on an initial subarc with the SPK initial state fixed; test
   prediction on withheld epochs or another quiet arc. SRP and Z recoil can
   be strongly correlated. A fit against SPK is a surrogate calibration,
   not independent validation with Doppler/range. The paper uses estimated
   forces, multi-arc states and/or stochastic accelerations, so its accuracy
   does not imply a 14-day uncorrected propagation should match exactly.

Priorities 2–4 are informed by [Di Benedetto, Iess & Roth (2009), §§1.2–2 and
Appendix](https://issfd.org/ISSFD_2009/OrbitDeterminationI/DiBENEDETTO.pdf).
The force implementations for priorities 1–2 remain available; priority 3
remains a proposal. A bounded first-seven-day implementation of priority 4
reduces full-arc RMS to 78.174 m, with no state resets. The remaining half was
excluded from coefficient estimation but used in candidate selection; another
arc is needed for independent validation. See [RECOIL_UPDATE.md](RECOIL_UPDATE.md).
Ablation with the original recoil identified indirect J2 as the dominant cause of the
recent RMS increase, adding 52–55 m RMS across matched configurations. The
rings-only addition gave 109.133 m RMS and 249.203 m endpoint error before recoil calibration.
Selection against this arc does not establish physical validity of omitting an
origin correction, or explain the original along-track drift. Numerical
tolerances remain low priority given the millimetre convergence checks.
