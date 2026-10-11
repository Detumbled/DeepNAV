# Bounded thermal recoil study

**The retained configuration uses 3220 kg and the bounded XYZ recoil calibration.**
The 9–23 October 2004 arc, exact SPK initial position/velocity, SAT427 gravity mask 2,
attitude, SRP geometry/optical coefficients and numerical settings are unchanged.
Only the signed body-frame recoil reference vector is adjusted.

## Candidates and selection

| Recoil candidate | First 7 days RMS (m) | Second 7 days RMS (m) | Full RMS (m) | Endpoint (m) |
|---|---:|---:|---:|---:|
| Nominal Z-only | 28.742 | 151.673 | 109.133 | 249.203 |
| X negative, 5% | 25.623 | 118.787 | 85.907 | 192.427 |
| X positive, 5% | 34.876 | 195.738 | 140.555 | 320.496 |
| Y negative, 5% | 23.006 | 148.441 | 106.191 | 249.867 |
| Y positive, 5% | 36.964 | 176.468 | 127.461 | 283.937 |
| Z magnitude +5% | 36.199 | 179.413 | 129.391 | 280.307 |
| Z magnitude -5% | 28.138 | 165.558 | 118.718 | 281.229 |
| Z fit, first 7 days | 27.410 | 154.545 | 110.960 | 260.399 |
| XYZ fit, first 7 days — retained | 14.023 | 109.688 | 78.174 | 196.992 |

The full-arc RMS decreases by 28.37% (109.133 → 78.174 m).
Second-half RMS decreases by 27.68% without resetting the state.
The negative-X trial has a smaller endpoint error than the selected XYZ vector,
but XYZ has the lowest whole-arc RMS, the user-requested selection criterion.

## Calibration and limits

The original reference vector is [0, 0, -5.45e-12] km/s² at 2781 kg.
Each trial changes one component by ±2.725e-13 km/s² (5% of the original Z magnitude).
The two estimates use only samples at elapsed times 0–604800 s inclusive
(1009 positions); the remaining 1008 samples are excluded from coefficient estimation.
Cartesian residuals receive equal weight in all three directions at every 600 s epoch.
Initial position and velocity are fixed, with no state updates or resets.

For each axis, central differences of the two fully integrated sensitivity trajectories
give a position-response column per unit bounded change. A bounded linear least-squares
solution minimizes the training residual vector with normalized component changes
between -1 and +1. The Z-only solve uses one column; the XYZ solve uses all three.
Both estimated candidates are then independently reintegrated with the nonlinear force
model through all 14 days. They are not synthetic sums of trajectories.

The retained reference vector, in body X/Y/Z, is:

```text
[-2.725e-13, -2.725e-13, -5.232237415805135e-12] km/s²
```

Reference mass 2781 kg, epoch 2005-08-13T00:00:00 and Pu-238 half-life 87.74 years
remain unchanged; the actual force uses 3220 kg and the prescribed CK attitude.
Z magnitude decreases by 3.996%; nonzero X/Y components change the recoil direction.
X and Y reach the negative limits of the allowed search. This is the best tested
configuration, not a demonstrated continuous optimum or an October thermal measurement.

Endpoint RTN becomes [-112.456, -160.893, 16.526] m.
Radial error worsens while along-track and cross-track errors improve; the vector
magnitude and full-arc RMS improve. Z-only tuning slightly improves training RMS but
worsens the later prediction, so it is not retained.

The second half was excluded from estimating coefficients, but its residuals were used
to select the lowest full-arc RMS candidate as requested. It is therefore a prediction
check, not an untouched final validation set. A separate quiet arc is still needed.
The calibration could absorb other force errors, including incomplete SRP geometry
and the disabled indirect J2 correction. It does not identify thermal recoil uniquely.

The paper estimates body-frame recoil components and discusses uncertain SRP optical
coefficients and variable attitude. Our procedure instead calibrates against SPK positions,
not independent tracking observations. [Di Benedetto et al., §§1.2–2](https://issfd.org/ISSFD_2009/OrbitDeterminationI/DiBENEDETTO.pdf).

## Retention and reproduction

The executable defaults to the selected vector. `build-clang/compare_cassini_spk`
reproduces the retained simulation; the plotting script regenerates its figures.
For the original nominal vector, use the diagnostic override:

```sh
build-clang/compare_cassini_spk --output /tmp/cassini-nominal-recoil --study-only --thermal-x-km-s2 0 --thermal-y-km-s2 0 --thermal-z-km-s2 -5.45e-12
```

Other trial coefficients, exact scores, fit results and sensitivity diagnostics are
preserved in [data/recoil_study.json](data/recoil_study.json). The selected model and
SPK reference are the only trajectories retained in `results/cassini_trajectory`;
`validation.csv` is the same model integrated with finer numerical settings.
Temporary candidate trajectories and study scripts are removed after validation.

The retained nominal/refined position difference is at most **1.384 mm**.
The nominal rerun gives **78.173772 m RMS**, reproducing the selected trial.
Force tests and plot/export validation pass.
