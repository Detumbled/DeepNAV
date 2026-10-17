# Optical pipeline review against Owen (2024)

Source: William M. Owen, Jr., *Spacecraft Optical Navigation*, JPL, DESCANSO,
October 2024. Reviewed the user-supplied 136-page PDF at
`/Users/adb/Library/CloudStorage/GoogleDrive-alessandrodebellis01@gmail.com/Il mio Drive/Books/Space/Descanso/Spacecraft-Optical-Navigation.pdf`.
Page numbers below are the document's printed page numbers, not PDF viewer indices.
This is a review and proposed validation plan; no generator or EKF behavior was
changed during this review.

## Architecture assessment

The proposed sequence is technically sound as an initial architecture:

1. Simulated spacecraft truth + target/Sun ephemerides + declared truth camera
   attitude/calibration produce synthetic detector images.
2. Image processing estimates geometric centers and their uncertainties using
   a body brightness model.
3. An optical observation model predicts geometric-center pixels from the EKF
   state, target ephemeris, and recorded camera information.
4. Optical center innovations are fused with radiometric innovations.

The monograph supports this separation and the underlying geometry. It does not
specifically validate our implementation, seven-parameter optimizer, covariance,
or EKF numerical results.

| Reference | Relevant content | Implication for DeepNAV |
|---|---|---|
| Sections 3.3-3.5, pp. 38-50 | Digital image formation, detector noise, flattening, readout smear | Truth detector and image calibration should be explicit. |
| Section 4.2, pp. 52-54; section 4.9, pp. 69-70 | TDB exposure midpoint, barycentric observer position/velocity, retarded target, stellar aberration, rotation, projection, pixel conversion | Our geometry order is appropriate; keep physical effects consistent in truth and prediction. |
| Section 5.1, p. 71 | Center extraction assumes known camera, corrected images, and calculated pixel noise | Raw EDR and exact synthetic variance maps must not be confused with realistic calibrated input products. |
| Section 5.3, pp. 75-77 | Lambert, Lommel-Seeliger, Minnaert, and Hapke reflectance laws | Lambert is one controlled brightness assumption, not a universal surface model. |
| Section 5.4, pp. 77-78 | Systematic photocenter displacement must be handled | Do not feed uncorrected brightness centroids into a geometric-center observation model. |
| Section 5.5, pp. 78-81 | Whole-body brightness matching and limb/terminator methods for resolved targets | Our illuminated-body fit is consistent with this class of methods; it minimizes weighted residuals rather than implementing Eq. 5.19's correlation score. |
| Section 5.5, p. 79 | Unmodeled albedo/topography can corrupt whole-image matching | Add truth/model mismatch tests and compare a limb-based method when resolution permits. |
| Sections 6.1, 6.3, pp. 85-91 | Dynamic and optical parameter derivatives | Validate the optical measurement Jacobian and uncertainty propagation before fusion. |
| Section 6.4, pp. 92-93 | Position/pointing correlations, multiple parameter and weighting scenarios | Known attitude is an optimistic baseline. Later model pointing uncertainty and persistent nuisance parameters. |
| Section 6.5, pp. 93-95 | In-flight calibration of focal length, pixel layout, distortion, pointing | Nominal FK/IK and a perfect synthetic camera do not establish calibrated flight accuracy. |

## Recommended truth-versus-estimator differences

Use an initial matched-model baseline, then introduce documented mismatches one
at a time before a combined case. These are engineering recommendations based
on the source and our code review, not a prescribed recipe from the monograph.

- **Detector statistics:** draw Poisson source/background electrons, Gaussian read
  noise, and digitization in truth. In the estimator, compute variance from noisy
  data and nominal gain/read-noise calibration, not from the noiseless truth image.
  Use appropriate count-dependent weighting; low-count tests may need a
  Poisson/read-noise likelihood. Add separate bad-pixel/outlier tests with masks.
- **PSF:** use a mildly elliptical or otherwise perturbed PSF in truth, with the
  nominal calibrated PSF in the estimator. Test width and shape errors separately.
- **Shape:** generate a triaxial ellipsoid using declared/PCK radii, fit a sphere
  initially, then compare an ellipsoid fitter. Terrain can follow at higher resolution.
- **Surface brightness:** introduce a body-fixed albedo pattern or a documented
  alternate reflectance law in truth; retain a simpler fitted brightness law. Avoid
  declaring arbitrary coefficients to be flight-calibrated Enceladus properties.
- **Camera calibration:** draw focal-scale, principal-point, distortion, and mounting
  errors from declared priors. Draw shared calibration errors once per Monte Carlo
  run, rather than as independent per-pixel or per-image noise.
- **Attitude:** retain a known-attitude baseline, then perturb true camera rotations
  relative to the supplied nominal rotations. Preserve temporal correlation where
  intended. Use independent attitude information or estimate/consider pointing
  nuisance parameters so the filter does not interpret every pointing error as
  spacecraft position error.
- **Timing/smear:** integrate truth over the shutter and any declared jitter; use
  exposure-aware prediction or a validated midpoint approximation. Timestamp offsets
  should have an explicit uncertainty model, separate from photon noise.
- **Brightness across dates:** replace identical full-phase counts at each exposure
  with a reference flux scaled by exposure duration and inverse-square range,
  together with the declared phase/reflectance law. This sets a physically more
  meaningful signal-to-noise evolution without changing observation geometry.

Keep reception light time and stellar aberration enabled in both truth and
estimation. Section 4.2 explicitly includes both and uses the Newtonian apparent
vector expression already implemented in DeepNAV. At our tested geometry, its
projection differs from direct CSPICE CN+S by about 0.000035 pixel. More accurate
truth geometry can be useful, but deliberately removing a known correction from
the estimator is an omission/sensitivity experiment, not a nominal realism model.
Optical aberrations (PSF/distortion) and stellar aberration are distinct effects.

## Covariance and navigation validation

The optimizer's local center covariance is conditional on its fixed image model,
PSF, calibration, and supplied pixel variances. It is not the entire optical
observation uncertainty. With independent nuisance errors, a first approximation
is R_effective = R_fit + J_nuisance P_nuisance J_nuisance^T + R_model.
This approximation must not double count nuisance effects already estimated in
the EKF or included in the fit covariance. Shared calibration/pointing errors
correlate measurements across images; represent them with augmented/consider
states or a correlated treatment rather than white-noise inflation alone.
Persistent center bias needs correction or modeling, not just larger weights.

Before EKF fusion:

1. Repeat each image over many independent detector-noise realizations; measure
   center bias, empirical covariance, confidence-ellipse coverage, and failure rate.
2. Repeat over nuisance draws and geometry changes: phase, apparent size, target
   position in the field, signal strength, and pointing conditions.
3. Compare empirical errors with the exported formal/augmented covariance. Inspect
   residual structure, parameter conditioning, and reduced chi-squared; none is an
   accuracy certificate by itself.
4. Run radiometric-only and radiometric-plus-optical EKFs on the same simulated
   truth, checking trajectory error and innovation consistency. Freeze the camera
   rotation while differentiating an observation with respect to spacecraft state.
5. Avoid treating a center fit strongly constrained by the EKF state as independent
   new information without handling that dependence.

The current 0.019 pixel geometric-center error is from one noisy synthetic frame
whose source and fitter are both uniform Lambert spheres with a known PSF. It is
an encouraging model-consistency check, not a demonstrated accuracy or bias floor.

## Apparent issue in photocenter Eq. 5.18

The displayed Eq. 5.18 on printed p. 78 appears to omit a factor of pi. This is
an independent mathematical check, not an official erratum. At phase pi/2, the
printed expression evaluates to 3r/16. For a unit Lambert sphere at quarter phase,
choose projected disk coordinates (x,y), illuminated half x > 0. Radiance is
proportional to x, so the integrated signal is integral(x dA) = 2/3 and its first
moment is integral(x^2 dA) = pi/8. Their ratio is 3pi/16 = 0.5890486 radii,
not 0.1875 radii. Direct numerical quadrature gives 0.5890484.

The current numerical renderer/photocenter calculation does not use Eq. 5.18.
Validate any future analytic photocenter correction against the moment definition
(Eqs. 5.16-5.17) and independent quadrature rather than copying that expression.
