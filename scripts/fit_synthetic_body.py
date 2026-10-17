"""Fit a Gaussian and a Lambert sphere to one synthetic image; truth is evaluation-only."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import least_squares
from scipy.special import erf

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path, nargs="?", default=ROOT / "results/cassini_synthetic_optical/simulated_sequence")
    parser.add_argument("--frame", type=int, default=0)
    parser.add_argument("--quadrature-order", type=int, default=24)
    parser.add_argument("--psf-sigma", type=float, default=.8, help="Independently calibrated Gaussian PSF width in pixels")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if not 12 <= args.quadrature_order <= 64:
        parser.error("--quadrature-order must be between 12 and 64")
    if not np.isfinite(args.psf_sigma) or args.psf_sigma <= 0:
        parser.error("--psf-sigma must be positive and finite")
    directory = args.directory.resolve()
    # The pixel file is the fitting input. Truth centers, phase, noiseless source,
    # and simulated spacecraft states are not used to initialize or fit the model.
    pixels = np.genfromtxt(directory / f"frame_{args.frame}.csv", delimiter=",", names=True)
    width, height = int(pixels["sample"].max())+1, int(pixels["line"].max())+1
    image = pixels["dn"].reshape(height, width)
    variance = pixels["variance_dn2"].reshape(height, width)
    mask = pixels["mask"].reshape(height, width).astype(bool)
    search = np.where(mask, -np.inf, image)
    peak_line, peak_sample = np.unravel_index(np.argmax(search), image.shape)
    x0, x1 = max(0, peak_sample-18), min(width, peak_sample+19)
    y0, y1 = max(0, peak_line-18), min(height, peak_line+19)
    observed = image[y0:y1, x0:x1]
    weights = 1/np.sqrt(np.maximum(variance[y0:y1, x0:x1], 1e-12))
    valid = ~mask[y0:y1, x0:x1]
    yy, xx = np.mgrid[y0:y1, x0:x1]
    border = np.r_[observed[0], observed[-1], observed[:, 0], observed[:, -1]]
    background = float(np.median(border))
    signal = np.maximum(observed-background-3*np.median(1/weights), 0)
    centroid = np.array([(xx*signal).sum(), (yy*signal).sum()])/signal.sum()
    centered = np.stack([xx-centroid[0], yy-centroid[1]], axis=-1)
    moments = np.einsum("abi,abj,ab->ij", centered, centered, signal)/signal.sum()
    _, eigenvectors = np.linalg.eigh(moments)
    minor = eigenvectors[:, 0]
    azimuth = np.arctan2(minor[1], minor[0])
    amplitude = float(observed.max()-background)

    def gaussian(parameters):
        x, y, height_dn, sigma, bg = parameters
        ex = erf((xx+.5-x)/(np.sqrt(2)*sigma))-erf((xx-.5-x)/(np.sqrt(2)*sigma))
        ey = erf((yy+.5-y)/(np.sqrt(2)*sigma))-erf((yy-.5-y)/(np.sqrt(2)*sigma))
        return bg+height_dn*np.pi*sigma*sigma/2*ex*ey

    def residual(model, parameters):
        return ((model(parameters)-observed)*weights)[valid]

    gaussian_fit = least_squares(lambda p: residual(gaussian, p), [*centroid, amplitude, 3, background],
        bounds=([x0, y0, 0, .2, 0], [x1, y1, 1e5, 15, 1e4]), x_scale="jac", max_nfev=300)
    # PSF width is assumed calibrated. Its value is a declared camera input,
    # not recovered from the source truth. Lighting angles are fitted freely.
    psf_sigma = args.psf_sigma
    # Integrate over a fixed unit-disk quadrature, then integrate the PSF
    # analytically over detector pixels. The fixed nodes keep center/size
    # derivatives smooth, unlike hard subpixel masks at a moving limb.
    order = args.quadrature_order
    radial, weights_r = np.polynomial.legendre.leggauss(order)
    radial = (radial+1)/2
    weights_r = weights_r/2*radial
    angles = (np.arange(4*order)+.5)*2*np.pi/(4*order)
    u = (radial[:, None]*np.cos(angles)).ravel()
    v = (radial[:, None]*np.sin(angles)).ravel()
    disk_weights = np.broadcast_to(weights_r[:, None]*2*np.pi/(4*order), (order, 4*order)).ravel()
    normal_z = np.sqrt(np.maximum(0, 1-u*u-v*v))
    samples, lines = np.arange(x0, x1)[:, None], np.arange(y0, y1)[:, None]

    def sphere(parameters):
        x, y, radius, phase, sun_angle, height_dn, bg = parameters
        brightness = np.maximum(0, np.sin(phase)*(u*np.cos(sun_angle)+v*np.sin(sun_angle))
                                + np.cos(phase)*normal_z)
        source_x, source_y = x+radius*u, y+radius*v
        ex = .5*(erf((samples+.5-source_x)/(np.sqrt(2)*psf_sigma))-erf((samples-.5-source_x)/(np.sqrt(2)*psf_sigma)))
        ey = .5*(erf((lines+.5-source_y)/(np.sqrt(2)*psf_sigma))-erf((lines-.5-source_y)/(np.sqrt(2)*psf_sigma)))
        return bg+(ey*(height_dn*radius*radius*disk_weights*brightness))@ex.T

    candidates = []
    radius_guess = max(3., min(12., .5*max(np.ptp(xx[signal > 0]), np.ptp(yy[signal > 0]))))
    for angle in (azimuth, azimuth+np.pi):
        center_guess = centroid-.35*radius_guess*np.array([np.cos(angle), np.sin(angle)])
        initial = [*center_guess, radius_guess, np.radians(65), angle, amplitude*1.2, background]
        fitted = least_squares(lambda p: residual(sphere, p), initial,
            bounds=([x0, y0, 2, .05, -4*np.pi, 0, 0], [x1, y1, 15, np.pi-.05, 4*np.pi, 1e5, 1e4]),
            x_scale=[1, 1, 1, .5, .5, amplitude, 10], diff_step=1e-5, max_nfev=250,
            ftol=1e-8, xtol=1e-8, gtol=1e-6)
        candidates.append(fitted)
    physical_fit = min(candidates, key=lambda fit: np.dot(fit.fun, fit.fun))
    if not gaussian_fit.success or not physical_fit.success:
        raise RuntimeError(f"Fit failed: Gaussian={gaussian_fit.message}; sphere={physical_fit.message}")
    gaussian_model, sphere_model = gaussian(gaussian_fit.x), sphere(physical_fit.x)
    _, singular_values, right_vectors = np.linalg.svd(physical_fit.jac, full_matrices=False)
    if singular_values[-1] <= singular_values[0]*1e-12:
        raise RuntimeError("Illuminated sphere fit is locally unobservable")
    parameter_covariance = (right_vectors.T/singular_values**2)@right_vectors
    # Only now read ground truth for evaluation and plotting.
    summary = json.loads((directory / "summary.json").read_text())
    truth = summary["frames"][args.frame]
    geometric = np.asarray(truth["geometric_center"])
    brightness_center = np.asarray(truth["brightness_center"])
    results = {
        "frame": args.frame, "utc": truth["utc"], "quadrature_order": order, "known_psf_sigma_pixels": psf_sigma,
        "truth_geometric_center": geometric.tolist(), "truth_brightness_center": brightness_center.tolist(),
        "gaussian_center": gaussian_fit.x[:2].tolist(), "sphere_geometric_center": physical_fit.x[:2].tolist(),
        "gaussian_geometric_error_pixels": float(np.linalg.norm(gaussian_fit.x[:2]-geometric)),
        "sphere_geometric_error_pixels": float(np.linalg.norm(physical_fit.x[:2]-geometric)),
        "gaussian_reduced_chi_squared": float(np.dot(gaussian_fit.fun, gaussian_fit.fun)/(valid.sum()-5)),
        "sphere_reduced_chi_squared": float(np.dot(physical_fit.fun, physical_fit.fun)/(valid.sum()-7)),
        "fitted_radius_pixels": float(physical_fit.x[2]), "fitted_phase_degrees": float(np.degrees(physical_fit.x[3])),
        "truth_phase_degrees": truth["phase_degrees"], "fitted_sun_image_angle_degrees": float(np.degrees(physical_fit.x[4]) % 360),
        "sphere_function_evaluations": physical_fit.nfev,
        "formal_center_covariance_pixels2": parameter_covariance[:2, :2].tolist(),
        "assumptions": "Uniform Lambert sphere; independently calibrated Gaussian PSF; fixed supplied detector variances; distant-body projection. Truth used only after fitting. Formal covariance excludes model/calibration error.",
    }
    output = args.output.resolve() if args.output else directory / f"body_fit_{args.frame}"
    output.mkdir(exist_ok=True, parents=True)
    (output / "comparison.json").write_text(json.dumps(results, indent=2)+"\n")
    np.savez_compressed(output / "fit_arrays.npz", observed=observed, gaussian=gaussian_model, sphere=sphere_model,
                        variance=variance[y0:y1, x0:x1], origin=np.array([x0, y0]))
    figure, axes = plt.subplots(2, 3, figsize=(14, 9), constrained_layout=True)
    extent = (x0-.5, x1-.5, y1-.5, y0-.5)
    for ax, data, title in zip(axes[0], (observed, gaussian_model, sphere_model),
            ("Observed synthetic image", "Fitted pixel-integrated Gaussian", "Fitted illuminated sphere + PSF")):
        raster = ax.imshow(data, origin="upper", extent=extent, cmap="inferno", vmin=background, vmax=observed.max(), interpolation="nearest")
        figure.colorbar(raster, ax=ax, label="DN", shrink=.8)
        ax.set_title(title)
    axes[1, 0].imshow(observed, origin="upper", extent=extent, cmap="gray", interpolation="nearest")
    for point, marker, color, name in zip((geometric, brightness_center, gaussian_fit.x[:2], physical_fit.x[:2]),
            ("+", "o", "x", "+"), ("cyan", "lime", "orange", "magenta"),
            ("True geometric center", "True brightness center", "Gaussian centroid", "Sphere geometric center")):
        axes[1, 0].scatter(*point, marker=marker, s=130, color=color, label=name,
                            facecolors="none" if marker == "o" else color)
    axes[1, 0].legend(fontsize=8)
    axes[1, 0].set_title("Center comparison")
    bound = max(float(np.max(np.abs(observed-gaussian_model))), float(np.max(np.abs(observed-sphere_model))))
    for ax, model, name, chi in zip(axes[1, 1:], (gaussian_model, sphere_model), ("Gaussian", "Sphere"),
            (results["gaussian_reduced_chi_squared"], results["sphere_reduced_chi_squared"])):
        raster = ax.imshow(observed-model, origin="upper", extent=extent, cmap="RdBu_r", vmin=-bound, vmax=bound, interpolation="nearest")
        figure.colorbar(raster, ax=ax, label="Observed − model [DN]", shrink=.8)
        ax.set_title(f"{name} residual · reduced χ²={chi:.2f}")
    for ax in axes.flat:
        ax.set(xlim=(geometric[0]-13, geometric[0]+13), ylim=(geometric[1]+13, geometric[1]-13), xlabel="Sample [px]", ylabel="Line [px]")
    figure.suptitle(f"Geometric-center recovery from a physical brightness model\n{truth['utc'][:19]} UTC · synthetic frame {args.frame}")
    figure.supxlabel(f"Geometric-center error: Gaussian {results['gaussian_geometric_error_pixels']:.4f} px; "
                      f"illuminated sphere {results['sphere_geometric_error_pixels']:.4f} px\n"
                      f"Center, radius, phase, Sun position angle, amplitude, and background fitted from pixels; PSF fixed at {psf_sigma:g} px.", fontsize=10)
    figure.savefig(output / "body_fit_comparison.png", dpi=160)
    plt.close(figure)
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
