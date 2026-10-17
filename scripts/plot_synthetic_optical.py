"""Plot SPK-driven synthetic images and direct pixel observations; no GUI required."""
import argparse
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path, nargs="?", default=ROOT / "results/cassini_synthetic_optical")
    args = parser.parse_args()
    directory = args.directory.resolve()
    required = [directory / name for name in ("frame.csv", "measurements.csv", "summary.json")]
    for path in required:
        if not path.is_file():
            parser.error(f"Missing {path}; run ./build-clang/cassini_synthetic_optical first")
    frame = np.genfromtxt(required[0], delimiter=",", names=True)
    measurements = np.genfromtxt(required[1], delimiter=",", names=True)
    summary = json.loads(required[2].read_text())
    width, height = int(frame["sample"].max())+1, int(frame["line"].max())+1
    def raster(key):
        return frame[key].reshape(height, width)
    source, observed = raster("source_electrons"), raster("dn")
    geometric = np.asarray(summary["geometric_center"])
    brightness = np.asarray(summary["brightness_center"])
    fitted = np.asarray(summary["fitted_center"])
    figure, axes = plt.subplots(2, 2, figsize=(12, 9), constrained_layout=True)
    points = np.array([geometric, brightness, fitted])
    for ax, data, title, unit in zip(axes[0], (source, observed),
            ("Expected source: illuminated sphere + PSF", "Observed image: photon/read noise + digitization"), ("Electrons", "DN")):
        image = ax.imshow(data, origin="upper", cmap="inferno", interpolation="nearest")
        figure.colorbar(image, ax=ax, label=unit, shrink=.85)
        ax.set_xlim(points[:, 0].min()-13, points[:, 0].max()+13)
        ax.set_ylim(points[:, 1].max()+13, points[:, 1].min()-13)
        ax.set(title=title, xlabel="Sample [px]", ylabel="Line [px]")
    axes[1, 0].imshow(source, origin="upper", cmap="gray", interpolation="nearest")
    for point, marker, color, label in zip(points, ("+", "o", "x"), ("cyan", "lime", "orange"),
            ("Geometric center (SPK truth)", "Noiseless brightness center", "Circular Gaussian fit")):
        axes[1, 0].scatter(*point, marker=marker, s=130, color=color, label=label,
                           facecolors="none" if marker == "o" else color)
    axes[1, 0].set_xlim(points[:, 0].min()-6, points[:, 0].max()+6)
    axes[1, 0].set_ylim(points[:, 1].max()+6, points[:, 1].min()-6)
    axes[1, 0].set(title="Geometric center and brightness center are distinct", xlabel="Sample [px]", ylabel="Line [px]")
    axes[1, 0].legend(fontsize=8)
    valid = (frame["masked"] == 0) & (frame["variance_dn2"] > 0)
    normalized = (frame["dn"][valid]-frame["expected_dn"][valid])/np.sqrt(frame["variance_dn2"][valid])
    residual_map = np.full(len(frame), np.nan)
    residual_map[valid] = normalized
    image = axes[1, 1].imshow(residual_map.reshape(height, width), origin="upper", cmap="RdBu_r", vmin=-3, vmax=3)
    figure.colorbar(image, ax=axes[1, 1], label="(Observed − expected) / predicted σ", shrink=.85)
    axes[1, 1].set(title=f"Digitized detector residuals: mean={np.mean(normalized):.3f}, std={np.std(normalized):.3f}",
                   xlabel="Sample [px]", ylabel="Line [px]")
    figure.suptitle(f"SPK-driven synthetic {summary['body']} image · {summary['mid_time_utc']} UTC\n"
                   f"Fixed synthetic pointing · Lambert sphere · phase {summary['phase_angle_degrees']:.2f}° · seed {summary['seed']}")
    photocenter_offset = float(np.linalg.norm(brightness-geometric))
    fit_offset = float(np.linalg.norm(fitted-brightness))
    figure.supxlabel(f"Brightness − geometric separation: {photocenter_offset:.3f} px; Gaussian fit − brightness separation: {fit_offset:.3f} px\n"
                     f"PSF σ={summary['psf_sigma_pixels']:.2f} px; exposure={summary['exposure_seconds']:.3f} s; "
                     "declared photon budget, not flight radiometric calibration.", fontsize=10)
    figure.savefig(directory / "synthetic_frame.png", dpi=160)
    plt.close(figure)

    figure, axes = plt.subplots(2, 2, figsize=(12, 7), constrained_layout=True)
    errors = []
    for column, coordinate in enumerate(("sample", "line")):
        t = measurements["offset_seconds"]
        truth = measurements[f"truth_{coordinate}"]
        observed = measurements[coordinate]
        sigma = np.sqrt(measurements[f"variance_{coordinate}"])
        residual = observed-truth
        errors.append(residual)
        axes[0, column].plot(t, truth, label="SPK-derived geometric center")
        axes[0, column].scatter(t, observed, s=12, label="Synthetic measurement")
        axes[0, column].set(ylabel=f"{coordinate.title()} [px]", xlabel="Seconds from reference exposure")
        axes[0, column].legend(fontsize=8)
        axes[1, column].scatter(t, residual, s=14, label="Measurement − truth")
        axes[1, column].plot(t, 3*sigma, "--", color="gray", label="±3σ")
        axes[1, column].plot(t, -3*sigma, "--", color="gray")
        axes[1, column].axhline(0, color="black", lw=.7)
        axes[1, column].set(ylabel=f"{coordinate.title()} residual [px]", xlabel="Seconds from reference exposure")
        axes[1, column].legend(fontsize=8)
    figure.suptitle("Direct synthetic pixel observations · same fixed camera attitude over the entire sequence\n"
                   "Geometric centers + independent 0.15 px Gaussian noise; no image or photocenter bias in this stream")
    figure.savefig(directory / "synthetic_measurements.png", dpi=160)
    plt.close(figure)
    diagnostics = {
        "photocenter_offset_pixels": photocenter_offset,
        "gaussian_fit_minus_brightness_pixels": fit_offset,
        "direct_pixel_noise_rms_sample": float(np.sqrt(np.mean(errors[0]**2))),
        "direct_pixel_noise_rms_line": float(np.sqrt(np.mean(errors[1]**2))),
        "detector_normalized_mean": float(np.mean(normalized)),
        "detector_normalized_std": float(np.std(normalized)),
    }
    diagnostics["kernel_sha256"] = {}
    for filename in summary["loaded_kernels"]:
        with Path(filename).open("rb") as stream:
            diagnostics["kernel_sha256"][filename] = hashlib.file_digest(stream, "sha256").hexdigest()
    (directory / "plot_diagnostics.json").write_text(json.dumps(diagnostics, indent=2)+"\n")
    print(json.dumps({key: value for key, value in diagnostics.items() if key != "kernel_sha256"}, indent=2))
    print(f"Saved synthetic_frame.png and synthetic_measurements.png in {directory}")


if __name__ == "__main__":
    main()
