"""Compare circular/elliptical fits, brightness centroid, and window sensitivity."""
import argparse
import json
from pathlib import Path

import matplotlib
import numpy as np

DEFAULT_DIRECTORY = Path(__file__).resolve().parent / "Output optical" / "cassini"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path, nargs="?", default=DEFAULT_DIRECTORY)
    parser.add_argument("--no-show", action="store_true", help="Save figures without opening windows")
    args = parser.parse_args()
    for name in ("image.csv", "fit.csv", "summary.json", "window_comparison.csv"):
        if not (args.directory / name).is_file():
            parser.error(f"Missing {name}; run ./build-clang/cassini_centroid_demo first.")
    if args.no_show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    image = np.loadtxt(args.directory / "image.csv", delimiter=",")
    pixels = np.genfromtxt(args.directory / "fit.csv", delimiter=",", names=True)
    metadata = json.loads((args.directory / "summary.json").read_text())
    if "elliptical" not in metadata or "elliptical_model_dn" not in pixels.dtype.names:
        parser.error("Old demo outputs; rebuild and run cassini_centroid_demo again.")
    xs, ys = np.unique(pixels["sample"]), np.unique(pixels["line"])
    shape = (len(ys), len(xs))
    extent = (xs[0]-.5, xs[-1]+.5, ys[-1]+.5, ys[0]-.5)
    observed, circular, elliptical = [pixels[key].reshape(shape) for key in ("dn", "model_dn", "elliptical_model_dn")]
    mask = pixels["masked"].reshape(shape).astype(bool)
    residuals = [np.ma.masked_array(pixels[key].reshape(shape), mask)
                 for key in ("residual_dn", "elliptical_residual_dn")]
    ellipse = metadata["elliptical"]
    centers = [(metadata["sample"], metadata["line"]),
               (ellipse["sample"], ellipse["line"]), metadata["brightness_centroid"]]
    names = ["Circular", "Elliptical", "Brightness moment"]
    colors, markers = ["cyan", "lime", "white"], ["x", "+", "o"]
    fig, axes = plt.subplots(2, 3, figsize=(15, 9), constrained_layout=True)
    fig.suptitle(f"Cassini ISS NAC — circular versus rotated elliptical Gaussian\n{metadata['image']}; {metadata['mid_time_utc']} UTC", fontsize=13)
    lo, hi = min(v.min() for v in (observed, circular, elliptical)), max(v.max() for v in (observed, circular, elliptical))
    titles = [f"Observed {len(xs)} × {len(ys)} window",
              f"Circular model: σ = {metadata['sigma_pixels']:.3f} px",
              f"Elliptical model: σ = ({ellipse['sigma_major_pixels']:.3f}, {ellipse['sigma_minor_pixels']:.3f}) px\n"
              f"Angle = {np.degrees(ellipse['angle_radians']):.1f}° toward +line"]
    for column, (values, title) in enumerate(zip((observed, circular, elliptical), titles)):
        ax = axes[0, column]
        raster = ax.imshow(values, origin="upper", extent=extent, cmap="inferno", vmin=lo, vmax=hi, interpolation="nearest")
        selected = range(3) if column == 0 else [column-1]
        for k in selected:
            ax.scatter(*centers[k], marker=markers[k], color=colors[k], s=90,
                       facecolors="none" if k == 2 else colors[k], label=names[k])
        ax.legend(fontsize=8, loc="upper right")
        ax.set_title(title, fontsize=10)
        fig.colorbar(raster, ax=ax, label="Raw DN", shrink=.8)
    # Stretch only the full-image display, leaving all fit data unchanged.
    low, high = np.percentile(image, [1, 99.9])
    raster = axes[1, 0].imshow(image, origin="upper", cmap="gray", vmin=low, vmax=max(low+1, high), interpolation="nearest")
    axes[1, 0].add_patch(Rectangle((xs[0]-.5, ys[0]-.5), len(xs), len(ys), fill=False, edgecolor="lime"))
    axes[1, 0].set_title("Original image — display contrast stretched", fontsize=10)
    fig.colorbar(raster, ax=axes[1, 0], label="Raw DN", shrink=.8)
    bound = max(1, *(float(np.ma.max(np.ma.abs(r))) for r in residuals))
    for column, (residual, fit) in enumerate(zip(residuals, (metadata, ellipse)), start=1):
        ax = axes[1, column]
        raster = ax.imshow(residual, origin="upper", extent=extent, cmap="RdBu_r", vmin=-bound, vmax=bound, interpolation="nearest")
        rms = float(np.sqrt(np.ma.mean(residual**2)))
        ax.set_title(f"{names[column-1]} residual: RMS = {rms:.1f} DN\nReduced χ² = {fit['reduced_chi_squared']:.2f}",
                     color="darkred" if fit["reduced_chi_squared"] > 3 else "black", fontsize=10)
        fig.colorbar(raster, ax=ax, label="Observed − model [DN]", shrink=.8)
    for ax in axes.flat:
        ax.set_xlabel("Sample [px; zero-based]")
        ax.set_ylabel("Line [px; zero-based]")
        ax.ticklabel_format(useOffset=False, style="plain")
    shift = np.linalg.norm(np.array(centers[1])-centers[0])
    fig.supxlabel(f"Circular = ({centers[0][0]:.4f}, {centers[0][1]:.4f}); elliptical = ({centers[1][0]:.4f}, {centers[1][1]:.4f}) px; separation = {shift:.3f} px\n"
                  "Uncalibrated EDR; estimated noise. Brightness centers are not verified geometric centers.", fontsize=10)
    destination = args.directory / "cassini_centroid.png"
    fig.savefig(destination, dpi=160)
    print(f"Saved {destination}")
    print(f"Circular centroid: {centers[0]}; elliptical centroid: {centers[1]}; brightness centroid: {centers[2]}")

    windows = np.atleast_1d(np.genfromtxt(args.directory / "window_comparison.csv", delimiter=",", names=True))
    stability, charts = plt.subplots(1, 2, figsize=(11, 4.5), constrained_layout=True)
    stability.suptitle("Centroid sensitivity to fitting-window size")
    for ax, coordinate in zip(charts, ("sample", "line")):
        for prefix, label, color in zip(("circular", "elliptical", "moment"), names, ("tab:blue", "tab:green", "tab:orange")):
            ax.plot(2*windows["radius"]+1, windows[f"{prefix}_{coordinate}"], "o-", label=label, color=color)
        ax.set_xlabel("Requested window side [px; clipped at image edges]")
        ax.set_ylabel(f"{coordinate.capitalize()} coordinate [px]")
        ax.ticklabel_format(useOffset=False, style="plain")
        ax.grid(alpha=.25)
        ax.legend(fontsize=9)
    stability.supxlabel("Same seed, background estimate, masks, and per-pixel noise formula across windows. Stability does not establish accuracy.", fontsize=9)
    destination = args.directory / "window_sensitivity.png"
    stability.savefig(destination, dpi=160)
    print(f"Saved {destination}")
    if not args.no_show:
        plt.show()
    plt.close(fig)
    plt.close(stability)


if __name__ == "__main__":
    main()
