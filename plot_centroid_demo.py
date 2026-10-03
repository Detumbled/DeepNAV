"""Save and display the centroid demo in Matplotlib; requires NumPy/Matplotlib."""
import argparse
import csv
from pathlib import Path

import matplotlib
from matplotlib.patches import Ellipse
import numpy as np

DEFAULT_DIRECTORY = Path(__file__).resolve().parent / "Output optical"

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv", type=Path, nargs="?", default=DEFAULT_DIRECTORY / "centroid_demo.csv")
    parser.add_argument("output", type=Path, nargs="?", help="Saved figure path; defaults beside the CSV")
    parser.add_argument("--no-show", action="store_true", help="Save without opening a window")
    args = parser.parse_args()
    if not args.csv.is_file():
        parser.error(f"Missing {args.csv}. Run ./build-clang/test_centroid_estimator first.")
    if args.no_show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    output = args.output if args.output else args.csv.parent / "centroid_demo.png"
    with args.csv.open(newline="") as source:
        rows = list(csv.DictReader(source))
    if not rows:
        parser.error("CSV contains no pixels")
    metadata = rows[0]
    width = 1 + max(int(row["sample"]) for row in rows)
    height = 1 + max(int(row["line"]) for row in rows)
    image = np.full((height, width), np.nan)
    for row in rows:
        image[int(row["line"]), int(row["sample"])] = float(row["dn"])
    truth = np.array([float(metadata["true_sample"]), float(metadata["true_line"])])
    measured = np.array([float(metadata["fit_sample"]), float(metadata["fit_line"])])
    covariance = np.array([
        [float(metadata["variance_sample"]), float(metadata["covariance_sample_line"])],
        [float(metadata["covariance_sample_line"]), float(metadata["variance_line"])]])
    values, vectors = np.linalg.eigh(covariance)
    angle = np.degrees(np.arctan2(vectors[1, 1], vectors[0, 1]))
    error = np.linalg.norm(measured - truth)

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.8), constrained_layout=True)
    fig.suptitle("Unresolved moon: geometry → noisy pixels → Gaussian centroid", fontsize=13)
    raster = axes[0].imshow(image, origin="upper", cmap="inferno", interpolation="nearest")
    axes[0].scatter(*truth, marker="+", s=110, linewidth=1.7, color="#66e0ff", label="Injected center")
    axes[0].scatter(*measured, marker="x", s=65, linewidth=1.3, color="white", label="Fitted center")
    sigma = metadata.get("psf_sigma_pixels")
    title = f"{width} × {height} pixel window"
    if sigma is not None:
        title += f"; circular PSF σ = {float(sigma):g} px"
    axes[0].set_title(title, fontsize=10)
    axes[0].set_xlabel("Sample [px]")
    axes[0].set_ylabel("Line [px]")
    axes[0].legend(loc="upper left", fontsize=8, facecolor="#222222", labelcolor="white")
    fig.colorbar(raster, ax=axes[0], label="Calibrated intensity [DN]", shrink=0.8)

    # The ellipse is a joint 95% region for a 2D Gaussian, not a 1D 2-sigma bar.
    ellipse = Ellipse(measured, 2*np.sqrt(5.991*values[1]), 2*np.sqrt(5.991*values[0]),
                      angle=angle, facecolor="#e8f2ff", edgecolor="#2867b2", linewidth=1.5,
                      label="Joint 95% fit region")
    axes[1].add_patch(ellipse)
    axes[1].scatter(*truth, marker="+", s=120, linewidth=2, color="#dd7622", label="Injected center")
    axes[1].scatter(*measured, marker="x", s=65, linewidth=1.6, color="#2867b2", label="Fitted center")
    extent = max(0.015, 3*np.sqrt(values[-1]), 2*error)
    axes[1].set_xlim(truth[0]-extent, truth[0]+extent)
    axes[1].set_ylim(truth[1]+extent, truth[1]-extent)
    axes[1].set_aspect("equal")
    axes[1].set_title(f"Subpixel detail: error = {error:.6f} px", fontsize=10)
    axes[1].set_xlabel("Sample [px]")
    axes[1].set_ylabel("Line [px]")
    axes[1].ticklabel_format(useOffset=False, style="plain")
    axes[1].grid(alpha=0.2)
    axes[1].legend(loc="upper right", fontsize=8)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=180)
    print(f"Saved {output}")
    if not args.no_show:
        plt.show()
    plt.close(fig)


if __name__ == "__main__":
    main()
