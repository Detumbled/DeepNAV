"""Run the real Cassini image/SPICE comparison and save an overlay and residual table."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import subprocess

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
ARCSEC = 180 * 3600 / np.pi


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def ray(pixel, geometry):
    """Invert the nominal IK pinhole model, in NAC coordinates (sample/line are -X/-Y)."""
    center = np.asarray(geometry["principal_point"])
    xy = -(np.asarray(pixel) - center) * geometry["pixel_pitch_mm"] / geometry["focal_length_mm"]
    vector = np.array([*xy, 1.0])
    return vector / np.linalg.norm(vector)


def compare(geometry, centroid, rows):
    if centroid["mid_time_utc"] != geometry["mid_time_utc"]:
        raise ValueError("Geometry and centroid exposure times differ")
    if centroid["image"] != Path(geometry["image"]).name:
        raise ValueError("Geometry and centroid image identifiers differ")
    midpoint = {r["correction"]: np.array([float(r["sample"]), float(r["line"])])
                for r in rows if float(r["offset_seconds"]) == 0}
    reference = midpoint["CN+S"]
    measured = {
        "Circular Gaussian": np.array([centroid["sample"], centroid["line"]]),
        "Elliptical Gaussian": np.array([centroid["elliptical"]["sample"], centroid["elliptical"]["line"]]),
        "Brightness moment": np.array(centroid["brightness_centroid"]),
    }
    records = []
    reference_ray = ray(reference, geometry)
    for name, pixel in measured.items():
        residual = pixel - reference
        observed_ray = ray(pixel, geometry)
        separation = np.arctan2(np.linalg.norm(np.cross(reference_ray, observed_ray)),
                                np.dot(reference_ray, observed_ray)) * ARCSEC
        records.append(dict(method=name, sample=float(pixel[0]), line=float(pixel[1]),
                            residual_sample_pixels=float(residual[0]), residual_line_pixels=float(residual[1]),
                            residual_norm_pixels=float(np.linalg.norm(residual)), angular_separation_arcsec=float(separation)))
    sweep = np.array([[float(r["sample"]), float(r["line"])] for r in rows if r["correction"] == "CN+S"])
    diagnostics = {
        "reference_pixel_cn_s": reference.tolist(),
        "deepnav_pixel": midpoint["DeepNAV"].tolist(),
        "exposure_endpoint_displacement_pixels": float(np.linalg.norm(sweep[-1] - sweep[0])),
        "centroid_method_spread_pixels": float(max(np.linalg.norm(a-b) for a in measured.values() for b in measured.values())),
        "residuals": records,
        "limitations": [
            "SPK/CK geometry is a reference, not independently exact flight truth.",
            "Nominal IK pinhole only; no calibrated distortion or independently measured pointing covariance.",
            "Centroids are brightness centers from uncalibrated EDR; no photocenter correction was applied.",
            "The brightest-source association with Enceladus has not been independently verified.",
            "Gaussian residuals are large; formal fit covariance is not the total measurement uncertainty.",
            "No pointing adjustment, orbit update, or navigation accuracy claim is made.",
        ],
    }
    return midpoint, measured, diagnostics


def plot(image, geometry, midpoint, measured, diagnostics, destination):
    figure, axes = plt.subplots(1, 2, figsize=(13, 6), constrained_layout=True)
    low, high = np.percentile(image, [1, 99.9])
    for ax in axes:
        ax.imshow(image, cmap="gray", vmin=low, vmax=max(low+1, high), origin="upper", interpolation="nearest")
        ax.scatter(*midpoint["CN+S"], marker="+", s=160, color="lime", label="SPICE CN+S geometric center")
        ax.scatter(*midpoint["DeepNAV"], marker="o", s=120, facecolors="none", edgecolors="cyan", label="DeepNAV geometric center")
        ax.scatter(*midpoint["NONE"], marker="x", s=55, color="yellow", label="SPICE simultaneous geometry (NONE)")
        for (name, pixel), marker, color in zip(measured.items(), ["x", "+", "o"], ["red", "orange", "magenta"]):
            ax.scatter(*pixel, marker=marker, s=70, color=color, label=name)
        ax.set(xlabel="Sample [px; zero-based]", ylabel="Line [px; zero-based]")
        ax.ticklabel_format(useOffset=False, style="plain")
    axes[0].set_title("Original image — display contrast stretched")
    points = np.array([*midpoint.values(), *measured.values()])
    x0, x1 = max(0, int(points[:, 0].min()-15)), min(image.shape[1], int(points[:, 0].max()+16))
    y0, y1 = max(0, int(points[:, 1].min()-15)), min(image.shape[0], int(points[:, 1].max()+16))
    local = image[y0:y1, x0:x1]
    local_low, local_high = np.percentile(local, [10, 99.9])
    axes[1].images[0].set_clim(local_low, max(local_low+1, local_high))
    axes[1].set_xlim(points[:, 0].min()-15, points[:, 0].max()+15)
    axes[1].set_ylim(points[:, 1].max()+15, points[:, 1].min()-15)
    axes[1].set_title("Independent CK pointing; local display contrast")
    axes[1].legend(fontsize=8, loc="upper left")
    figure.suptitle(f"Cassini Enceladus image versus SPK/CK reference\n{Path(geometry['image']).name} · {geometry['mid_time_utc']} UTC")
    first = diagnostics["residuals"][0]
    figure.supxlabel(f"Circular centroid − reference: ({first['residual_sample_pixels']:.3f}, {first['residual_line_pixels']:.3f}) px; "
                      f"separation {first['angular_separation_arcsec']:.2f} arcsec\n"
                      "Nominal IK pinhole; uncalibrated brightness centers; photocenter and pointing errors unresolved.", fontsize=10)
    figure.savefig(destination, dpi=160)
    plt.close(figure)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build-dir", type=Path, default=ROOT / "build-clang")
    parser.add_argument("--image", type=Path, default=ROOT / "Input Optical/N1476124698_2.IMG")
    parser.add_argument("--kernels", type=Path, default=ROOT / "Kernels/Cassini/cassini_2004_opnav.tm")
    parser.add_argument("--output", type=Path, default=ROOT / "results/cassini_optical_validation")
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    # Recompute both paths for this exposure; do not reuse potentially stale centroid outputs.
    subprocess.run([str(args.build_dir.resolve() / "cassini_optical_validation"), "--image", str(args.image.resolve()),
                    "--kernels", str(args.kernels.resolve()), "--output", str(output)], check=True)
    centroid_dir = output / "centroid"
    subprocess.run([str(args.build_dir.resolve() / "cassini_centroid_demo"), str(args.image.resolve()),
                    "--output", str(centroid_dir)], check=True)
    geometry = json.loads((output / "geometry.json").read_text())
    centroid = json.loads((centroid_dir / "summary.json").read_text())
    with (output / "predictions.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    midpoint, measured, diagnostics = compare(geometry, centroid, rows)
    # Pin the actual kernel bytes and image used, not only their filenames.
    files = [args.image.resolve(), args.image.resolve().with_suffix(".LBL"),
             *map(Path, geometry["loaded_kernels"])]
    diagnostics["sha256"] = {str(path): digest(path) for path in files}
    (output / "comparison.json").write_text(json.dumps(diagnostics, indent=2) + "\n")
    with (output / "residuals.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(diagnostics["residuals"][0]))
        writer.writeheader()
        writer.writerows(diagnostics["residuals"])
    plot(np.loadtxt(centroid_dir / "image.csv", delimiter=","), geometry, midpoint, measured,
         diagnostics, output / "image_spice_overlay.png")
    for record in diagnostics["residuals"]:
        print(f"{record['method']}: Δsample={record['residual_sample_pixels']:.4f} px, "
              f"Δline={record['residual_line_pixels']:.4f} px, separation={record['angular_separation_arcsec']:.3f} arcsec")
    print(f"Exposure endpoint displacement: {diagnostics['exposure_endpoint_displacement_pixels']:.6f} px")
    print(f"Centroid method spread: {diagnostics['centroid_method_spread_pixels']:.6f} px")
    print(f"Saved comparison and overlay in {output}")


if __name__ == "__main__":
    main()
