"""Plot four synthetic images generated from the retained EKF simulation truth."""
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
    parser.add_argument("directory", type=Path, nargs="?", default=ROOT / "results/cassini_synthetic_optical/simulated_sequence")
    args = parser.parse_args()
    directory = args.directory.resolve()
    if not (directory / "summary.json").is_file():
        parser.error("Run ./build-clang/cassini_synthetic_sequence first")
    summary = json.loads((directory / "summary.json").read_text())
    frames = [np.genfromtxt(directory / frame["file"], delimiter=",", names=True) for frame in summary["frames"]]
    # Shared scaling makes brightness changes visible instead of stretching every frame separately.
    vmin = min(float(pixels["dn"].min()) for pixels in frames)
    vmax = max(float(pixels["dn"].max()) for pixels in frames)
    figure, axes = plt.subplots(2, 2, figsize=(12, 10), constrained_layout=True)
    diagnostics = []
    for ax, frame, pixels in zip(axes.flat, summary["frames"], frames):
        width, height = int(pixels["sample"].max())+1, int(pixels["line"].max())+1
        source = pixels["dn"].reshape(height, width)
        geometric, brightness, fitted = [np.asarray(frame[key]) for key in ("geometric_center", "brightness_center", "fitted_center")]
        raster = ax.imshow(source, origin="upper", cmap="inferno", interpolation="nearest", vmin=vmin, vmax=vmax)
        figure.colorbar(raster, ax=ax, label="DN", shrink=.8)
        for point, marker, color, label in zip((geometric, brightness, fitted), ("+", "o", "x"), ("cyan", "lime", "white"),
                ("Simulated geometric center", "Noiseless brightness center", "Gaussian fit")):
            ax.scatter(*point, marker=marker, s=110, color=color, label=label, facecolors="none" if marker == "o" else color)
        ax.set_xlim(geometric[0]-14, geometric[0]+14)
        ax.set_ylim(geometric[1]+14, geometric[1]-14)
        ax.set(title=f"{frame['utc'][:19]} UTC\nPhase {frame['phase_degrees']:.2f}° · independent seed {frame['seed']}",
               xlabel="Sample [px]", ylabel="Line [px]")
        diagnostics.append({"utc": frame["utc"], "photocenter_offset_pixels": float(np.linalg.norm(brightness-geometric)),
                            "fit_minus_brightness_pixels": float(np.linalg.norm(fitted-brightness)),
                            "source_electrons": float(pixels["source_electrons"].sum()),
                            "peak_dn": float(pixels["dn"].max()),
                            "simulated_vs_spk_spacecraft_km": frame["spk_position_difference_km"]})
    axes[0, 0].legend(fontsize=8, loc="upper right")
    figure.suptitle("Four Enceladus images · simulated Cassini truth shared with the radio EKF\n"
                   "One UTC day apart; target re-centered at each exposure; shared DN color scale")
    figure.savefig(directory / "simulated_sequence.png", dpi=160)
    plt.close(figure)
    reference = frames[0]["source_electrons"]
    differences = [pixels["source_electrons"]-reference for pixels in frames[1:]]
    bound = max(float(np.max(np.abs(difference))) for difference in differences)
    figure, axes = plt.subplots(1, 3, figsize=(13, 5), constrained_layout=True)
    for index, (ax, difference) in enumerate(zip(axes, differences), start=1):
        raster = ax.imshow(difference.reshape(height, width), origin="upper", cmap="RdBu_r",
                           vmin=-bound, vmax=bound, interpolation="nearest")
        ax.set_xlim(geometric[0]-14, geometric[0]+14)
        ax.set_ylim(geometric[1]+14, geometric[1]-14)
        ax.set(title=f"{summary['frames'][index]['utc'][:10]} − first exposure", xlabel="Sample [px]", ylabel="Line [px]")
        diagnostics[index]["source_relative_l2_difference_from_first"] = float(np.linalg.norm(difference)/np.linalg.norm(reference))
    figure.colorbar(raster, ax=axes.tolist(), label="Difference in expected source electrons", shrink=.8)
    figure.suptitle("Changes between independently rendered scenes\nNoiseless source differences: geometry and illumination, excluding random detector noise")
    figure.savefig(directory / "scene_differences.png", dpi=160)
    plt.close(figure)
    observations = np.genfromtxt(directory / "observations.csv", delimiter=",", names=True, dtype=None, encoding="utf-8")
    figure, axes = plt.subplots(1, 2, figsize=(11, 4), constrained_layout=True)
    days = (observations["tdb_seconds"]-observations["tdb_seconds"][0])/86400
    center = np.array(summary["camera"]["principal_point"])
    for column, coordinate in enumerate(("sample", "line")):
        sigma = np.sqrt(observations["cov_ss" if column == 0 else "cov_ll"])
        axes[column].scatter(days, observations[coordinate]-center[column], label="Direct pixel measurement − truth")
        axes[column].plot(days, 3*sigma, "--", color="gray", label="±3σ (σ = 0.15 px)")
        axes[column].plot(days, -3*sigma, "--", color="gray")
        axes[column].axhline(0, color="black", lw=.7)
        axes[column].set(xlabel="Days from first exposure", ylabel=f"{coordinate.title()} residual [px]", xticks=days)
        axes[column].legend(fontsize=8)
    figure.suptitle("Four direct geometric-center observations for a future EKF baseline\n"
                   "Separate from image-derived centroids, which contain illumination/model bias")
    figure.savefig(directory / "direct_observations.png", dpi=160)
    plt.close(figure)
    hashes = {}
    for filename in [summary["truth_csv"], *summary["loaded_kernels"]]:
        with Path(filename).open("rb") as stream:
            hashes[filename] = hashlib.file_digest(stream, "sha256").hexdigest()
    (directory / "diagnostics.json").write_text(json.dumps({"frames": diagnostics, "input_sha256": hashes}, indent=2)+"\n")
    print(json.dumps(diagnostics, indent=2))
    print(f"Saved sequence and direct-observation plots in {directory}")


if __name__ == "__main__":
    main()
