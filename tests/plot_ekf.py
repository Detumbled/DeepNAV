"""Plot single-run EKF accuracy and innovations; ensemble studies belong in montecarlo/."""

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def load(path):
    data = np.genfromtxt(path, delimiter=",", names=True)
    if data.size < 2:
        raise ValueError(f"Expected at least two epochs in {path}")
    return data


def shade_tracking(axes, data):
    time = data["time_s"] / 3600
    tracking = data["tracking"].astype(bool)
    edges = np.r_[0, np.flatnonzero(np.diff(tracking)) + 1, len(time)]
    for start, end in zip(edges[:-1], edges[1:]):
        if tracking[start]:
            left = time[0] if start == 0 else (time[start - 1] + time[start]) / 2
            right = time[-1] if end == len(time) else (time[end - 1] + time[end]) / 2
            for axis in np.atleast_1d(axes).flat:
                axis.axvspan(left, right, color="#43aa8b", alpha=0.10, linewidth=0)


def finish(fig, axes, title, path, scenario):
    for axis in np.atleast_1d(axes).flat:
        axis.grid(alpha=0.22)
        axis.set_xlabel("Elapsed time [h]")
    fig.suptitle(title + ("\n" + scenario if scenario else ""))
    fig.savefig(path, dpi=170)


def state_errors(data, n, output, scenario):
    time = data["time_s"] / 3600
    fig, axes = plt.subplots(3, 2, figsize=(12, 8), sharex=True, layout="constrained")
    for j in range(6):
        axis = axes[j % 3, j // 3]
        scale = 1000 if j < 3 else 1e6
        axis.plot(time, data[f"error_{j}"] * scale, color="#215b8f", lw=1, label="Estimate − truth")
        sigma = data[f"sigma_{j}"] * scale
        axis.fill_between(time, -3 * sigma, 3 * sigma, color="#f4a261", alpha=0.28, label="±3σ")
        axis.set_ylabel(f"{'Position' if j < 3 else 'Velocity'} {'XYZ'[j % 3]} [{'m' if j < 3 else 'mm/s'}]")
    shade_tracking(axes, data)
    axes[0, 0].legend(loc="upper right", fontsize=8)
    finish(fig, axes, f"{n}-state EKF • tracking gaps • green bands: measurements available",
           output / f"ekf_{n}_state_errors.png", scenario)

    if n == 8:
        fig, axes = plt.subplots(2, 1, figsize=(11, 6), sharex=True, layout="constrained")
        for axis, j, scale, label in zip(axes, [6, 7], [1e9, 1e11],
                                        ["Clock bias error [ns]", "Fractional-frequency error [×10⁻¹¹]"]):
            axis.plot(time, data[f"error_{j}"] * scale, color="#215b8f", label="Estimate − truth")
            sigma = data[f"sigma_{j}"] * scale
            axis.fill_between(time, -3 * sigma, 3 * sigma, color="#f4a261", alpha=.28, label="±3σ")
            axis.set_ylabel(label)
        shade_tracking(axes, data)
        axes[0].legend(fontsize=8)
        finish(fig, axes, "8-state EKF • jointly estimated clock errors", output / "ekf_8_clock_errors.png", scenario)


def position_accuracy(data, continuous, n, threshold, output, scenario):
    fig, axis = plt.subplots(figsize=(11, 5), layout="constrained")
    time = data["time_s"] / 3600
    axis.plot(time, data["position_error_km"] * 1000, color="#215b8f", label="Error norm, tracking gaps")
    axis.plot(time, data["position_bound_3sigma_km"] * 1000, color="#e48b35", ls="--",
              label="3√λmax(Prr), gaps (uncertainty reference)")
    if continuous is not None:
        axis.plot(continuous["time_s"] / 3600, continuous["position_error_km"] * 1000,
                  color="#79858d", alpha=.8, label="Error norm, all visible opportunities")
    axis.axhline(threshold * 1000, color="#bd4141", ls=":", label=f"Example threshold: {threshold * 1000:g} m")
    shade_tracking(axis, data)
    axis.set_ylabel("Position [m]")
    axis.set_yscale("log")
    axis.legend(fontsize=8, loc="best")
    finish(fig, axis, f"{n}-state EKF • position accuracy through contact gaps",
           output / f"ekf_{n}_position_accuracy.png", scenario)


def innovations(data, n, output, scenario):
    time = data["time_s"] / 3600
    fig, axes = plt.subplots(2, 2, figsize=(12, 7), sharex=True, layout="constrained")
    for row, key, scale, label in [(0, "range", 1000, "Range innovation [m]"),
                                    (1, "rate", 1e6, "Range-rate innovation [mm/s]")]:
        unit = "km" if key == "range" else "km_s"
        axis = axes[row, 0]
        axis.plot(time, data[f"{key}_innovation_{unit}"] * scale, ".", ms=2.5, color="#215b8f")
        sigma = data[f"{key}_innovation_sigma_{unit}"] * scale
        axis.plot(time, 3 * sigma, color="#e48b35", lw=1, label="±3 innovation σ")
        axis.plot(time, -3 * sigma, color="#e48b35", lw=1)
        axis.set_yscale("symlog", linthresh=30 if key == "range" else .3)
        axis.set_ylabel(label + " (symmetric log)")
        axis = axes[row, 1]
        axis.plot(time, data[f"white_{key}"], ".", ms=2.5, color="#215b8f")
        axis.axhline(3, color="#e48b35", lw=1)
        axis.axhline(-3, color="#e48b35", lw=1)
        axis.set_ylabel(f"Whitened component {row + 1}")
    shade_tracking(axes, data)
    axes[0, 0].legend(fontsize=8)
    finish(fig, axes, f"{n}-state EKF • pre-update innovations (missing points: no tracking)",
           output / f"ekf_{n}_innovations.png", scenario)


def eclipse_flags(data, n, output, scenario):
    time = data["time_s"] / 3600
    fig, axes = plt.subplots(2, 1, figsize=(11, 5), sharex=True, layout="constrained")
    for source, label, color in [("truth", "Simulation truth", "#215b8f"),
                                  ("estimate", "Estimated geometry", "#e48b35")]:
        axes[0].plot(time, data[f"illumination_{source}"], color=color, lw=1.3,
                     ls="-" if source == "truth" else "--", label=label)
        axes[1].step(time, data[f"eclipse_{source}"], where="mid", color=color,
                     ls="-" if source == "truth" else "--", lw=1.3)
    axes[0].set_ylabel("Illuminated solar fraction")
    axes[0].set_ylim(-.05, 1.05)
    axes[0].legend(fontsize=8)
    axes[1].set_ylabel("Eclipse flag")
    axes[1].set_yticks([0, 1, 2], ["Sunlit", "Penumbra", "Umbra"])
    axes[1].set_ylim(-.2, 2.2)
    shade_tracking(axes, data)
    finish(fig, axes, f"{n}-state EKF • Earth eclipses • 60 s diagnostics sampling",
           output / f"ekf_{n}_eclipse_flags.png", scenario)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=Path("Output EKF"))
    parser.add_argument("--threshold-km", type=float, default=1.0)
    parser.add_argument("--no-show", action="store_true")
    args = parser.parse_args()
    if not np.isfinite(args.threshold_km) or args.threshold_km <= 0:
        parser.error("--threshold-km must be positive and finite")
    found = False
    scenario = ""
    metadata_path = args.input / "scenario.json"
    if metadata_path.exists():
        metadata = json.loads(metadata_path.read_text())
        scenario = f"{metadata['integrator']} • {metadata['force_model']}"
    for n in [6, 8]:
        path = args.input / f"ekf_{n}_gaps.csv"
        if not path.exists():
            continue
        found = True
        data = load(path)
        continuous_path = args.input / f"ekf_{n}_continuous.csv"
        continuous = load(continuous_path) if continuous_path.exists() else None
        state_errors(data, n, args.input, scenario)
        position_accuracy(data, continuous, n, args.threshold_km, args.input, scenario)
        innovations(data, n, args.input, scenario)
        eclipse_flags(data, n, args.input, scenario)
        # Retire the superseded diagnostic when regenerating an existing output folder.
        (args.input / f"ekf_{n}_innovation_autocorrelation.png").unlink(missing_ok=True)
    if not found:
        parser.error(f"No EKF gap-run CSVs in {args.input}; run ekf_demo first")
    if args.no_show:
        plt.close("all")
    else:
        plt.show()


if __name__ == "__main__":
    main()
