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

    if n >= 8:
        fig, axes = plt.subplots(2, 1, figsize=(11, 6), sharex=True, layout="constrained")
        for axis, j, scale, label in zip(axes, [6, 7], [1e9, 1e11],
                                        ["Clock bias error [ns]", "Fractional-frequency error [×10⁻¹¹]"]):
            axis.plot(time, data[f"error_{j}"] * scale, color="#215b8f", label="Estimate − truth")
            sigma = data[f"sigma_{j}"] * scale
            axis.fill_between(time, -3 * sigma, 3 * sigma, color="#f4a261", alpha=.28, label="±3σ")
            axis.set_ylabel(label)
        shade_tracking(axes, data)
        axes[0].legend(fontsize=8)
        finish(fig, axes, f"{n}-state EKF • jointly estimated clock errors", output / f"ekf_{n}_clock_errors.png", scenario)


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


def model_comparison(baseline, n, directories, output, baseline_metadata):
    cases = [(baseline_metadata.get("scenario_label", "baseline"), baseline)]
    # A comparison must isolate estimator assumptions, with identical truth and measurements.
    keys = (["time_s", "tracking", "station", "observed_range_km", "observed_rate_km_s"]
            + [f"truth_{j}" for j in range(n)])
    for directory in directories:
        path = directory / f"ekf_{n}_gaps.csv"
        if not path.exists():
            raise ValueError(f"Missing comparison run: {path}")
        data = load(path)
        metadata = json.loads((directory / "scenario.json").read_text())
        if not all(key in baseline.dtype.names and key in data.dtype.names and
                   np.array_equal(baseline[key], data[key], equal_nan=True) for key in keys):
            raise ValueError(f"Comparison truth, epoch, tracking or observations differ: {directory}")
        cases.append((metadata.get("scenario_label", directory.name), data))

    tracking = baseline["tracking"].astype(bool)
    last_contact = baseline["time_s"][tracking][-1]
    fig, axes = plt.subplots(3 if n >= 8 else 1, 1, figsize=(11, 8 if n >= 8 else 4),
                             sharex=True, layout="constrained", squeeze=False)
    axes = axes[:, 0]
    for (label, data), color in zip(cases, plt.get_cmap("tab10").colors * 2):
        coast = data["time_s"] >= last_contact
        time = data["time_s"][coast] / 3600
        axes[0].plot(time, data["position_error_km"][coast] * 1000, color=color, label=label)
        axes[0].plot(time, data["position_bound_3sigma_km"][coast] * 1000, color=color, ls=":", alpha=.7)
        if n >= 8:
            for axis, j, scale in zip(axes[1:], [6, 7], [1e9, 1e14]):
                axis.plot(time, data[f"error_{j}"][coast] * scale, color=color)
                sigma = data[f"sigma_{j}"][coast] * scale
                axis.fill_between(time, -3 * sigma, 3 * sigma, color=color, alpha=.08)
    axes[0].set_ylabel("Position error / uncertainty [m]")
    axes[0].set_yscale("log")
    axes[0].legend(fontsize=8)
    axes[0].set_title("Solid: actual error norm • dotted: 3√λmax(Prr) uncertainty reference", fontsize=10)
    if n >= 8:
        axes[1].set_ylabel("Clock bias error [ns]")
        axes[2].set_ylabel("Frequency error [×10⁻¹⁴]")
        axes[1].set_title("Clock shading: each estimator's ±3σ • illustrative parameter mismatches", fontsize=10)
    finish(fig, axes, f"{n}-state EKF • model mismatch during the final tracking outage",
           output / f"ekf_{n}_model_comparison.png", baseline_metadata.get("force_model", ""))


def filter_comparison(eight, nine, output, scenario):
    keys = (["time_s", "tracking", "station", "observed_range_km", "observed_rate_km_s"]
            + [f"truth_{j}" for j in range(8)])
    if not all(np.array_equal(eight[k], nine[k], equal_nan=True) for k in keys):
        raise ValueError("8/9-state comparison requires identical truth, measurements and tracking")
    time = eight["time_s"] / 3600
    last = time[eight["tracking"].astype(bool)][-1]
    fig, axes = plt.subplots(3, 1, figsize=(11, 10), layout="constrained")
    masks = [np.ones(len(time), dtype=bool), (time >= 2) & (time < 4), time >= last]
    for axis, mask, title in zip(axes, masks, ["Full arc", "Scheduled gap: 2–4 h",
                                              f"After last contact: {last:.2f} h"]):
        for n, data, color in [(8, eight, "#215b8f"), (9, nine, "#d47720")]:
            axis.plot(time[mask], data["position_error_km"][mask] * 1000,
                      color=color, label=f"{n} states: error norm")
            axis.plot(time[mask], data["position_bound_3sigma_km"][mask] * 1000,
                      color=color, ls="--", label=f"{n} states: 3√λmax(Prr)")
        axis.set_title(title, fontsize=10)
        axis.set_ylabel("Position [m]")
    axes[0].set_yscale("log")
    shade_tracking(axes[0], eight)
    axes[0].legend(fontsize=8)
    finish(fig, axes, "8/9-state EKF • position error and uncertainty reference • tracking gaps",
           output / "ekf_8_9_position_accuracy.png", scenario)

    fig, axes = plt.subplots(3, 1, figsize=(11, 8), sharex=True, layout="constrained")
    for j, axis in enumerate(axes):
        for n, data, color in [(8, eight, "#215b8f"), (9, nine, "#d47720")]:
            axis.plot(time, data[f"error_{j}"] * 1000, color=color, label=f"{n} states: error")
            sigma = data[f"sigma_{j}"] * 1000
            axis.fill_between(time, -3 * sigma, 3 * sigma, color=color, alpha=.12,
                              label=f"{n} states: ±3σ")
        axis.set_ylabel(f"Position {'XYZ'[j]} [m]")
        axis.set_yscale("symlog", linthresh=10)
    shade_tracking(axes, eight)
    axes[0].legend(fontsize=8, ncol=2)
    finish(fig, axes, "8/9-state EKF • signed position errors and ±3σ • symmetric log scale",
           output / "ekf_8_9_position_errors.png", scenario)


def srp_scale(data, output, scenario):
    time = data["time_s"] / 3600
    fig, axis = plt.subplots(figsize=(11, 4), layout="constrained")
    axis.plot(time, data["truth_8"], color="black", ls="--", label="Truth")
    axis.plot(time, data["estimate_8"], color="#215b8f", label="Estimated SRP scale")
    sigma = data["sigma_8"]
    axis.fill_between(time, data["estimate_8"] - 3 * sigma, data["estimate_8"] + 3 * sigma,
                      color="#215b8f", alpha=.15, label="Estimate ±3σ")
    axis.set_ylabel("Effective SRP scale α [dimensionless]")
    axis.legend(fontsize=8)
    shade_tracking(axis, data)
    finish(fig, axis, "9-state EKF • constant SRP parameter • no parameter process noise",
           output / "ekf_9_srp_scale.png", scenario)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=Path("Output EKF"))
    parser.add_argument("--threshold-km", type=float, default=1.0)
    parser.add_argument("--no-show", action="store_true")
    parser.add_argument("--compare", type=Path, action="append", default=[],
                        help="Mismatch output folder; repeat to compare cases against --input")
    args = parser.parse_args()
    if not np.isfinite(args.threshold_km) or args.threshold_km <= 0:
        parser.error("--threshold-km must be positive and finite")
    found = False
    compared = False
    scenario = ""
    metadata = {}
    metadata_path = args.input / "scenario.json"
    if metadata_path.exists():
        metadata = json.loads(metadata_path.read_text())
        scenario = f"{metadata['integrator']} • {metadata['force_model']}"
        if metadata.get("scenario_label"):
            scenario += f" • {metadata['scenario_label']}"
    runs = {}
    # The six-state implementation remains available, but is retired from these figures.
    for path in args.input.glob("ekf_6_*.png"):
        path.unlink()
    for n in [8, 9]:
        path = args.input / f"ekf_{n}_gaps.csv"
        if not path.exists():
            continue
        found = True
        data = load(path)
        runs[n] = data
        continuous_path = args.input / f"ekf_{n}_continuous.csv"
        continuous = load(continuous_path) if continuous_path.exists() else None
        state_errors(data, n, args.input, scenario)
        position_accuracy(data, continuous, n, args.threshold_km, args.input, scenario)
        innovations(data, n, args.input, scenario)
        eclipse_flags(data, n, args.input, scenario)
        if n == 9:
            srp_scale(data, args.input, scenario)
        if args.compare and all((directory / f"ekf_{n}_gaps.csv").exists() for directory in args.compare):
            try:
                model_comparison(data, n, args.compare, args.input, metadata)
            except ValueError as error:
                parser.error(str(error))
            compared = True
        # Retire the superseded diagnostic when regenerating an existing output folder.
        (args.input / f"ekf_{n}_innovation_autocorrelation.png").unlink(missing_ok=True)
    if 8 in runs and 9 in runs:
        try:
            filter_comparison(runs[8], runs[9], args.input, scenario)
        except ValueError as error:
            parser.error(str(error))
    if not found:
        parser.error(f"No EKF gap-run CSVs in {args.input}; run ekf_demo first")
    if args.compare and not compared:
        parser.error("No state layout is shared by --input and every --compare folder")
    if args.no_show:
        plt.close("all")
    else:
        plt.show()


if __name__ == "__main__":
    main()
