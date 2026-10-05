"""Paired 8/9-state ensemble diagnostics; NumPy, Matplotlib and SciPy only."""

import argparse
import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import chi2

COLORS = {8: "#215b8f", 9: "#d47720"}


def load(path, runs):
    data = np.genfromtxt(path, delimiter=",", names=True)
    if data.size == 0 or data.size % runs:
        raise ValueError(f"Incomplete ensemble: {path}")
    data = data.reshape(runs, -1)
    if not np.all(data["run"] == np.arange(runs)[:, None]):
        raise ValueError(f"Unexpected run ordering: {path}")
    if not np.all(data["time_s"] == data["time_s"][0]):
        raise ValueError(f"Different time grids: {path}")
    return data


def mean_bounds(runs, dof):
    # Independent realizations at one epoch; no assumption of independence over time.
    return chi2.ppf([.025, .975], runs * dof) / runs


def summarize(data):
    summary = {}
    position = data["position_error_m"]
    summary["position_median"] = np.median(position, axis=0)
    summary["position_p05"], summary["position_p95"] = np.percentile(position, [5, 95], axis=0)
    summary["position_rmse"] = np.sqrt(np.mean(position**2, axis=0))
    summary["position_reference"] = 3 * np.median(data["position_sigma_max_m"], axis=0)
    summary["anees"] = np.mean(data["nees"], axis=0)
    summary["anees_orbit"] = np.mean(data["nees_orbit"], axis=0)
    tracking = data["tracking"][0].astype(bool)
    summary["anis"] = np.full(tracking.size, np.nan)
    summary["anis"][tracking] = np.mean(data["nis"][:, tracking], axis=0)
    for key, components in [("position", range(3)), ("velocity", range(3, 6)),
                            ("clock", range(6, 8))]:
        inside = np.stack([abs(data[f"error_{j}"]) <= 3 * data[f"sigma_{j}"] for j in components])
        summary[f"coverage_{key}"] = np.mean(inside, axis=(0, 1))
    summary["srp_rmse"] = np.sqrt(np.mean((data["srp_estimate"] - data["srp_true"])**2, axis=0))
    summary["srp_sigma_rms"] = np.sqrt(np.mean(data["srp_sigma"]**2, axis=0))
    return summary


def save(fig, axes, title, path):
    for axis in np.asarray(axes).flat:
        axis.grid(alpha=.22)
        axis.set_xlabel("Elapsed time [h]")
    fig.suptitle(title)
    fig.savefig(path, dpi=160)


def plots(data, summaries, study, output):
    runs = study["runs"]
    prior = study["initial_sigma"]
    density = study.get("acceleration_noise_density_m_s32", 0)
    label = (f"{runs} paired runs • initial σr={1000*prior[0]:g} m, σv={1000*prior[3]:g} m/s"
             f" • SNC √q={density:g} m/s^(3/2)")
    time = data["matched", 8]["time_s"][0] / 3600
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), layout="constrained")
    for col, case in enumerate(("matched", "srp")):
        for row, mask in enumerate((np.ones(time.size, dtype=bool), (time >= 2) & (time <= 4))):
            axis = axes[row, col]
            for n in (8, 9):
                s = summaries[case, n]
                axis.plot(time[mask], s["position_median"][mask], color=COLORS[n], label=f"{n} states: median")
                axis.fill_between(time[mask], s["position_p05"][mask], s["position_p95"][mask],
                                  color=COLORS[n], alpha=.16, label=f"{n} states: 5–95%")
                axis.plot(time[mask], s["position_reference"][mask], color=COLORS[n], ls="--",
                          label=f"{n} states: median 3√λmax(Prr)")
            axis.set_ylabel("Position [m]")
            axis.set_title(f"{case}: {'full arc' if row == 0 else 'first gap'}")
            if row == 0:
                axis.set_yscale("log")
                axis.legend(fontsize=7)
    save(fig, axes, label + "\nPosition-error distribution and uncertainty reference",
         output / "position_ensemble.png")

    fig, axes = plt.subplots(3, 2, figsize=(12, 10), layout="constrained", sharex=True)
    for col, case in enumerate(("matched", "srp")):
        for row, key in enumerate(("anees", "anees_orbit", "anis")):
            axis = axes[row, col]
            for n in (8, 9):
                dof = (9 if case == "srp" and n == 9 else 8) if row == 0 else (6 if row == 1 else 2)
                axis.plot(time, summaries[case, n][key] / dof, color=COLORS[n], label=f"{n} states")
                low, high = mean_bounds(runs, dof) / dof
                axis.fill_between(time, low, high, color=COLORS[n], alpha=.09)
            axis.axhline(1, color="black", ls="--", lw=1)
            axis.set_ylabel(["Mean NEES / state dof", "Mean orbital NEES / 6", "Mean NIS / 2"][row])
            axis.set_yscale("log")
            axis.legend(fontsize=8)
        axes[0, col].set_title(case)
    save(fig, axes, label + "\nEnsemble consistency • shaded: 95% pointwise χ² intervals • expected ratio: 1",
         output / "consistency.png")

    fig, axes = plt.subplots(3, 2, figsize=(12, 9), layout="constrained", sharex=True)
    expected = 99.7300204
    for col, case in enumerate(("matched", "srp")):
        for row, key in enumerate(("position", "velocity", "clock")):
            axis = axes[row, col]
            for n in (8, 9):
                axis.plot(time, 100 * summaries[case, n][f"coverage_{key}"],
                          color=COLORS[n], label=f"{n} states")
            axis.axhline(expected, color="black", ls="--", lw=1, label="Gaussian 3σ: 99.73%")
            axis.set_ylabel(f"{key.capitalize()} marginal coverage [%]")
            axis.set_ylim(0, 101)
            axis.legend(fontsize=8)
        axes[0, col].set_title(case)
    save(fig, axes, label + "\n±3σ marginal coverage • fraction over runs and components in each block",
         output / "coverage.png")

    fig, axes = plt.subplots(1, 2, figsize=(12, 4), layout="constrained")
    gap = (time >= 2) & (time <= 4)
    start = np.flatnonzero(gap)[0]
    for axis, case in zip(axes, ("matched", "srp")):
        for n in (8, 9):
            d = data[case, n]
            ratio = d["position_error_m"][:, gap] / d["position_error_m"][:, start, None]
            low, high = np.percentile(ratio, [5, 95], axis=0)
            axis.plot(time[gap], np.median(ratio, axis=0), color=COLORS[n], label=f"{n} states: median")
            axis.fill_between(time[gap], low, high, color=COLORS[n], alpha=.16, label=f"{n} states: 5–95%")
        axis.axhline(1, color="black", ls="--", lw=1)
        axis.set_title(case)
        axis.set_ylabel("Error norm / own error norm at 2 h")
        axis.legend(fontsize=8)
    save(fig, axes, label + "\nFirst-gap evolution • normalize each run separately • ratio 1: same error norm",
         output / "first_gap_growth.png")

    uncertain = data["srp", 9]
    fig, axes = plt.subplots(2, 1, figsize=(11, 7), layout="constrained")
    errors = uncertain["srp_estimate"] - uncertain["srp_true"]
    low, high = np.percentile(errors, [5, 95], axis=0)
    axes[0].plot(time, np.median(errors, axis=0), color=COLORS[9], label="Median error")
    axes[0].fill_between(time, low, high, color=COLORS[9], alpha=.2, label="Error 5–95%")
    sigma = np.median(uncertain["srp_sigma"], axis=0)
    axes[0].plot(time, 3 * sigma, color="black", ls="--", label="±3 median σ")
    axes[0].plot(time, -3 * sigma, color="black", ls="--")
    axes[0].set_ylabel("SRP-scale error")
    axes[1].plot(time, summaries["srp", 9]["srp_rmse"], color=COLORS[9], label="Empirical RMSE")
    axes[1].plot(time, summaries["srp", 9]["srp_sigma_rms"], color="black", ls="--", label="√mean(Pαα)")
    axes[1].set_ylabel("SRP-scale error / 1σ")
    for axis in axes:
        axis.legend(fontsize=8)
    save(fig, axes, label + "\n9-state SRP estimation • constant parameter, 20% prior • no parameter process noise",
         output / "srp_estimation.png")


def compare_baseline(data, summaries, study, baseline, output):
    reference = json.loads((baseline / "study.json").read_text())
    for key in ("runs", "seed", "start_utc", "step_s", "duration_s", "initial_sigma", "srp_sigma"):
        if reference[key] != study[key]:
            raise ValueError(f"Baseline differs in {key}")
    if reference.get("acceleration_noise_density_m_s32", 0) != 0:
        raise ValueError("Baseline must have zero orbital SNC noise")
    time = data["matched", 8]["time_s"][0] / 3600
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), layout="constrained", sharex=True)
    report = {}
    for col, case in enumerate(("matched", "srp")):
        for n in (8, 9):
            old = load(baseline / f"{case}_{n}.csv", study["runs"])
            for key in ("run", "time_s", "tracking", "station", "srp_true", "observed_range_km", "observed_rate_km_s"):
                if not np.array_equal(old[key], data[case, n][key], equal_nan=True):
                    raise ValueError(f"Baseline has different realizations: {case}/{n}/{key}")
            old_summary = summarize(old)
            dof = 9 if case == "srp" and n == 9 else 8
            for source, style, label in ((old_summary, ":", "Q orbit=0"),
                                         (summaries[case, n], "-", "SNC")):
                axes[0, col].plot(time, source["anees"] / dof, color=COLORS[n], ls=style,
                                  label=f"{n} states: {label}")
                axes[1, col].plot(time, source["position_median"], color=COLORS[n], ls=style,
                                  label=f"{n} states: {label}")
            low, high = mean_bounds(study["runs"], dof) / dof
            axes[0, col].fill_between(time, low, high, color=COLORS[n], alpha=.08)
            report[f"{case}_{n}"] = {
                "baseline_final_anees": float(old_summary["anees"][-1]),
                "snc_final_anees": float(summaries[case, n]["anees"][-1]),
                "baseline_final_position_median_m": float(old_summary["position_median"][-1]),
                "snc_final_position_median_m": float(summaries[case, n]["position_median"][-1]),
            }
        axes[0, col].axhline(1, color="black", ls="--", lw=1)
        axes[0, col].set_title(case)
        axes[0, col].set_ylabel("Mean NEES / state dof")
        axes[1, col].set_ylabel("Median position error [m]")
        for axis in axes[:, col]:
            axis.set_yscale("log")
            axis.legend(fontsize=8)
    save(fig, axes, f"{study['runs']} identical realizations • EKF with/without orbital SNC",
         output / "process_noise_comparison.png")
    (output / "process_noise_comparison.json").write_text(json.dumps(report, indent=2) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=Path("Output_montecarlo"))
    parser.add_argument("--no-show", action="store_true")
    parser.add_argument("--baseline", type=Path, help="Paired zero-SNC study for comparison")
    args = parser.parse_args()
    study = json.loads((args.input / "study.json").read_text())
    data, summaries, report = {}, {}, {}
    for case in ("matched", "srp"):
        for n in (8, 9):
            d = load(args.input / f"{case}_{n}.csv", study["runs"])
            data[case, n] = d
            summaries[case, n] = summarize(d)
            s = summaries[case, n]
            dof = 9 if case == "srp" and n == 9 else 8
            report[f"{case}_{n}"] = {
                "nees_dof": dof, "final_position_median_m": float(s["position_median"][-1]),
                "final_position_p95_m": float(s["position_p95"][-1]),
                "final_anees": float(s["anees"][-1]),
                "anees_95_pointwise_interval": mean_bounds(study["runs"], dof).tolist(),
                "mean_nis_over_tracking": float(np.mean(s["anis"][np.isfinite(s["anis"])])),
                "first_gap_end_start_ratio_median": float(np.median(
                    d["position_error_m"][:, d["time_s"][0] == 14400].ravel() /
                    d["position_error_m"][:, d["time_s"][0] == 7200].ravel())),
                "final_srp_rmse": float(s["srp_rmse"][-1]),
                "final_srp_sigma_rms": float(s["srp_sigma_rms"][-1]),
            }
        keys = ("run", "time_s", "tracking", "station", "srp_true", "observed_range_km", "observed_rate_km_s")
        if not all(np.array_equal(data[case, 8][key], data[case, 9][key], equal_nan=True) for key in keys):
            raise ValueError(f"{case}: 8/9-state truth, observations or tracking differ")
    plots(data, summaries, study, args.input)
    if args.baseline:
        compare_baseline(data, summaries, study, args.baseline, args.input)
    (args.input / "summary.json").write_text(json.dumps(report, indent=2) + "\n")
    keys = list(summaries["matched", 8])
    with (args.input / "ensemble.csv").open("w", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(["case", "states", "time_s", *keys])
        for (case, n), s in summaries.items():
            for k, time in enumerate(data[case, n]["time_s"][0]):
                writer.writerow([case, n, time, *[s[key][k] for key in keys]])
    print(json.dumps(report, indent=2))
    if args.no_show:
        plt.close("all")
    else:
        plt.show()


if __name__ == "__main__":
    main()
