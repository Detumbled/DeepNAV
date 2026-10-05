"""Paired 8/9-state ensemble diagnostics; NumPy, Matplotlib and SciPy only."""

import argparse
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


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=Path("montecarlo/data"))
    parser.add_argument("--output", type=Path, default=Path("Output_montecarlo"))
    parser.add_argument("--no-show", action="store_true")
    args = parser.parse_args()
    if args.input.resolve() == args.output.resolve():
        raise ValueError("Input data and output plots must use different directories")
    study = json.loads((args.input / "study.json").read_text())
    data, summaries = {}, {}
    for case in ("matched", "srp"):
        for n in (8, 9):
            d = load(args.input / f"{case}_{n}.csv", study["runs"])
            data[case, n] = d
            summaries[case, n] = summarize(d)
        keys = ("run", "time_s", "tracking", "station", "srp_true",
                "observed_range_km", "observed_rate_km_s")
        if not all(np.array_equal(data[case, 8][key], data[case, 9][key], equal_nan=True)
                   for key in keys):
            raise ValueError(f"{case}: 8/9-state truth, observations or tracking differ")
    args.output.mkdir(parents=True, exist_ok=True)
    plots(data, summaries, study, args.output)
    print(f"Saved four 8/9-state comparison plots to {args.output}")
    if args.no_show:
        plt.close("all")
    else:
        plt.show()


if __name__ == "__main__":
    main()
