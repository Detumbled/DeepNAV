"""Evaluate C++ rejection decisions against evaluation-only truth; never select thresholds here."""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from plot_optical_validation_sweep import read, metrics, MODELS, LABELS

ROOT = Path(__file__).resolve().parents[1]


def summarize(data, gates):
    old = data["accepted"] == 1
    keep = gates["accepted"] == 1
    errors = np.hypot(data["extract_dx"], data["extract_dy"])
    if np.any(keep & ~old) or not np.isfinite(errors[old]).all():
        raise ValueError("Invalid retained or baseline-accepted fit")
    bad = old & (errors > 1.)
    good = old & (errors <= 1.)
    reasons = Counter(reason for row in gates[~keep]["reasons"] for reason in str(row).split(";") if reason)
    return {"attempted": len(data), "old_accepted": int(old.sum()), "retained": int(keep.sum()),
            "retained_fraction_all_attempts": float(keep.mean()),
            "bad_reference_threshold_pixels": 1., "previously_accepted_bad": int(bad.sum()),
            "caught_bad": int(np.sum(bad & ~keep)), "remaining_accepted_bad": int(np.sum(bad & keep)),
            "previously_accepted_at_most_one_pixel": int(good.sum()),
            "rejected_at_most_one_pixel": int(np.sum(good & ~keep)),
            "retention_among_at_most_one_pixel": float(np.sum(good & keep)/good.sum()) if good.sum() else None,
            "rejection_reasons_overlapping": dict(reasons),
            "before": metrics(data[old]) if old.sum() else None,
            "after": metrics(data[keep]) if keep.sum() else None}


def load(directory):
    config = json.loads((directory/"config.json").read_text())
    data, gates = read(directory/"observations.csv"), read(directory/"rejections.csv")
    expected = 27*config["draws_per_cell"]*3+config["tracks"]*8
    if len(data) != expected or len(gates) != expected:
        raise ValueError(f"Incomplete rejection run in {directory}")
    if not np.array_equal(data["scene"], gates["scene"]) or not np.array_equal(data["model"], gates["model"]):
        raise ValueError("Rejection ledger is not aligned with observation records")
    if np.any(gates["ambiguity_check"] != 0) or np.any(gates["navigation_check"] != 0):
        raise ValueError("This experiment must use image-only rejection")
    return config, data, gates


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path, nargs="?", default=ROOT/"results/cassini_synthetic_optical/validation_rejection")
    parser.add_argument("--holdout", type=Path)
    args = parser.parse_args()
    directories = [args.directory.resolve()]
    if args.holdout:
        directories.append(args.holdout.resolve())
    runs = [load(path) for path in directories]
    thresholds = [json.loads((path/"rejection_options.json").read_text()) for path in directories]
    if len(runs)>1 and (thresholds[0] != thresholds[1] or runs[0][0]["seed"] == runs[1][0]["seed"]):
        raise ValueError("Holdout must use identical thresholds and a distinct random seed")
    output = directories[0]
    summary = {"thresholds": thresholds[0], "runs": [],
               "note": "Truth is used only for scoring. >1 px is a reference failure threshold, not a flight accuracy requirement. Reasons overlap. Shared calibration uncertainty is not corrected by selection."}
    labels = ["Original sweep", "Independent holdout"][:len(runs)]
    baseline = ROOT/"results/cassini_synthetic_optical/validation/observations.csv"
    if baseline.exists() and runs[0][0]["seed"] == 20261017:
        old = read(baseline); current = runs[0][1]
        if len(old) == len(current):
            for name in old.dtype.names:
                if old[name].dtype.kind in "f":
                    if not np.allclose(old[name], current[name], atol=1e-8, rtol=1e-10, equal_nan=True):
                        raise ValueError(f"Baseline fit changed: {name}")
                elif not np.array_equal(old[name], current[name]):
                    raise ValueError(f"Baseline fit changed: {name}")
            summary["original_fits_unchanged"] = True
    for label, (config, data, gates) in zip(labels, runs):
        result = {"label": label, "config": config, "grid": {}, "sequence": {}, "cells": []}
        for experiment in ("grid", "sequence"):
            for model in MODELS:
                selected = (data["experiment"] == experiment) & (data["model"] == model)
                if np.any(selected):
                    result[experiment][model] = summarize(data[selected], gates[selected])
        for phase in config["phases_degrees"]:
            for radius in config["radii_pixels"]:
                for budget in config["full_phase_electrons"]:
                    for model in MODELS:
                        selected = (data["experiment"] == "grid") & (data["phase_requested"] == phase) & (data["radius_requested"] == radius) & (data["budget"] == budget) & (data["model"] == model)
                        result["cells"].append({"phase": phase, "radius": radius, "budget": budget, "model": model,
                                                **summarize(data[selected], gates[selected])})
        summary["runs"].append(result)

    fig, axes = plt.subplots(len(runs), 3, figsize=(13, 4*len(runs)), squeeze=False, constrained_layout=True)
    x = np.arange(3)
    for row, result in enumerate(summary["runs"]):
        groups = [result["grid"][model] for model in MODELS]
        before = np.array([g["old_accepted"]/g["attempted"] for g in groups])*100
        after = np.array([g["retained_fraction_all_attempts"] for g in groups])*100
        axes[row, 0].bar(x-.18, before, .36, label="Previous fit checks")
        axes[row, 0].bar(x+.18, after, .36, label="Additional image-only rejection")
        axes[row, 0].set(ylabel="Accepted / attempted [%]", ylim=(0, 110), title=result["label"]+" · retention")
        axes[row, 1].bar(x-.18, [g["previously_accepted_bad"] for g in groups], .36)
        axes[row, 1].bar(x+.18, [g["remaining_accepted_bad"] for g in groups], .36)
        axes[row, 1].set(ylabel="Accepted center errors >1 px", title="Reference failures surviving selection")
        for i, g in enumerate(groups):
            axes[row, 1].text(i-.18, g["previously_accepted_bad"]+.1, str(g["previously_accepted_bad"]), ha="center", fontsize=9)
            axes[row, 1].text(i+.18, g["remaining_accepted_bad"]+.1, str(g["remaining_accepted_bad"]), ha="center", fontsize=9)
        rms_before = [g["before"]["extraction_rms_pixels"] if g["before"] else np.nan for g in groups]
        rms_after = [g["after"]["extraction_rms_pixels"] if g["after"] else np.nan for g in groups]
        axes[row, 2].bar(x-.18, rms_before, .36)
        axes[row, 2].bar(x+.18, rms_after, .36)
        maximum = float(np.nanmax([rms_before, rms_after]))
        for i in range(3):
            for offset, value in ((-.18, rms_before[i]), (.18, rms_after[i])):
                if np.isfinite(value):
                    axes[row, 2].text(i+offset, value+.02*maximum, f"{value:.3f}", ha="center", fontsize=8)
        axes[row, 2].set(ylabel="Accepted extraction RMS [px]", ylim=(0, maximum*1.18), title="Accuracy of retained observations")
        for ax in axes[row]:
            ax.set(xticks=x, xticklabels=["Matched", "Realistic sphere", "Realistic ellipsoid"])
            ax.tick_params(axis="x", labelrotation=15); ax.grid(axis="y", alpha=.2)
    axes[0, 0].legend(fontsize=8)
    fig.suptitle("Fixed rejection thresholds · original sweep and independent nuisance/noise draws\nNo navigation prediction, truth center or extra fitting is used for selection")
    fig.savefig(output/"rejection_summary.png", dpi=160); plt.close(fig)

    config, data, gates = runs[0]
    fig, axes = plt.subplots(3, 3, figsize=(12, 9), constrained_layout=True)
    for row, phase in enumerate(config["phases_degrees"]):
        for col, model in enumerate(MODELS):
            fractions = np.zeros((3, 3))
            for i, radius in enumerate(config["radii_pixels"]):
                for j, budget in enumerate(config["full_phase_electrons"]):
                    selected = (data["experiment"] == "grid") & (data["model"] == model) & (data["phase_requested"] == phase) & (data["radius_requested"] == radius) & (data["budget"] == budget)
                    m = summarize(data[selected], gates[selected]); fractions[i, j] = m["retained_fraction_all_attempts"]*100
                    axes[row, col].text(j, i, f"{m['retained']}/{m['attempted']} retained\n{m['remaining_accepted_bad']} error >1 px", ha="center", va="center", fontsize=8,
                                        color="white" if fractions[i,j]<40 else "black")
            image = axes[row, col].imshow(fractions, origin="lower", vmin=0, vmax=100, cmap="viridis", aspect="auto")
            axes[row, col].set(title=f"{LABELS[model]} · {phase}°", xticks=range(3), yticks=range(3),
                               xticklabels=[f"{v:.0e}" for v in config["full_phase_electrons"]], yticklabels=config["radii_pixels"],
                               xlabel="Full-phase photon budget [e−]", ylabel="Radius [px]")
    fig.colorbar(image, ax=axes, label="Retained / all attempts [%]", shrink=.75)
    fig.suptitle("Original sweep · coverage lost through conservative image-quality selection\nTruth only scores remaining >1 px errors; rejected images are never silently removed from the denominator")
    fig.savefig(output/"rejection_retention.png", dpi=160); plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(13, 4), constrained_layout=True)
    for ax, model in zip(axes, MODELS):
        selected = (data["experiment"] == "grid") & (data["model"] == model) & (data["accepted"] == 1)
        records, decisions = data[selected], gates[selected]
        finite = np.isfinite(decisions["peak_snr"]) & np.isfinite(decisions["formal_sigma"])
        keep = decisions["accepted"] == 1
        for accepted, color, label in ((True, "tab:green", "Retained"), (False, "tab:red", "Rejected")):
            subset = finite & (keep == accepted)
            ax.scatter(decisions["peak_snr"][subset], decisions["formal_sigma"][subset], s=15, color=color, alpha=.5, label=label)
        bad = finite & (np.hypot(records["extract_dx"], records["extract_dy"]) > 1)
        ax.scatter(decisions["peak_snr"][bad], decisions["formal_sigma"][bad], s=75, facecolors="none", edgecolors="black", label=">1 px reference error")
        ax.axvline(thresholds[0]["minimum_peak_snr"], color="black", ls="--", lw=.8)
        ax.axhline(thresholds[0]["maximum_formal_sigma"], color="black", ls=":", lw=.8)
        ax.set(xscale="log", yscale="log", xlabel="Observed peak SNR", ylabel="Formal center σ, largest axis [px]",
               title=f"{LABELS[model]}\n{int(finite.sum())}/{len(records)} prior-accepted fits plotted")
        ax.grid(alpha=.2)
    axes[0].legend(fontsize=7)
    fig.suptitle("Diagnostics from the original sweep · additional source, masking and boundary gates also apply")
    fig.savefig(output/"rejection_diagnostics.png", dpi=160); plt.close(fig)

    summary["sha256"] = {}
    for directory in directories:
        for filename in ("observations.csv", "rejections.csv", "rejection_options.json", "config.json"):
            with (directory/filename).open("rb") as stream:
                summary["sha256"][str(directory/filename)] = hashlib.file_digest(stream, "sha256").hexdigest()
    for filename in ("include/opnav/estimation/EstimatorRejection.hpp", "src/opnav/estimation/EstimatorRejection.cpp", "scripts/plot_optical_rejection.py"):
        with (ROOT/filename).open("rb") as stream:
            summary["sha256"][filename] = hashlib.file_digest(stream, "sha256").hexdigest()
    (output/"rejection_metrics.json").write_text(json.dumps(summary, indent=2, allow_nan=False)+"\n")
    for run in summary["runs"]:
        print(run["label"])
        for model, m in run["grid"].items():
            print(model, "retained", m["retained"], "/", m["attempted"], "caught bad", m["caught_bad"], "/", m["previously_accepted_bad"], "remaining", m["remaining_accepted_bad"],
                  "RMS", m["after"]["extraction_rms_pixels"] if m["after"] else None)
    print("Saved rejection plots and metrics to", output)


if __name__ == "__main__":
    main()
