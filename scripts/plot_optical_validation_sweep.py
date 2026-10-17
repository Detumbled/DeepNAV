"""Plot the C++ phase/size/signal sweep and shared-nuisance sequences; no image fitting."""
import argparse
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
MODELS = ("matched", "sphere", "ellipsoid")
LABELS = {"matched": "Matched sphere control", "sphere": "Realistic → sphere", "ellipsoid": "Realistic → ellipsoid"}
COLORS = {"matched": "tab:blue", "sphere": "tab:orange", "ellipsoid": "tab:green"}
THRESHOLD_95 = 5.991464547


def read(path):
    return np.atleast_1d(np.genfromtxt(path, names=True, delimiter=",", dtype=None, encoding="utf-8"))


def matrices(data, prefix):
    result = np.empty((len(data), 2, 2))
    result[:, 0, 0] = data[prefix+"ss"]
    result[:, 0, 1] = result[:, 1, 0] = data[prefix+"sl"]
    result[:, 1, 1] = data[prefix+"ll"]
    return result


def squared_mahalanobis(errors, covariance):
    return np.einsum("ni,ni->n", errors, np.linalg.solve(covariance, errors[..., None])[..., 0])


def wilson(successes, total):
    if total == 0:
        return [None, None]
    z = 1.9599639845
    p = successes/total
    center = (p+z*z/(2*total))/(1+z*z/total)
    half = z*np.sqrt(p*(1-p)/total+z*z/(4*total*total))/(1+z*z/total)
    return [float(center-half), float(center+half)]


def metrics(data):
    good = data[data["converged"] == 1]
    total = len(data)
    result = {"attempted": total, "converged": len(good), "accepted": int(np.sum(data["accepted"])),
              "convergence_rate": len(good)/total if total else None,
              "acceptance_rate": float(np.mean(data["accepted"])) if total else None}
    if not len(good):
        return result
    extraction = np.column_stack((good["extract_dx"], good["extract_dy"]))
    navigation = np.column_stack((good["nav_dx"], good["nav_dy"]))
    result.update({"extraction_rms_pixels": float(np.sqrt(np.mean(np.sum(extraction**2, axis=1)))),
                   "navigation_rms_pixels": float(np.sqrt(np.mean(np.sum(navigation**2, axis=1)))),
                   "mean_extraction_error_pixels": extraction.mean(axis=0).tolist(),
                   "mean_navigation_error_pixels": navigation.mean(axis=0).tolist(),
                   "median_reduced_chi_squared": float(np.median(good["chi2"]))})
    norms = np.linalg.norm(extraction, axis=1)
    result["extraction_error_quantiles_pixels"] = {str(q): float(np.quantile(norms, q)) for q in (.5, .9, .95, 1.)}
    result["extraction_errors_over_one_pixel"] = int(np.sum(norms > 1))
    result["accepted_extraction_errors_over_one_pixel"] = int(np.sum((norms > 1) & (good["accepted"] == 1)))
    for name, errors, prefix in (("extraction_formal", extraction, "formal_"),
                                  ("navigation_formal", navigation, "formal_"),
                                  ("navigation_effective", navigation, "effective_")):
        distances = squared_mahalanobis(errors, matrices(good, prefix))
        covered = distances <= THRESHOLD_95
        successes = int(np.sum(covered))
        result[name+"_95_coverage"] = float(np.mean(covered))
        result[name+"_wilson_95"] = wilson(successes, len(good))
        result[name+"_covered_and_converged_fraction"] = successes/total
        result[name+"_mean_normalized_squared_error"] = float(np.mean(distances))
    accepted = good[good["accepted"] == 1]
    if len(accepted):
        accepted_errors = np.column_stack((accepted["nav_dx"], accepted["nav_dy"]))
        result["accepted_navigation_effective_95_coverage"] = float(np.mean(
            squared_mahalanobis(accepted_errors, matrices(accepted, "effective_")) <= THRESHOLD_95))
        result["accepted_extraction_rms_pixels"] = float(np.sqrt(np.mean(accepted["extract_dx"]**2+accepted["extract_dy"]**2)))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path, nargs="?", default=ROOT/"results/cassini_synthetic_optical/validation")
    args = parser.parse_args()
    directory = args.directory.resolve()
    config = json.loads((directory/"config.json").read_text())
    data = read(directory/"observations.csv")
    grid = data[data["experiment"] == "grid"]
    sequences = data[data["experiment"] == "sequence"]
    expected = 27*config["draws_per_cell"]*3+config["tracks"]*4*2
    if len(data) != expected:
        raise ValueError(f"Incomplete run: {len(data)} fit records, expected {expected}")
    if not np.allclose(grid["phase_actual"], grid["phase_requested"], atol=1e-4, rtol=0):
        raise ValueError("Rendered phase does not match the controlled sweep")
    if not np.allclose(grid["radius_actual"], grid["radius_requested"], atol=0, rtol=1e-3):
        raise ValueError("Apparent-size sweep was not realized in the rendered geometry")
    if len(np.unique(sequences["tdb"])) != 4:
        raise ValueError("Shared-error sequences must contain four distinct epochs")
    if len(grid) != len(np.unique(grid["scene"]))*3 or len(sequences) != len(np.unique(sequences["scene"]))*2:
        raise ValueError("Missing or duplicated paired-model fitting records")
    nuisances = read(directory/"nuisances.csv")
    if not np.array_equal(np.sort(nuisances["scene"]), np.unique(data["scene"])):
        raise ValueError("Nuisance ledger must contain exactly one row per scene")
    for scene in np.unique(data["scene"]):
        records = data[data["scene"] == scene]
        expected_models = MODELS if records["experiment"][0] == "grid" else ("sphere", "ellipsoid")
        if sorted(records["model"].tolist()) != sorted(expected_models):
            raise ValueError("A scene is missing a paired-model comparison")
    for track in range(config["tracks"]):
        scenes = np.unique(sequences[sequences["track"] == track]["scene"])
        ledger = nuisances[np.isin(nuisances["scene"], scenes)]
        for name in ("principal_x", "principal_y", "focal", "k1", "psf_x", "psf_y", "albedo", "ls_mix", "gain", "read", "bias"):
            if np.ptp(ledger[name]) > 1e-12:
                raise ValueError(f"Expected shared sequence nuisance changed: {name}")
    grid_ledger = nuisances[np.isin(nuisances["scene"], grid["scene"])]
    if any(np.ptp(grid_ledger[name]) == 0 for name in ("principal_x", "focal", "psf_x", "albedo", "pointing_x")):
        raise ValueError("Independent nuisance draws were not varied in the grid")
    # Fail loudly on corrupt successful covariances; never silently drop them.
    for prefix in ("formal_", "effective_"):
        covariance = matrices(data[data["converged"] == 1], prefix)
        if not np.isfinite(covariance).all() or np.any(np.linalg.eigvalsh(covariance) <= 0):
            raise ValueError("Non-positive or nonfinite converged covariance")
    phases, radii, signals = (config[k] for k in ("phases_degrees", "radii_pixels", "full_phase_electrons"))
    summary = {"config": config, "grid": {}, "cells": [], "sequences": {},
               "coverage_note": "Per-image coverage is conditional on convergence. Also report covered-and-converged fractions over all attempts. Wilson intervals are descriptive binomial approximations; cells have few draws."}
    for model in MODELS:
        summary["grid"][model] = metrics(grid[grid["model"] == model])
    fig, axes = plt.subplots(3, 3, figsize=(12, 10), constrained_layout=True)
    cmap = plt.get_cmap("viridis").copy()
    cmap.set_bad("#dddddd")
    for row, phase in enumerate(phases):
        for col, model in enumerate(MODELS):
            errors = np.full((3, 3), np.nan)
            annotations = {}
            for i, radius in enumerate(radii):
                for j, signal in enumerate(signals):
                    records = grid[(grid["model"] == model) & (grid["phase_requested"] == phase) &
                                   (grid["radius_requested"] == radius) & (grid["budget"] == signal)]
                    m = metrics(records)
                    summary["cells"].append({"model": model, "phase_degrees": phase, "radius_pixels": radius, "full_phase_electrons": signal, **m})
                    errors[i, j] = m.get("extraction_rms_pixels", np.nan)
                    annotations[i, j] = f"{errors[i,j]:.3f} px\n{m['converged']}/{m['attempted']} converged" if m["converged"] else "No converged fit"
            image = axes[row, col].imshow(errors, norm=LogNorm(.003, 3), cmap=cmap, origin="lower", aspect="auto")
            for (i, j), text in annotations.items():
                axes[row, col].text(j, i, text, ha="center", va="center", fontsize=8,
                                    color="white" if np.isfinite(errors[i,j]) and errors[i,j] < .15 else "black")
            axes[row, col].set(xticks=range(3), xticklabels=[f"{v:.0e}" for v in signals],
                               yticks=range(3), yticklabels=radii, xlabel="Full-phase photon budget [e−]",
                               ylabel="Requested radius [px]", title=f"{LABELS[model]} · phase {phase}°")
    fig.colorbar(image, ax=axes, label="Extraction RMS among converged fits [px] · color clipped to 0.003–3", shrink=.8)
    fig.suptitle(f"Independent phase × size × signal sweep · {config['draws_per_cell']} nuisance/noise draws per cell\n"
                 "Controlled Sun and focal scale; four SPK/PCK epochs with simulated observer truth")
    fig.savefig(directory/"sweep_accuracy.png", dpi=160); plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(14, 5), constrained_layout=True)
    x = np.arange(3)
    for j, model in enumerate(MODELS):
        records = grid[grid["model"] == model]
        groups = [metrics(records[records["phase_requested"] == p]) for p in phases]
        positions = x+(j-1)*.24
        axes[0].bar(positions, [g["convergence_rate"]*100 for g in groups], .23, color=COLORS[model], label=LABELS[model])
        axes[1].bar(positions, [g["acceptance_rate"]*100 for g in groups], .23, color=COLORS[model])
        values = np.array([g.get("navigation_effective_95_coverage", np.nan) for g in groups])*100
        interval = np.array([g.get("navigation_effective_wilson_95", [np.nan, np.nan]) for g in groups], float)*100
        axes[2].errorbar(positions, values, yerr=np.vstack((values-interval[:, 0], interval[:, 1]-values)),
                         fmt="o", color=COLORS[model], capsize=3)
        axes[2].scatter(positions, [g.get("navigation_formal_95_coverage", np.nan)*100 for g in groups],
                         marker="x", color=COLORS[model], alpha=.7)
    for ax, title in zip(axes, ("Converged / attempted", "Pass convergence + χ² checks / attempted", "95% navigation coverage · converged fits")):
        ax.set(xticks=x, xticklabels=[f"{p}°" for p in phases], xlabel="Phase", ylabel="Fraction [%]", ylim=(-2, 103), title=title)
        ax.grid(axis="y", alpha=.25)
    axes[0].legend(fontsize=8); axes[2].axhline(95, color="black", ls="--", lw=1)
    axes[2].plot([], [], "ko", label="With nuisance budgets")
    axes[2].plot([], [], "kx", label="Formal fit covariance only")
    axes[2].legend(fontsize=7, loc="lower right")
    fig.suptitle("Fit reliability and declared-covariance coverage\nCoverage uses formal covariance + declared pointing/calibration budgets; those budgets are not retuned")
    fig.savefig(directory/"sweep_reliability.png", dpi=160); plt.close(fig)

    # Averaging residuals is an uncertainty diagnostic, not an orbit solution.
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)
    for col, model in enumerate(("sphere", "ellipsoid")):
        records = sequences[sequences["model"] == model]
        summary["sequences"][model] = metrics(records)
        errors, covariances = [], []
        accepted_tracks = 0
        for track in range(config["tracks"]):
            entries = records[records["track"] == track]
            entries = np.sort(entries, order="day")
            if len(entries) != 4 or not np.all(entries["converged"]):
                continue
            errors.append(np.column_stack((entries["nav_dx"], entries["nav_dy"])))
            covariances.append(matrices(entries, "effective_"))
            accepted_tracks += int(np.all(entries["accepted"]))
        errors, covariances = np.array(errors), np.array(covariances)
        if not len(errors):
            summary["sequences"][model]["four_image_mean"] = {"complete_converged_tracks": 0}
            for ax in axes[:, col]:
                ax.text(.5, .5, "No complete converged tracks", transform=ax.transAxes, ha="center")
            continue
        means = errors.mean(axis=1)
        independent = covariances.sum(axis=1)/16
        normalized = squared_mahalanobis(means, independent)
        covered = int(np.sum(normalized <= THRESHOLD_95))
        correlations = [np.corrcoef(errors[:, :, axis], rowvar=False) for axis in range(2)] if len(errors)>1 else None
        average_correlation = [float(c[np.triu_indices(4, 1)].mean()) for c in correlations] if correlations else None
        empirical = np.cov(means.T) if len(errors)>1 else None
        predicted = independent.mean(axis=0)
        summary["sequences"][model]["four_image_mean"] = {
            "attempted_tracks": config["tracks"], "complete_converged_tracks": len(means), "all_four_accepted_tracks": accepted_tracks,
            "mean_error_pixels": means.mean(axis=0).tolist(), "rms_error_pixels": float(np.sqrt(np.mean(np.sum(means**2, axis=1)))),
            "naive_independent_95_coverage": covered/len(means), "wilson_95": wilson(covered, len(means)),
            "covered_complete_fraction_all_tracks": covered/config["tracks"],
            "mean_normalized_squared_error": float(normalized.mean()), "empirical_mean_covariance": empirical.tolist() if empirical is not None else None,
            "mean_naive_predicted_covariance": predicted.tolist(), "empirical_to_predicted_covariance_trace_ratio": float(np.trace(empirical)/np.trace(predicted)) if empirical is not None else None,
            "mean_cross_day_correlation_sample_line": average_correlation}
        for track in errors:
            axes[0, col].plot(range(4), track[:, 0], alpha=.3, lw=1, color=COLORS[model])
        axes[0, col].axhline(0, color="black", lw=.7)
        axes[0, col].set(xticks=range(4), xticklabels=["Oct 10", "Oct 11", "Oct 12", "Oct 13"],
                         ylabel="Sample navigation error [px]", title=f"{LABELS[model]} · {len(means)} complete tracks")
        axes[1, col].hist(normalized, bins=12, color=COLORS[model], alpha=.75)
        axes[1, col].axvline(THRESHOLD_95, color="black", ls="--", label="Nominal 95% boundary (2D)")
        axes[1, col].set(xlabel="Four-image mean errorᵀ (ΣP / 16)⁻¹ error", ylabel="Tracks",
                         title=f"Independent-covariance coverage: {covered}/{len(means)} = {covered/len(means):.0%}")
        axes[1, col].legend(fontsize=8)
    for ax in axes.flat:
        ax.grid(alpha=.2)
    fig.suptitle("Shared calibration, PSF and alignment across four days\nSame nominal per-image budgets; averaging residuals does not remove persistent errors")
    fig.savefig(directory/"sequence_correlations.png", dpi=160); plt.close(fig)

    fig, axes = plt.subplots(3, 6, figsize=(16, 8), constrained_layout=True)
    for row, phase in enumerate(phases):
        for col, signal in enumerate(signals):
            record = grid[(grid["model"] == "ellipsoid") & (grid["draw"] == 0) &
                          (grid["phase_requested"] == phase) & (grid["radius_requested"] == 6) & (grid["budget"] == signal)][0]
            pixels = read(directory/f"scene_{record['scene']}.csv")
            image = axes[row, 2*col].imshow(np.maximum(pixels["dn"].reshape(64, 64)-100, .05),
                    cmap="inferno", norm=LogNorm(.5, 5000), origin="upper", extent=(-.5, 63.5, 63.5, -.5))
            residual = axes[row, 2*col+1].imshow(pixels["residual_sigma"].reshape(64, 64),
                    cmap="RdBu_r", vmin=-5, vmax=5, origin="upper", extent=(-.5, 63.5, 63.5, -.5))
            axes[row, 2*col].set_title(f"{phase}° · {signal:.0e} e−\nObserved DN − nominal bias", fontsize=9)
            axes[row, 2*col+1].set_title(f"Ellipsoid fit residual / σ\nχ²={record['chi2']:.2f}", fontsize=9)
            for ax in axes[row, 2*col:2*col+2]:
                ax.set(xlim=(24, 51), ylim=(48, 20), xlabel="Sample", ylabel="Line")
    fig.colorbar(image, ax=axes[:, ::2], shrink=.7, label="DN above nominal bias (shared log scale)")
    fig.colorbar(residual, ax=axes[:, 1::2], shrink=.7, label="Residual / inferred σ (shared scale)")
    fig.suptitle("Representative realistic images · requested radius 6 px\nSynthetic phase and signal vary independently; noise, PSF, albedo and pointing are resampled")
    fig.savefig(directory/"sweep_images.png", dpi=160); plt.close(fig)

    checks = read(directory/"resolution.csv")
    def scalar(value):
        value = value.item()
        return None if isinstance(value, float) and not np.isfinite(value) else value
    summary["resolution_checks"] = [{k: scalar(row[k]) for k in checks.dtype.names} for row in checks]
    summary["sha256"] = {}
    for name in (config["truth_csv"], *config["loaded_kernels"], str(directory/"observations.csv"), str(directory/"nuisances.csv"),
                 str(ROOT/"benchmarks/cassini_synthetic_optical/validation.cpp"), str(ROOT/"src/opnav/estimation/BodyRenderer.cpp"),
                 str(ROOT/"src/opnav/estimation/BodyEstimator.cpp"), str(Path(__file__).resolve())):
        with Path(name).open("rb") as stream:
            summary["sha256"][name] = hashlib.file_digest(stream, "sha256").hexdigest()
    (directory/"metrics.json").write_text(json.dumps(summary, indent=2, allow_nan=False)+"\n")
    print(json.dumps({"grid": summary["grid"], "sequences": summary["sequences"], "resolution": summary["resolution_checks"]}, indent=2))
    print(f"Saved four plots and metrics.json to {directory}")


if __name__ == "__main__":
    main()
