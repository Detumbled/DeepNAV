"""Plot C++ body-estimation robustness results; this script does not fit images."""
import argparse
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
LABELS = {"matched":"Matched", "psf":"PSF mismatch", "shape_sphere":"Ellipsoid → sphere",
          "shape_ellipsoid":"Ellipsoid → ellipsoid", "surface":"Surface mismatch", "calibration":"Calibration",
          "attitude":"Attitude", "exposure_midpoint":"Smear → midpoint", "exposure_integrated":"Smear → integrated",
          "combined_sphere":"Combined → sphere", "combined_ellipsoid":"Combined → ellipsoid"}


def matrix(record, prefix):
    return np.array([[record[prefix+"ss"], record[prefix+"sl"]], [record[prefix+"sl"], record[prefix+"ll"]]], float)


def ellipse(ax, covariance, color, label):
    values, vectors = np.linalg.eigh(covariance)
    if values.min() <= 0:
        return
    angle = np.degrees(np.arctan2(vectors[1, 1], vectors[0, 1]))
    patch = Ellipse((0, 0), 2*np.sqrt(5.991*values[1]), 2*np.sqrt(5.991*values[0]), angle=angle,
                    fill=False, edgecolor=color, linewidth=1.5, label=label)
    ax.add_patch(patch)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path, nargs="?", default=ROOT / "results/cassini_synthetic_optical/realism")
    args = parser.parse_args()
    directory = args.directory.resolve()
    if not (directory/"cases.csv").is_file():
        parser.error("Run ./build-clang/cassini_optical_realism first")
    cases = np.atleast_1d(np.genfromtxt(directory/"cases.csv", delimiter=",", names=True, dtype=None, encoding="utf-8"))
    trials = np.atleast_1d(np.genfromtxt(directory/"trials.csv", delimiter=",", names=True, dtype=None, encoding="utf-8"))
    config = json.loads((directory/"config.json").read_text())
    lookup = {row["case"]: row for row in cases}
    # Quantities are kept separate: source center extraction vs nominal EKF projection.
    extraction = np.hypot(cases["extract_dx"], cases["extract_dy"])
    navigation = np.hypot(cases["nav_dx"], cases["nav_dy"])
    figure, axes = plt.subplots(2, 1, figsize=(13, 9), constrained_layout=True)
    x = np.arange(len(cases)); labels = [LABELS[name] for name in cases["case"]]
    axes[0].bar(x-.2, extraction, .4, label="Fitted − actual image geometric center")
    axes[0].bar(x+.2, navigation, .4, label="Fitted − nominal navigation prediction")
    axes[0].set(yscale="log", ylabel="Center error [px]", xticks=x, xticklabels=labels,
                title="One common noisy realization per controlled scenario")
    axes[0].tick_params(axis="x", labelrotation=35); axes[0].legend(fontsize=9)
    axes[1].bar(x, cases["chi2"], color=["tab:red" if v else "tab:green" for v in cases["model_mismatch"]])
    axes[1].axhline(1, color="black", ls="--", label="Unit reduced χ²")
    axes[1].axhline(3, color="tab:red", ls=":", label="Configured mismatch flag threshold")
    axes[1].set(ylabel="Reduced χ²", xticks=x, xticklabels=labels, title="Small residuals alone do not establish an unbiased center")
    axes[1].tick_params(axis="x", labelrotation=35); axes[1].legend(fontsize=9)
    for ax in axes:
        ax.grid(axis="y", alpha=.25)
    figure.suptitle("C++ BodyEstimator · inferred pixel variances · known simulated Cassini trajectory\n"
                   "Illustrative truth/model differences, not flight-calibrated Cassini error budgets")
    figure.savefig(directory/"robustness_summary.png", dpi=160); plt.close(figure)

    selected = ["matched", "shape_sphere", "combined_sphere", "combined_ellipsoid"]
    pixels = {name: np.genfromtxt(directory/(name+".csv"), delimiter=",", names=True) for name in selected}
    maximum = max(float(data["dn"].max()) for data in pixels.values())
    figure, axes = plt.subplots(len(selected), 3, figsize=(12, 13), constrained_layout=True)
    for row, name in enumerate(selected):
        data = pixels[name]; record = lookup[name]
        width = len(np.unique(data["sample"])); height = len(np.unique(data["line"]))
        extent = (data["sample"].min()-.5, data["sample"].max()+.5, data["line"].max()+.5, data["line"].min()-.5)
        geometric = [data["truth_sample"][0], data["truth_line"][0]]
        nominal = [data["nominal_sample"][0], data["nominal_line"][0]]
        fitted = [record["center_sample"], record["center_line"]]
        normalized = data["residual_dn"]/np.sqrt(data["variance_estimated"])
        for col, (values, title) in enumerate(zip((data["dn"], data["model_dn"], normalized),
                ("Observed", "C++ fitted model", f"Residual / σ · χ²={record['chi2']:.2f}"))):
            raster = axes[row, col].imshow(values.reshape(height, width), extent=extent, origin="upper", interpolation="nearest",
                                           cmap="inferno" if col<2 else "RdBu_r", vmin=100 if col<2 else -5, vmax=maximum if col<2 else 5)
            axes[row, col].set_title(f"{LABELS[name]}\n{title}", fontsize=10)
            axes[row, col].set(xlim=(geometric[0]-13, geometric[0]+13), ylim=(geometric[1]+13, geometric[1]-13), xlabel="Sample [px]", ylabel="Line [px]")
            if col == 0:
                for point, marker, color, label in zip((geometric, nominal, fitted), ("+", "o", "x"), ("cyan", "lime", "white"),
                        ("Actual geometric center", "Nominal prediction", "C++ fit")):
                    axes[row, col].scatter(*point, marker=marker, s=60, color=color, label=label,
                                           facecolors="none" if marker=="o" else color)
                if row == 0: axes[row, col].legend(fontsize=6, loc="upper left")
            figure.colorbar(raster, ax=axes[row, col], label="DN" if col<2 else "σ", shrink=.75)
    figure.suptitle("Pixel reconstruction and residual structure · shared DN and residual scales")
    figure.savefig(directory/"image_comparison.png", dpi=160); plt.close(figure)

    figure, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
    summary = {"cases": [], "monte_carlo": {}}
    for row, ex, nav in zip(cases, extraction, navigation):
        summary["cases"].append({"case": row["case"], "extraction_error_pixels": float(ex), "navigation_error_pixels": float(nav),
                                  "reduced_chi_squared": float(row["chi2"]), "passes_fit_checks": bool(row["fit_checks_pass"])})
    for ax, name in zip(axes, ("matched", "combined_ellipsoid")):
        data = trials[trials["case"] == name]
        good = data[(data["status"] == 0) & np.isfinite(data["effective_ss"])]
        errors = np.column_stack((good["nav_dx"], good["nav_dy"]))
        extract = np.column_stack((good["extract_dx"], good["extract_dy"]))
        formal = np.mean([matrix(row, "formal_") for row in good], axis=0)
        effective = np.mean([matrix(row, "effective_") for row in good], axis=0)
        coverage = {"formal": [], "effective": []}
        for record, error in zip(good, errors):
            for kind in coverage:
                coverage[kind].append(error@np.linalg.solve(matrix(record, kind+"_"), error) <= 5.991)
        ax.scatter(*errors.T, s=18, alpha=.75, label="Repeated noisy-image center errors")
        ax.scatter(0, 0, marker="+", color="black", s=100, label="Nominal prediction")
        ellipse(ax, formal, "tab:orange", "Formal fit-only 95% ellipse")
        ellipse(ax, effective, "tab:green", "With declared nuisance budgets")
        limit = max(float(np.abs(errors).max())*1.2, float(np.sqrt(np.linalg.eigvalsh(effective).max()*5.991))*1.15, .02)
        ax.set(xlim=(-limit, limit), ylim=(-limit, limit), aspect="equal", xlabel="Sample error [px]", ylabel="Line error [px]",
               title=f"{LABELS[name]} · {len(good)}/{len(data)} converged\nFit-only coverage {np.mean(coverage['formal'])*100:.0f}%; with budgets {np.mean(coverage['effective'])*100:.0f}%")
        ax.grid(alpha=.2); ax.legend(fontsize=7, loc="lower left")
        summary["monte_carlo"][name] = {
            "trials": len(data), "converged": len(good), "mean_navigation_error_pixels": errors.mean(axis=0).tolist(),
            "mean_extraction_error_pixels": extract.mean(axis=0).tolist(), "extraction_rms_pixels": float(np.sqrt(np.mean(np.sum(extract**2, axis=1)))),
            "empirical_noise_covariance_pixels2": np.cov(errors.T).tolist(), "mean_formal_covariance_pixels2": formal.tolist(),
            "fit_only_95_coverage": float(np.mean(coverage["formal"])), "with_declared_budgets_95_coverage": float(np.mean(coverage["effective"]))}
    figure.suptitle("Conditional detector-noise Monte Carlo · fixed nuisance offsets\n"
                   "Pilot sample; nuisance budgets are declared synthetic priors, not calibrated measurement accuracy")
    figure.savefig(directory/"uncertainty_comparison.png", dpi=160); plt.close(figure)
    summary["sha256"] = {}
    for name in [config["truth_csv"], *config["loaded_kernels"]]:
        with Path(name).open("rb") as stream:
            summary["sha256"][name] = hashlib.file_digest(stream, "sha256").hexdigest()
    (directory/"metrics.json").write_text(json.dumps(summary, indent=2)+"\n")
    print(json.dumps(summary["monte_carlo"], indent=2))
    print(f"Saved three diagnostic plots and metrics.json in {directory}")


if __name__ == "__main__":
    main()
