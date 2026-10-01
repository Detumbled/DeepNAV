#!/usr/bin/env python3
"""Plot clock history CSVs. Requires Pandas and Matplotlib; no repository imports.

Run python3 plot_clocks.py with no arguments to use the editable SETUP below.
Explicit command-line arguments are also supported and use their own defaults.
Display reduction affects curves only and may miss narrow peaks. Never use it
for Allan analysis or first threshold crossing calculations.
"""

import argparse
import csv
from pathlib import Path
import sys

try:
    import numpy as np
    import pandas as pd
except ImportError as error:
    raise SystemExit("Plotter dependencies are missing; install pandas matplotlib with pip.") from error


C = 299_792_458.0
REQUIRED = ("time_s", "bias_s", "fractional_frequency")
SIGMAS = ("sigma_bias_s", "sigma_fractional_frequency")

# ============================== SETUP ==============================
# Used when running: python3 plot_clocks.py
# Relative filenames are resolved inside output_directory, independently of
# the terminal's working directory. Generate the CSVs with test_clock_history.
SETUP = {
    "output_directory": Path(__file__).resolve().parent / "Output clocks",
    "csv_files": ["local.csv", "dsac.csv"],
    "mean_csv": ["local_mean.csv", "dsac_mean.csv"],  # None disables total budget.
    "allan_csv": ["local_allan.csv", "dsac_allan.csv"],  # None uses auto-discovery.
    "labels": ["USO (aging + white FM)", "DSAC (day-matched white FM)"],
    "time_unit": "days",  # "seconds", "hours", or "days"; Allan tau stays in s.
    "k": 3.0,
    "threshold_m": 1.0,  # None hides the range threshold.
    "max_points": 10000,  # Display only; analysis always uses complete histories.
    "q_bias_s": [2.5e-25, 7.776e-25],
    "drift_per_s": [1e-10 / 86400, 3e-16 / 86400],
    "uso_inset_index": 0,  # None disables the first-two-hours range inset.
    "output": "clocks.png",  # PNG, PDF or SVG; other figures add suffixes.
    "show": True,  # Open all four windows; False only saves the figures.
}
# ===================================================================


def setup_arguments():
    """Route SETUP through the same parser/validation as explicit CLI inputs."""
    directory = Path(SETUP["output_directory"])
    arguments = [str(directory / name) for name in SETUP["csv_files"]]
    for key, option in (("mean_csv", "--mean-csv"), ("allan_csv", "--allan-csv")):
        if SETUP[key] is not None:
            arguments.extend([option, *(str(directory / name) for name in SETUP[key])])
    if SETUP["labels"] is not None:
        arguments.extend(["--labels", *SETUP["labels"]])
    for key, option in (("time_unit", "--time-unit"), ("k", "--k"),
                        ("max_points", "--max-points")):
        arguments.extend([option, str(SETUP[key])])
    if SETUP["threshold_m"] is not None:
        arguments.extend(["--threshold-m", str(SETUP["threshold_m"])])
    arguments.extend(["--output", str(directory / SETUP["output"])])
    for key, option in (("q_bias_s", "--q-b"), ("drift_per_s", "--drift-per-s")):
        if SETUP[key] is not None:
            arguments.extend([option, *(str(value) for value in SETUP[key])])
    if SETUP["uso_inset_index"] is not None:
        arguments.extend(["--uso-inset-index", str(SETUP["uso_inset_index"])])
    if SETUP["show"]:
        arguments.append("--show")
    return arguments


def load_clock_csv(path):
    """Validate without sorting, filling missing values or aligning clocks."""
    path = Path(path)
    try:
        with path.open(newline="", encoding="utf-8-sig") as stream:
            header = next(csv.reader(stream), [])
        if len(header) != len(set(header)):
            raise ValueError("duplicate CSV column names")
        missing = set(REQUIRED) - set(header)
        if missing:
            raise ValueError("missing required columns: " + ", ".join(sorted(missing)))
        optional = set(SIGMAS) & set(header)
        if optional and optional != set(SIGMAS):
            raise ValueError("both sigma columns must be supplied together")
        if set(header) - set(REQUIRED + SIGMAS):
            raise ValueError("unexpected columns; expected the clock history CSV contract")
        frame = pd.read_csv(path, dtype="float64", engine="c", float_precision="round_trip",
                            encoding="utf-8-sig", skip_blank_lines=False)
        if not isinstance(frame.index, pd.RangeIndex):
            raise ValueError("data rows have more fields than the CSV header")
        if len(frame) < 2:
            raise ValueError("at least two data rows are required")
        if not np.isfinite(frame.to_numpy()).all():
            raise ValueError("missing or nonfinite numeric data")
        times = frame["time_s"]
        if (times < 0).any() or not (times.diff().iloc[1:] > 0).all():
            raise ValueError("timestamps must be nonnegative and strictly increasing")
        if optional and (frame[list(SIGMAS)] < 0).any().any():
            raise ValueError("sigma values must be nonnegative")
        return frame
    except (OSError, ValueError, pd.errors.ParserError) as error:
        raise ValueError(f"{path}: {error}") from error


def display_frame(frame, max_points):
    """Return positional display selection; leave the full DataFrame untouched."""
    if max_points == 0 or len(frame) <= max_points:
        return frame
    indices = np.linspace(0, len(frame) - 1, max_points, dtype=np.int64)
    return frame.iloc[indices]


def first_recorded_contact(frame, values, threshold):
    """First contact in the full recorded series, before display reduction.

    This does not claim to locate an unobserved crossing between samples.
    """
    indices = np.flatnonzero(np.abs(np.asarray(values)) >= threshold)
    if not len(indices):
        return None
    index = int(indices[0])
    previous = float(frame["time_s"].iloc[index - 1]) if index else None
    return index, float(frame["time_s"].iloc[index]), previous


def overlapping_allan_curve(frame):
    """Fractional-frequency ADEV from full uniform time-error data; drift retained.

    Limit tau to roughly one tenth of the record; long-tau estimates are noisy.
    """
    times = frame["time_s"].to_numpy()
    if len(times) < 3:
        raise ValueError("Allan deviation requires at least three samples")
    intervals = np.diff(times)
    dt = float(np.median(intervals))
    tolerance = 8 * abs(np.spacing(np.max(np.abs(times))))
    if dt <= 0 or not np.allclose(intervals, dt, rtol=1e-9, atol=tolerance):
        raise ValueError("Allan deviation requires uniformly sampled timestamps; no interpolation is applied")
    max_m = max(1, (len(times) - 1) // 10)
    factors = np.unique(np.geomspace(1, max_m, min(40, max_m)).astype(np.int64))
    bias = frame["bias_s"].to_numpy(dtype=np.longdouble)
    deviations = []
    for m in factors:
        tau = np.longdouble(m) * dt
        differences = bias[2*m:] - 2 * bias[m:-m] + bias[:-2*m]
        deviations.append(float(np.sqrt(np.mean((differences / tau)**2) / 2)))
    deviations = np.asarray(deviations)
    if not np.isfinite(deviations).all():
        raise ValueError("Allan deviation calculation overflowed")
    return factors * dt, deviations


def white_fm_adev(tau, q_bias_s, drift_per_s):
    """Simplified white FM plus retained constant aging, not a device spectrum fit."""
    tau = np.asarray(tau, dtype=float)
    return np.sqrt(q_bias_s / tau + 0.5 * (drift_per_s * tau)**2)


def analytical_budget_crossing(threshold, q_bias_s, drift_per_s, k=3.0):
    """Ideal zero-initial-state pointwise budget crossing; no first-passage claim."""
    if q_bias_s == 0 and drift_per_s == 0:
        return None
    def budget(t):
        return C * (0.5 * abs(drift_per_s) * t*t + k * np.sqrt(q_bias_s * t))
    if drift_per_s == 0 and (q_bias_s == 0 or k == 0):
        return None
    upper = 1.0
    while budget(upper) < threshold:
        upper *= 2
    lower = 0.0
    for _ in range(80):
        middle = (lower + upper) / 2
        if budget(middle) < threshold:
            lower = middle
        else:
            upper = middle
    return upper


def plot_histories(frames, labels, output, time_unit="hours", k=3.0,
                   threshold_m=None, max_points=10000, show=False, mean_frames=None,
                   allan_frames=None, q_bias_s=None, drift_per_s=None, uso_inset_index=None):
    import matplotlib
    if not show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    if not frames or len(frames) != len(labels):
        raise ValueError("provide one label per clock history")
    if allan_frames is not None and len(allan_frames) != len(frames):
        raise ValueError("provide one full Allan history per clock")
    allan_curves = [overlapping_allan_curve(frame) for frame in allan_frames] if allan_frames is not None else None
    if mean_frames is not None:
        if len(mean_frames) != len(frames):
            raise ValueError("provide one deterministic mean CSV per clock history")
        for frame, mean, label in zip(frames, mean_frames, labels):
            if "sigma_bias_s" not in frame:
                raise ValueError(f"{label}: total budget requires covariance")
            if not np.array_equal(frame["time_s"].to_numpy(), mean["time_s"].to_numpy()):
                raise ValueError(f"{label}: deterministic mean timestamps must match history exactly")

    if (q_bias_s is None) != (drift_per_s is None):
        raise ValueError("provide both q_b and drift for theoretical white-FM curves")
    if q_bias_s is not None:
        if len(q_bias_s) != len(frames) or len(drift_per_s) != len(frames):
            raise ValueError("provide one q_b and drift per clock")
        if not np.isfinite(q_bias_s).all() or min(q_bias_s) < 0 or not np.isfinite(drift_per_s).all():
            raise ValueError("q_b must be finite/nonnegative and drift must be finite")
    if uso_inset_index is not None and not 0 <= uso_inset_index < len(frames):
        raise ValueError("USO inset index must identify an input clock")
    horizon = max(float(frame["time_s"].iloc[-1]) for frame in frames)
    assumptions = ("Ideal initial bias/frequency calibration and ground reference; constant aging uncompensated; simplified white FM.\n"
                   "Flicker/other long-term noise omitted; simplified holdover comparison, not hardware validation.") if q_bias_s is not None else ""
    scale, unit = {"seconds": (1.0, "s"), "hours": (3600.0, "h"),
                   "days": (86400.0, "days")}[time_unit]
    output = Path(output)
    paths = [output,
             output.with_name(output.stem + "_uncertainty" + output.suffix),
             output.with_name(output.stem + "_diagnostics" + output.suffix)]
    if allan_curves is not None:
        paths.append(output.with_name(output.stem + "_allan" + output.suffix))
    figures = []

    def make_figure(rows, title):
        fig, axes = plt.subplots(rows, len(frames), squeeze=False,
                                 figsize=(6 * len(frames), 3.5 * rows),
                                 sharex="col", layout="constrained")
        figures.append(fig)
        fig.suptitle(title)
        if getattr(fig.canvas, "manager", None):
            fig.canvas.manager.set_window_title(title)
        return fig, axes

    def draw(ax, frame, values, label, ylabel, color, signed_threshold=False, envelope=False):
        shown = display_frame(frame, max_points)
        ax.plot(shown["time_s"] / scale, values.loc[shown.index], color=color, label=label)
        ax.set(xlabel=f"Elapsed time since calibration [{unit}]", ylabel=ylabel, xlim=(0, horizon / scale))
        ax.grid(True, alpha=0.25)
        if threshold_m is not None:
            ax.axhline(threshold_m, color="0.35", linestyle="--", label="Range threshold")
            if signed_threshold:
                ax.axhline(-threshold_m, color="0.35", linestyle="--")
            contact = first_recorded_contact(frame, values, threshold_m)
            if contact is None:
                status = f"No threshold crossing within {frame['time_s'].iloc[-1] / scale:.9g} {unit}"
            else:
                index, time, previous = contact
                ax.scatter([time / scale], [values.iloc[index]], color=color, marker="x", zorder=5)
                status = f"First sampled threshold crossing: {time / scale:.9f} {unit}"
                if previous is not None:
                    status += f"\nPrevious sample: {previous / scale:.9f} {unit}"
            if envelope:
                status = status.replace("First sampled threshold crossing", "Sampled stochastic envelope crossing")
            ax.text(0.02, 0.97, status, transform=ax.transAxes, va="top", fontsize=9,
                    bbox={"facecolor": "white", "alpha": 0.85, "edgecolor": "0.8"})
        ax.legend(loc="lower right", fontsize=9)

    try:
        primary, main_axes = make_figure(2, "Simplified clock holdover: one-way range error and pointwise budget")
        _, sigma_axes = make_figure(1, f"Stochastic range uncertainty ({k:g} sigma)")
        _, diagnostic_axes = make_figure(2, "Clock diagnostics: bias and fractional frequency")
        for column, (frame, label) in enumerate(zip(frames, labels)):
            color = f"C{column % 10}"
            bias_ns = 1e9 * frame["bias_s"]
            range_error = C * frame["bias_s"]
            uncertainty = k * C * frame["sigma_bias_s"] if "sigma_bias_s" in frame else None
            budget = None
            if mean_frames is not None:
                budget = C * mean_frames[column]["bias_s"].abs() + uncertainty
            for values in (bias_ns, range_error, uncertainty, budget):
                if values is not None and not np.isfinite(values.to_numpy()).all():
                    raise ValueError(f"{label}: derived plotting quantities overflow")
            draw(main_axes[0, column], frame, range_error, label,
                 "Signed one-way range error [m]", color, signed_threshold=True)
            main_axes[0, column].set_title(label + " — independent scale")
            if budget is not None:
                draw(main_axes[1, column], frame, budget, label,
                     f"Total budget c (|mean bias| + {k:g} sigma) [m]", color)
            else:
                main_axes[1, column].text(0.5, 0.5, "Total budget requires deterministic mean CSV\nand covariance; use --mean-csv",
                                               ha="center", va="center", transform=main_axes[1, column].transAxes)
                main_axes[1, column].set(xlabel=f"Elapsed time since calibration [{unit}]",
                                            ylabel="Total pointwise clock budget [m]")
            if uncertainty is not None:
                draw(sigma_axes[0, column], frame, uncertainty, label,
                     "Stochastic range uncertainty [m]", color, envelope=True)
            else:
                sigma_axes[0, column].text(0.5, 0.5, "No covariance supplied", ha="center", va="center",
                                          transform=sigma_axes[0, column].transAxes)
            sigma_axes[0, column].set_title(label + " — excludes deterministic bias")
            # Diagnostics deliberately have no range thresholds.
            shown = display_frame(frame, max_points)
            for row, (values, ylabel) in enumerate(((bias_ns, "Clock bias [ns]"),
                                                    (frame["fractional_frequency"], "Fractional frequency state y [1]"))):
                ax = diagnostic_axes[row, column]
                ax.plot(shown["time_s"] / scale, values.loc[shown.index], color=color, label=label)
                ax.set(xlabel=f"Elapsed time since calibration [{unit}]", ylabel=ylabel, xlim=(0, horizon / scale))
                ax.grid(True, alpha=0.25)
                ax.legend(loc="best")
            diagnostic_axes[0, column].set_title(label + " — independent scale")
            diagnostic_axes[1, column].set_title("Fractional frequency state — white FM excluded", fontsize=11)
            if budget is not None and q_bias_s is not None and threshold_m is not None:
                # The closed form applies only to these explicitly checked ideal initial conditions.
                mean = mean_frames[column]
                qb, drift = q_bias_s[column], drift_per_s[column]
                times = frame["time_s"].to_numpy()
                matches = (times[0] == 0 and
                           np.allclose(mean["bias_s"], 0.5 * drift * times**2, rtol=1e-7, atol=1e-22) and
                           np.allclose(frame["sigma_bias_s"], np.sqrt(qb * times), rtol=1e-7, atol=1e-22))
                if matches:
                    crossing = analytical_budget_crossing(threshold_m, qb, drift, k)
                    text = "Analytical budget: no crossing" if crossing is None else f"Analytical budget crossing: {crossing / scale:.8f} {unit}"
                    main_axes[1, column].text(0.02, 0.74, text, transform=main_axes[1, column].transAxes,
                        fontsize=9, bbox={"facecolor": "white", "alpha": 0.85, "edgecolor": "0.8"})
            if column == uso_inset_index:
                inset = main_axes[0, column].inset_axes([0.48, 0.32, 0.49, 0.37])
                early = frame[frame["time_s"] <= 7200]
                selected = display_frame(early, max_points)
                inset.plot(selected["time_s"] / 3600, C * selected["bias_s"], color=color)
                inset.set(xlim=(0, min(7200, horizon) / 3600), xlabel="Time [h]", title="First two hours")
                inset.tick_params(labelsize=8)
                inset.grid(True, alpha=0.25)
                if threshold_m is not None:
                    inset.axhline(threshold_m, color="0.35", linestyle="--")
                    inset.axhline(-threshold_m, color="0.35", linestyle="--")
                    contact = first_recorded_contact(early, C * early["bias_s"], threshold_m)
                    if contact is not None:
                        index, time, _ = contact
                        inset.scatter([time / 3600], [C * early["bias_s"].iloc[index]], marker="x", color=color)
        uncertainty_max = max([threshold_m or 0] + [float((k*C*frame["sigma_bias_s"]).max())
                              for frame in frames if "sigma_bias_s" in frame])
        for ax in sigma_axes[0]:
            ax.set_ylim(0, 1.12 * uncertainty_max if uncertainty_max else 1)
        figures[1].supxlabel("Deterministic aging excluded.\n" + assumptions, fontsize=8)
        figures[2].supxlabel(assumptions, fontsize=8)
        if threshold_m is not None:
            primary.supxlabel("Sampled crossings use full histories; excursions between samples may be missed. The chosen allocation is illustrative; no mandatory contacts.\n" + assumptions, fontsize=8)
        if allan_curves is not None:
            allan_fig, allan_axes = make_figure(1, "Overlapping Allan deviation — drift retained")
            for column, ((tau, deviation), label) in enumerate(zip(allan_curves, labels)):
                ax = allan_axes[0, column]
                ax.set_xscale("log")
                ax.set_yscale("log")
                positive = deviation > 0
                if positive.any():
                    ax.plot(tau[positive], deviation[positive], "o-", color=f"C{column % 10}", label=label)
                else:
                    ax.text(0.5, 0.5, "ADEV = 0; cannot display zero on a log axis", ha="center",
                            va="center", transform=ax.transAxes)
                if q_bias_s is not None:
                    theory = white_fm_adev(tau, q_bias_s[column], drift_per_s[column])
                    nonzero = theory > 0
                    ax.plot(tau[nonzero], theory[nonzero], "--", color="0.25", label="White FM + aging theory")
                if ax.lines:
                    ax.legend()
                ax.set(title=label, xlabel="Averaging time tau [s]",
                       ylabel="Allan deviation sigma_y(tau) [1]")
                ax.grid(True, which="both", alpha=0.25)
            positive_values = [adev[adev > 0] for _, adev in allan_curves if (adev > 0).any()]
            if positive_values:
                shared = np.concatenate(positive_values)
                for ax in allan_axes[0]:
                    ax.set_ylim(shared.min() / 1.3, shared.max() * 1.3)
            allan_fig.supxlabel("Full uniform bias histories; drift retained; tau ≤ record duration / 10; long-tau estimates remain uncertain.\n" + assumptions, fontsize=8)
        for fig, path in zip(figures, paths):
            fig.savefig(path, dpi=180)
        if show:
            plt.show() # All windows are opened together.
        return paths
    finally:
        for fig in figures:
            plt.close(fig)


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv:
        argv = setup_arguments()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv_files", nargs="+", type=Path)
    parser.add_argument("--labels", nargs="+")
    parser.add_argument("--mean-csv", nargs="+", type=Path,
                        help="Deterministic mean histories, in the same order and at matching timestamps")
    parser.add_argument("--allan-csv", nargs="+", type=Path,
                        help="Full uniform histories for Allan analysis; defaults to sibling *_allan.csv files")
    parser.add_argument("--output", type=Path, default=Path("clocks.png"))
    parser.add_argument("--time-unit", choices=("seconds", "hours", "days"), default="hours")
    parser.add_argument("--k", type=float, default=3.0, help="Pointwise sigma multiplier (default: 3)")
    parser.add_argument("--threshold-m", type=float)
    parser.add_argument("--max-points", type=int, default=10000,
                        help="Display limit per curve (default 10000); 0 keeps every sample")
    parser.add_argument("--q-b", nargs="+", type=float)
    parser.add_argument("--drift-per-s", nargs="+", type=float)
    parser.add_argument("--uso-inset-index", type=int)
    parser.add_argument("--show", action="store_true")
    args = parser.parse_args(argv)
    if args.labels is not None and len(args.labels) != len(args.csv_files):
        parser.error("--labels must provide one label per CSV file")
    if args.mean_csv is not None and len(args.mean_csv) != len(args.csv_files):
        parser.error("--mean-csv must provide one deterministic mean file per CSV file")
    if args.allan_csv is not None and len(args.allan_csv) != len(args.csv_files):
        parser.error("--allan-csv must provide one full history per CSV file")
    if args.max_points < 0 or args.max_points == 1:
        parser.error("--max-points must be 0 or >= 2 to preserve both endpoints")
    if not np.isfinite(args.k) or args.k < 0:
        parser.error("--k must be finite and nonnegative")
    if args.threshold_m is not None and (not np.isfinite(args.threshold_m) or args.threshold_m <= 0):
        parser.error("--threshold-m must be finite and positive")
    if args.output.suffix.lower() not in (".png", ".pdf", ".svg"):
        parser.error("--output must have a PNG, PDF or SVG extension")
    labels = args.labels or [path.stem for path in args.csv_files]
    try:
        frames = [load_clock_csv(path) for path in args.csv_files]
        means = [load_clock_csv(path) for path in args.mean_csv] if args.mean_csv is not None else None
        allan_paths = args.allan_csv
        if allan_paths is None:
            companions = [path.with_name(path.stem + "_allan.csv") for path in args.csv_files]
            if all(path.is_file() for path in companions):
                allan_paths = companions
            else:
                print("Allan plot skipped: provide full uniform histories with --allan-csv "
                      "or generate sibling *_allan.csv files.")
        allan_frames = [load_clock_csv(path) for path in allan_paths] if allan_paths is not None else None
        outputs = plot_histories(frames, labels, args.output, args.time_unit, args.k,
                                args.threshold_m, args.max_points, args.show, means, allan_frames,
                                args.q_b, args.drift_per_s, args.uso_inset_index)
    except (OSError, ValueError, ImportError) as error:
        print(f"Error: {error}", file=sys.stderr)
        return 1
    for output in outputs:
        print(f"Saved {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
