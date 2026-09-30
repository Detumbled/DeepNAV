#!/usr/bin/env python3
"""Plot clock history CSVs. Requires Pandas and Matplotlib; no repository imports.

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


def plot_histories(frames, labels, output, time_unit="hours", k=3.0,
                   threshold_m=None, max_points=20_000, show=False, mean_frames=None):
    import matplotlib
    if not show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    if not frames or len(frames) != len(labels):
        raise ValueError("provide one label per clock history")
    if mean_frames is not None:
        if len(mean_frames) != len(frames):
            raise ValueError("provide one deterministic mean CSV per clock history")
        for frame, mean, label in zip(frames, mean_frames, labels):
            if "sigma_bias_s" not in frame:
                raise ValueError(f"{label}: total budget requires covariance")
            if not np.array_equal(frame["time_s"].to_numpy(), mean["time_s"].to_numpy()):
                raise ValueError(f"{label}: deterministic mean timestamps must match history exactly")

    scale, unit = {"seconds": (1.0, "s"), "hours": (3600.0, "h"),
                   "days": (86400.0, "days")}[time_unit]
    output = Path(output)
    paths = [output,
             output.with_name(output.stem + "_uncertainty" + output.suffix),
             output.with_name(output.stem + "_diagnostics" + output.suffix)]
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

    def draw(ax, frame, values, label, ylabel, color, signed_threshold=False):
        shown = display_frame(frame, max_points)
        ax.plot(shown["time_s"] / scale, values.loc[shown.index], color=color, label=label)
        ax.set(xlabel=f"Elapsed time since calibration [{unit}]", ylabel=ylabel)
        ax.grid(True, alpha=0.25)
        if threshold_m is not None:
            ax.axhline(threshold_m, color="0.35", linestyle="--", label="Range threshold")
            if signed_threshold:
                ax.axhline(-threshold_m, color="0.35", linestyle="--")
            contact = first_recorded_contact(frame, values, threshold_m)
            if contact is None:
                status = f"No recorded contact through {frame['time_s'].iloc[-1] / scale:.4g} {unit}"
            else:
                index, time, previous = contact
                ax.scatter([time / scale], [values.iloc[index]], color=color, marker="x", zorder=5)
                status = f"First recorded contact: {time / scale:.4g} {unit}"
                if previous is not None:
                    status += f"\nPrevious sample: {previous / scale:.4g} {unit}"
            ax.text(0.02, 0.97, status, transform=ax.transAxes, va="top", fontsize=9,
                    bbox={"facecolor": "white", "alpha": 0.85, "edgecolor": "0.8"})
        ax.legend(loc="lower right", fontsize=9)

    try:
        primary, main_axes = make_figure(2, "One-way range error and total pointwise clock budget")
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
                     f"Stochastic range uncertainty ({k:g} sigma) [m]", color)
            else:
                sigma_axes[0, column].text(0.5, 0.5, "No covariance supplied", ha="center", va="center",
                                          transform=sigma_axes[0, column].transAxes)
            sigma_axes[0, column].set_title(label + " — excludes deterministic bias")
            # Diagnostics deliberately have no range thresholds.
            shown = display_frame(frame, max_points)
            for row, (values, ylabel) in enumerate(((bias_ns, "Clock bias [ns]"),
                                                    (frame["fractional_frequency"], "Fractional frequency [1]"))):
                ax = diagnostic_axes[row, column]
                ax.plot(shown["time_s"] / scale, values.loc[shown.index], color=color, label=label)
                ax.set(xlabel=f"Elapsed time since calibration [{unit}]", ylabel=ylabel)
                ax.grid(True, alpha=0.25)
                ax.legend(loc="best")
            diagnostic_axes[0, column].set_title(label + " — independent scale")
        if threshold_m is not None:
            primary.supxlabel("Contact annotations use all recorded samples; crossings between samples may be missed.", fontsize=9)
        for fig, path in zip(figures, paths):
            fig.savefig(path, dpi=180)
        if show:
            plt.show() # All three windows are opened together.
        return paths
    finally:
        for fig in figures:
            plt.close(fig)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv_files", nargs="+", type=Path)
    parser.add_argument("--labels", nargs="+")
    parser.add_argument("--mean-csv", nargs="+", type=Path,
                        help="Deterministic mean histories, in the same order and at matching timestamps")
    parser.add_argument("--output", type=Path, default=Path("clocks.png"))
    parser.add_argument("--time-unit", choices=("seconds", "hours", "days"), default="hours")
    parser.add_argument("--k", type=float, default=3.0, help="Pointwise sigma multiplier (default: 3)")
    parser.add_argument("--threshold-m", type=float)
    parser.add_argument("--max-points", type=int, default=20_000, help="Display limit per curve; 0 keeps all")
    parser.add_argument("--show", action="store_true")
    args = parser.parse_args(argv)
    if args.labels is not None and len(args.labels) != len(args.csv_files):
        parser.error("--labels must provide one label per CSV file")
    if args.mean_csv is not None and len(args.mean_csv) != len(args.csv_files):
        parser.error("--mean-csv must provide one deterministic mean file per CSV file")
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
        outputs = plot_histories(frames, labels, args.output, args.time_unit, args.k,
                                args.threshold_m, args.max_points, args.show, means)
    except (OSError, ValueError, ImportError) as error:
        print(f"Error: {error}", file=sys.stderr)
        return 1
    for output in outputs:
        print(f"Saved {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
