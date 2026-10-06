"""Plot only running 3D position RMS and prefit NIS for the two-way EKF arc."""
import argparse
import json
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=Path("Output EKF/Cassini two-way"))
    args = parser.parse_args()
    data = np.genfromtxt(args.input / "diagnostics.csv", delimiter=",", names=True,
                         dtype=None, encoding="utf-8", ndmin=1)
    summary = json.loads((args.input / "summary.json").read_text())
    elapsed = data["reception_elapsed_s"] / 3600
    tracking = data["tracking"].astype(bool)
    rms = data["running_position_rms_m"]
    # Verify the statistic against the saved pointwise errors, including gaps.
    expected = np.sqrt(np.cumsum(data["position_error_m"] ** 2) / np.arange(1, len(data) + 1))
    if not np.allclose(rms, expected, rtol=1e-10):
        raise ValueError("Saved position RMS does not match the diagnostic errors")
    if summary["measurement_dimension"] != 2:
        raise ValueError("This NIS reference is for two measurement components")
    figure, axes = plt.subplots(2, 1, figsize=(10, 7), sharex=True, layout="constrained")
    axes[0].plot(elapsed, rms, color="#1864ab", linewidth=2)
    axes[0].set_ylabel("Running 3D position RMS [m]")
    axes[0].set_title("Cassini two-way EKF · 24-hour arc" if summary["arc_hours"] == 24
                      else f"Cassini two-way EKF · {summary['arc_hours']:g}-hour arc")
    axes[1].plot(elapsed[tracking], data["nis"][tracking], ".", color="#0b7285",
                 markersize=4, label="Prefit NIS")
    axes[1].axhline(2, color="#868e96", linewidth=1, label="Expected mean (2 components)")
    # chi-square(2) has CDF 1-exp(-x/2); this is a single-update 95th percentile.
    axes[1].axhline(-2 * np.log(.05), color="#e67700", linestyle="--", linewidth=1,
                    label="95th percentile (single update)")
    axes[1].set_ylabel("Normalized innovation squared")
    axes[1].set_xlabel("Earth reception time since arc start [h]")
    axes[1].legend(loc="upper right", frameon=False, fontsize=9)
    for axis in axes:
        axis.grid(alpha=.2)
        axis.spines[["top", "right"]].set_visible(False)
    output = args.input / "ekf_rms_nis.png"
    figure.savefig(output, dpi=180)
    plt.close(figure)
    print(output)


if __name__ == "__main__":
    main()
