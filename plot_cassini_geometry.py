"""Extract SPICE positions and display Cassini and Enceladus as 3D points."""

import argparse
import csv
import json
import os
from pathlib import Path
import webbrowser

import numpy as np
import plotly.graph_objects as go
import spiceypy as spice

ROOT = Path(__file__).resolve().parent
DEFAULT_EPOCH = "2004-10-10T18:13:06.046"
COLORS = {"SATURN": "#b88b38", "ENCELADUS": "#1686b0", "CASSINI": "#d55e00"}


def extract_positions(kernel, epoch):
    """Return simultaneous geometric states in Saturn-centered J2000, in km."""
    previous_directory = Path.cwd()
    try:
        # The example meta-kernel resolves paths relative to the repository root.
        os.chdir(ROOT)
        spice.furnsh(str(kernel))
        et = spice.str2et(epoch)
        states = {
            body: np.asarray(spice.spkezr(body, et, "J2000", "NONE", "SATURN")[0])
            for body in COLORS
        }
        relative = states["CASSINI"] - states["ENCELADUS"]
        _, light_time = spice.spkezr("ENCELADUS", et, "J2000", "CN", "CASSINI")
        return {
            "epoch_utc": spice.et2utc(et, "ISOC", 3),
            "epoch_tdb_seconds": et,
            "frame": "J2000",
            "origin": "SATURN",
            "aberration_correction": "NONE",
            "meta_kernel": str(kernel),
            "states_km_km_s": {body: state.tolist() for body, state in states.items()},
            "cassini_minus_enceladus_km": relative[:3].tolist(),
            "separation_km": float(np.linalg.norm(relative[:3])),
            "optical_reception_light_time_s": float(light_time),
        }
    finally:
        spice.kclear()
        os.chdir(previous_directory)


def make_figure(data):
    figure = go.Figure()
    for body, state in data["states_km_km_s"].items():
        x, y, z = state[:3]
        figure.add_trace(go.Scatter3d(
            x=[x], y=[y], z=[z], mode="markers+text", name=body.title(),
            text=[body.title()],
            textposition="bottom center" if body == "SATURN" else "top center",
            marker={"size": 8, "color": COLORS[body]},
            hovertemplate=(f"{body.title()}<br>X: %{{x:,.1f}} km<br>"
                           "Y: %{y:,.1f} km<br>Z: %{z:,.1f} km<extra></extra>"),
        ))
    figure.update_layout(
        title={"text": f"Cassini and Enceladus — {data['epoch_utc']} UTC<br>"
               f"<sup>Separation: {data['separation_km']:,.1f} km · "
               "Saturn-centered J2000 · marker sizes are symbolic</sup>"},
        scene={"xaxis_title": "X [km]", "yaxis_title": "Y [km]",
               "zaxis_title": "Z [km]", "aspectmode": "data",
               "camera": {"eye": {"x": 1.5, "y": 1.5, "z": 0.9}}},
        template="plotly_white", margin={"l": 0, "r": 0, "b": 0, "t": 90},
        legend={"x": 0.01, "y": 0.98}, height=780,
    )
    return figure


def save_preview(data, destination):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FuncFormatter

    fig = plt.figure(figsize=(10, 8), constrained_layout=True)
    ax = fig.add_subplot(projection="3d")
    points = np.array([state[:3] for state in data["states_km_km_s"].values()])
    for body, state in data["states_km_km_s"].items():
        ax.scatter(*state[:3], color=COLORS[body], s=65, depthshade=False,
                   label=body.title())
    # Equal axis ranges keep Euclidean distances visually comparable.
    middle = (points.min(axis=0) + points.max(axis=0)) / 2
    half_width = max(float(np.ptp(points, axis=0).max()) * 0.6, 1.0)
    for body, state in data["states_km_km_s"].items():
        offset = -0.15 if body == "SATURN" else 0.15
        ax.text(state[0], state[1], state[2] + offset * half_width,
                body.title(), fontsize=10)
    for index, axis in enumerate((ax.xaxis, ax.yaxis, ax.zaxis)):
        axis.set_major_formatter(FuncFormatter(lambda value, _: f"{value / 1e6:.1f}"))
        (ax.set_xlim, ax.set_ylim, ax.set_zlim)[index](
            middle[index] - half_width, middle[index] + half_width)
    ax.set_box_aspect((1, 1, 1))
    ax.set(xlabel="X [million km]", ylabel="Y [million km]", zlabel="Z [million km]")
    ax.view_init(elev=24, azim=40)
    ax.legend(loc="upper left")
    ax.set_title(f"Cassini and Enceladus — {data['epoch_utc']} UTC\n"
                 f"Separation: {data['separation_km']:,.0f} km", pad=20)
    fig.supxlabel("Saturn-centered J2000 · simultaneous geometric positions · symbolic point sizes",
                  fontsize=10)
    fig.savefig(destination, dpi=160)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--epoch", default=DEFAULT_EPOCH, help="UTC epoch")
    parser.add_argument("--kernels", type=Path,
                        default=ROOT / "Kernels/Cassini/cassini_2004_opnav.tm")
    parser.add_argument("--output", type=Path,
                        default=ROOT / "Output optical/cassini_geometry")
    parser.add_argument("--no-show", action="store_true", help="Do not open the HTML viewer")
    args = parser.parse_args()
    kernel, output = args.kernels.resolve(), args.output.resolve()
    if not kernel.is_file():
        parser.error(f"Missing meta-kernel: {kernel}")
    try:
        data = extract_positions(kernel, args.epoch)
    except spice.utils.exceptions.SpiceyError as error:
        parser.exit(1, f"SPICE geometry extraction failed: {error}\n")
    output.mkdir(parents=True, exist_ok=True)
    (output / "geometry.json").write_text(json.dumps(data, indent=2) + "\n")
    with (output / "positions.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["body", "x_km", "y_km", "z_km", "vx_km_s", "vy_km_s", "vz_km_s"])
        writer.writerows((body, *state) for body, state in data["states_km_km_s"].items())
    html = output / "cassini_enceladus_3d.html"
    make_figure(data).write_html(html, include_plotlyjs=True, config={"displaylogo": False})
    save_preview(data, output / "cassini_enceladus_3d.png")
    print(f"Cassini–Enceladus separation: {data['separation_km']:,.3f} km")
    print(f"Enceladus reception light-time: {data['optical_reception_light_time_s']:.6f} s")
    print(f"Saved viewer and extracted data in {output}")
    if not args.no_show:
        webbrowser.open(html.as_uri())


if __name__ == "__main__":
    main()
