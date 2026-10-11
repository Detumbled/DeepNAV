"""Plot saved C++ EKF prefit innovations, without using GMAT observations."""
import argparse
import csv
import json
import math
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

COLORS = {'DSS-14': '#0072B2', 'DSS-43': '#D55E00', 'DSS-63': '#009E73'}
COMPONENTS = [('range_innovation_m', 'Range innovation [m]'),
              ('doppler_innovation_hz', 'Counted Doppler innovation [Hz]')]


def panels(axes, rows, station, zoom):
    selected = [r for r in rows if r['station'] == station and r['tracking'] == '1']
    for index, (ax, (column, label)) in enumerate(zip(axes, COMPONENTS)):
        ax.scatter([float(r['reception_elapsed_s']) / 3600 for r in selected],
                   [float(r[column]) for r in selected], s=19, color=COLORS[station])
        ax.axhline(0, color='0.4', linewidth=0.7)
        ax.set_ylabel(label)
        ax.set_xlim(0, 24)
        ax.grid(alpha=0.25)
        if zoom:
            ax.set_ylim((-35, 35) if index == 0 else (-0.05, 0.05))
            clipped = [r for r in selected if abs(float(r[column])) > (35 if index == 0 else 0.05)]
            if clipped:
                details = ', '.join(f"t={float(r['reception_elapsed_s'])/3600:g} h: {float(r[column]):+.3f}" for r in clipped)
                ax.text(0.03, 0.96, 'Outside zoom: ' + details, transform=ax.transAxes,
                        va='top', fontsize=9, bbox={'facecolor': 'white', 'alpha': 0.9, 'edgecolor': '0.8'})
    axes[0].set_title(f'{station}: {len(selected)} range/Doppler pairs')
    axes[1].set_xlabel('Earth reception time since arc start [h]')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, default=Path('Output EKF/Cassini two-way'))
    args = parser.parse_args()
    with (args.input / 'diagnostics.csv').open() as file:
        rows = list(csv.DictReader(file))
    metadata = json.loads((args.input / 'summary.json').read_text())
    tracking = [r for r in rows if r['tracking'] == '1']
    if len(tracking) != metadata['updates']:
        raise ValueError('Tracking count disagrees with the C++ summary')
    for row in tracking:
        if any(not math.isfinite(float(row[column])) for column, _ in COMPONENTS):
            raise ValueError('Non-finite innovation at a tracking epoch')
    # Full versions preserve every innovation, including the initial transient.
    # Zoom versions use common linear axes and explicitly label clipped points.
    for zoom in (False, True):
        suffix = '_zoom' if zoom else ''
        figure, axes = plt.subplots(2, 3, figsize=(15, 7), constrained_layout=True)
        title = f"Cassini C++ EKF: {len(tracking)} range + {len(tracking)} Doppler prefit innovations"
        figure.suptitle(title + '\nSynthetic observation − prediction; no GMAT' + ('; common-axis zoom' if zoom else '; full data range'))
        for col, station in enumerate(COLORS):
            panels(axes[:, col], rows, station, zoom)
            individual, individual_axes = plt.subplots(2, 1, figsize=(10, 6), constrained_layout=True)
            panels(individual_axes, rows, station, zoom)
            individual.suptitle('Cassini C++ EKF: synthetic observation − prediction\n'
                                + ('Common-axis zoom; initial clipped points labelled' if zoom else 'All prefit innovations; no GMAT'))
            for extension in ('png', 'pdf'):
                individual.savefig(args.input / f'cpp_innovations_{station}{suffix}.{extension}', dpi=180)
            plt.close(individual)
        for extension in ('png', 'pdf'):
            figure.savefig(args.input / f'cpp_innovations_by_station{suffix}.{extension}', dpi=180)
        plt.close(figure)
    print(args.input / 'cpp_innovations_by_station.png')
    print(args.input / 'cpp_innovations_by_station_zoom.png')

    # This independently reconstructed SPICE visibility table verifies the saved
    # C++ station choices. It is optional and never supplies innovation values.
    visibility = args.input / 'visibility_check.csv'
    if visibility.exists():
        with visibility.open() as file:
            geometry = list(csv.DictReader(file))
        figure, ax = plt.subplots(figsize=(11, 5), constrained_layout=True)
        for station, color in COLORS.items():
            column = station + '_min_elevation_deg'
            ax.plot([float(r['elapsed_s'])/3600 for r in geometry],
                    [float(r[column]) for r in geometry], color=color, label=station, linewidth=1.5)
            selected = [r for r in geometry if r['selected_station'] == station and r['tracking'] == '1']
            ax.scatter([float(r['elapsed_s'])/3600 for r in selected],
                       [float(r[column]) for r in selected], color=color, s=14)
        ax.axhline(10, color='0.3', linestyle='--', label='10° tracking limit')
        ax.set(xlim=(0, 24), ylim=(-90, 90), xlabel='Earth reception time since arc start [h]',
               ylabel='Minimum elevation over both legs and count endpoints [deg]',
               title='C++ tracking schedule: dots mark the selected station\n'
                     'DSS-14 unselected at sampled epochs from 1.83 h until 17.83 h')
        ax.legend(ncol=4)
        ax.grid(alpha=0.25)
        for extension in ('png', 'pdf'):
            figure.savefig(args.input / f'cpp_station_visibility.{extension}', dpi=180)
        plt.close(figure)
        print(args.input / 'cpp_station_visibility.png')


if __name__ == '__main__':
    main()
