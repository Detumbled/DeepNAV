#!/usr/bin/env python3
"""Plot the deterministic 260-observation GMAT comparison."""
import argparse
import csv
import math
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('directory', type=Path)
args = parser.parse_args()
with (args.directory / 'comparison_260.csv').open() as file:
    rows = list(csv.DictReader(file))
if len(rows) != 260 or any(not math.isfinite(float(row['error_m_or_hz'])) for row in rows):
    raise ValueError('Expected 260 finite benchmark residuals')
fig, axes = plt.subplots(2, 1, figsize=(10, 6), sharex=True, constrained_layout=True)
colors = {'DSS-14': '#0072B2', 'DSS-43': '#D55E00', 'DSS-63': '#009E73'}
for ax, kind, label in zip(axes, ('Range', 'DSN_TCP'), ('Range residual [m]', 'Counted Doppler residual [Hz]')):
    for station, color in colors.items():
        subset = [r for r in rows if r['type'] == kind and r['station'] == station]
        ax.scatter([float(r['elapsed_s']) / 3600 for r in subset],
                   [float(r['error_m_or_hz']) for r in subset], s=18, label=station, color=color)
    ax.axhline(0, color='0.4', linewidth=0.7)
    ax.set_ylabel(label)
    ax.grid(alpha=0.25)
axes[0].legend(ncol=3)
axes[0].set_title('Cassini: DeepNAV − GMAT, 130 range + 130 Doppler observations\nShared trajectories; 1 µs transponder delay; corrections and noise disabled')
axes[1].set_xlabel('Hours from 2004-10-10 18:13:06.046 UTC')
axes[1].set_xlim(0, 24)
fig.savefig(args.directory / 'residuals_24h.png', dpi=180)
fig.savefig(args.directory / 'residuals_24h.pdf')
print(args.directory / 'residuals_24h.png')

# Compare stations on the same time and residual scales. Each point is an
# actual C++ TwoWayLink prediction minus the converted native GMAT observable.
limits = [ax.get_ylim() for ax in axes]
labels = ('Range residual [m]', 'Counted Doppler residual [Hz]')
kinds = ('Range', 'DSN_TCP')


def station_panels(station_axes, station):
    station_rows = [row for row in rows if row['station'] == station]
    counts = []
    for index, (ax, kind, label) in enumerate(zip(station_axes, kinds, labels)):
        subset = [row for row in station_rows if row['type'] == kind]
        counts.append(len(subset))
        ax.scatter([float(row['elapsed_s']) / 3600 for row in subset],
                   [float(row['error_m_or_hz']) for row in subset],
                   s=18, color=colors[station])
        ax.axhline(0, color='0.4', linewidth=0.7)
        ax.set_ylabel(label)
        ax.set_xlim(0, 24)
        ax.set_ylim(limits[index])
        ax.grid(alpha=0.25)
    if counts[0] != counts[1]:
        raise ValueError(f'Unpaired observations for {station}')
    station_axes[0].set_title(f'{station}: {counts[0]} ranges + {counts[1]} Doppler counts')
    station_axes[1].set_xlabel('Elapsed time [h]')


station_fig, station_axes = plt.subplots(2, 3, figsize=(15, 7), sharex=True,
                                       constrained_layout=True)
station_fig.suptitle('Cassini: DeepNAV C++ − GMAT, all 260 benchmark observations\n'
                    'Shared trajectories; 1 µs transponder delay; corrections and noise disabled')
for column, station in enumerate(colors):
    station_panels(station_axes[:, column], station)
    individual_fig, individual_axes = plt.subplots(2, 1, figsize=(10, 6), sharex=True,
                                                  constrained_layout=True)
    station_panels(individual_axes, station)
    individual_fig.suptitle('Cassini: DeepNAV C++ − GMAT\n'
                           '1 µs transponder delay; corrections and noise disabled')
    individual_axes[1].set_xlabel('Hours from 2004-10-10 18:13:06.046 UTC')
    for extension in ('png', 'pdf'):
        individual_fig.savefig(args.directory / f'residuals_24h_{station}.{extension}', dpi=180)
    plt.close(individual_fig)
for extension in ('png', 'pdf'):
    station_fig.savefig(args.directory / f'residuals_24h_by_station.{extension}', dpi=180)
print(args.directory / 'residuals_24h_by_station.png')
plt.close(station_fig)
plt.close(fig)
