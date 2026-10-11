#!/usr/bin/env python3
"""Prepare native GMAT scripts and compare their observations with DeepNAV."""
import argparse
import csv
import hashlib
import json
import statistics
import math
from pathlib import Path
import re
import subprocess

ROOT = Path(__file__).resolve().parents[2]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--gmat', type=Path, default=Path('/Applications/GMAT R2026a'))
    parser.add_argument('--build', type=Path, default=ROOT / 'build-clang')
    parser.add_argument('--prepare-only', action='store_true', help='Generate scripts without running GMAT or CMake')
    parser.add_argument('--output', type=Path, help='Output directory (default: BUILD/gmat-benchmark)')
    parser.add_argument('--full-day', action='store_true', help='Compare the 130 visible updates (260 scalar observations) in the saved Cassini EKF schedule')
    parser.add_argument('--schedule', type=Path, default=ROOT / 'Output EKF/Cassini two-way/diagnostics.csv')
    parser.add_argument('--state-step', type=float, default=1., help='Trajectory sampling step in seconds; use 0.5 for a refinement check')
    args = parser.parse_args()
    install = args.gmat.resolve()
    build = args.build.resolve()
    if not math.isfinite(args.state_step) or args.state_step <= 0 or args.state_step > 10:
        parser.error('--state-step must be finite and in (0, 10] seconds')
    out = args.output.resolve() if args.output else build / ('gmat-benchmark-24h' if args.full_day else 'gmat-benchmark')
    out.mkdir(parents=True, exist_ok=True)
    # Preserve GMAT's shipped paths, but write logs/output inside this project.
    # Only load the plugins needed for SPK propagation and navigation.
    shipped = install / 'bin/gmat_startup_file.txt'
    selected = {'libEphemPropagator', 'libStation', 'libGmatEstimation', 'libProductionPropagators'}
    lines = []
    for line in shipped.read_text().splitlines():
        if line.startswith('PLUGIN') and not any(plugin in line for plugin in selected):
            continue
        line = line.replace('../', str(install) + '/')
        if re.match(r'OUTPUT_PATH\s*=', line):
            line = f'OUTPUT_PATH = {out}/'
        elif re.match(r'LSK_FILE\s*=', line):
            line = f'LSK_FILE = {ROOT}/Kernels/naif0012.tls'
        lines.append(line)
    startup = out / 'startup.txt'
    startup.write_text('\n'.join(lines) + '\n')
    template = (ROOT / 'benchmarks/gmat/cassini_two_way.script.in').read_text()
    cases = [('zero_delay', '0'), ('delay_1us', '1e-6'), ('delay_1ms', '1e-3')]
    schedule = None
    if args.full_day:
        with args.schedule.open() as file:
            all_rows = list(csv.DictReader(file))
        schedule = {int(float(row['reception_elapsed_s'])): row['station'] for row in all_rows if row['tracking'] == '1'}
        if len(all_rows) != 145 or len(schedule) != 130 or set(schedule.values()) != {'DSS-14', 'DSS-43', 'DSS-63'}:
            raise ValueError('Expected the original 24-hour, 145-epoch, 130-update Cassini schedule')
        (out / 'schedule.csv').write_text(args.schedule.read_text())
        cases = [('DSS-14', '1e-6'), ('DSS-43', '1e-6'), ('DSS-63', '1e-6')]
        template = template.replace("Sim.FinalEpoch = '10 Oct 2004 18:33:06.046';", "Sim.FinalEpoch = '11 Oct 2004 18:13:06.046';")
    required = [ROOT / 'Kernels/de442.bsp', ROOT / 'Kernels/naif0012.tls',
                ROOT / 'Kernels/Cassini/spk/200128RU_SCPSE_04251_04303.bsp']
    for file in required:
        if not file.is_file():
            raise FileNotFoundError(file)
    for name, delay in cases:
        case = out / name
        case.mkdir(exist_ok=True)
        script = (template.replace('@ROOT@', str(ROOT)).replace('@OUTPUT@', str(case))
                  .replace('@DELAY@', delay).replace('@STEP@', str(args.state_step))
                  .replace('@SAMPLES@', str(math.ceil((98000 if args.full_day else 12000) / args.state_step))))
        if args.full_day:
            locations = {'DSS-14': (-2353.621251, -4641.341542, 3677.052370),
                         'DSS-43': (-4460.894917, 2682.361507, -3674.748152),
                         'DSS-63': (4849.092518, -360.180348, 4115.109251)}
            for axis, value in enumerate(locations[name], 1):
                script = re.sub(rf'DSN.Location{axis} = [^;]+;', f'DSN.Location{axis} = {value};', script)
            script = script.replace("DSN.Id = '14';", f"DSN.Id = '{name[-2:]}';")
        (case / 'cassini_two_way.script').write_text(script)
    print(f'Generated standalone GMAT scripts in {out}', flush=True)
    if args.prepare_only:
        return
    subprocess.run(['cmake', '--build', str(build), '--target', 'compare_gmat_two_way', '-j', '4'], check=True)
    for name, delay in cases:
        case = out / name
        command = [str(install / 'bin/GmatConsole'), '-s', str(startup), '-r', str(case / 'cassini_two_way.script')]
        with (case / 'console.log').open('w') as log:
            subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True, cwd=out)
        if 'Mission run completed.' not in (case / 'console.log').read_text():
            raise RuntimeError(f'GMAT did not complete {name}; inspect console.log')
        with (case / 'comparison.csv').open('w') as comparison_file, (case / 'comparison.log').open('w') as log:
            subprocess.run([str(build / 'compare_gmat_two_way'), str(case / 'states.csv'),
                            str(case / 'measurements.gmd'), delay, str(required[1])],
                           stdout=comparison_file, stderr=log, check=True)
        print(f'{name}: {(case / "comparison.log").read_text().strip()}', flush=True)
    if args.full_day:
        summarize_full_day(out, cases, schedule, args.schedule, args.state_step)
        return
    ranges = {}
    for name, _ in cases:
        with (out / name / 'comparison.csv').open() as file:
            ranges[name] = [row for row in csv.DictReader(file) if row['type'] == 'Range']
        if len(ranges[name]) != 3:
            raise ValueError(f'Expected three range/Doppler pairs in short case {name}')
    with (out / 'delay_comparison.csv').open('w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['case', 'tai_mjd', 'gmat_range_increase_m', 'deepnav_range_increase_m', 'difference_m'])
        for name, _ in cases[1:]:
            for base, delayed in zip(ranges['zero_delay'], ranges[name]):
                if base['tai_mjd'] != delayed['tai_mjd']:
                    raise RuntimeError('Delay cases have different reception epochs')
                gmat = (float(delayed['gmat_converted']) - float(base['gmat_converted'])) * 1000
                deepnav = (float(delayed['deepnav']) - float(base['deepnav'])) * 1000
                writer.writerow([name, base['tai_mjd'], gmat, deepnav, deepnav - gmat])
    print('Residuals are diagnostic; this runner does not declare a tolerance-based validation pass.')


def summarize_full_day(out, cases, schedule, source, state_step):
    # GMD uses TAI tags; UTC + 32 seconds at the fixed 2004 arc start.
    start_tai_mjd = 23289.259468125
    selected = {}
    for station, _ in cases:
        with (out / station / 'comparison.csv').open() as file:
            for row in csv.DictReader(file):
                elapsed = (float(row['tai_mjd']) - start_tai_mjd) * 86400
                grid = round(elapsed / 600) * 600
                if abs(elapsed - grid) > 0.001:
                    raise ValueError('GMAT observation is outside the requested epoch grid')
                if schedule.get(grid) != station:
                    continue
                key = (grid, row['type'])
                if key in selected:
                    raise ValueError(f'Duplicate observation: {key}')
                row['elapsed_s'] = grid
                row['epoch_offset_s'] = float(row['elapsed_tdb_s']) - grid
                if abs(row['epoch_offset_s']) > 0.001:
                    raise ValueError('GMAT TDB tag differs from the saved schedule by more than 1 ms')
                row['station'] = station
                selected[key] = row
    expected = {(epoch, kind) for epoch in schedule for kind in ('Range', 'DSN_TCP')}
    if set(selected) != expected:
        missing = sorted(expected - set(selected))
        raise ValueError(f'Missing requested observations: {missing}; GMAT visibility differs from saved schedule')
    rows = [selected[key] for key in sorted(selected)]
    with (out / 'comparison_260.csv').open('w', newline='') as file:
        writer = csv.DictWriter(file, fieldnames=['elapsed_s', 'station', 'tai_mjd', 'elapsed_tdb_s', 'epoch_offset_s', 'type', 'gmat_converted', 'deepnav', 'error_m_or_hz'])
        writer.writeheader()
        writer.writerows(rows)
    summary = {'arc_hours': 24, 'scheduled_epochs': 145, 'visible_updates': 130,
               'scalar_observations_compared': len(rows), 'transponder_delay_s': 1e-6,
               'state_step_s': state_step, 'schedule_source': str(source.resolve()),
               'schedule_sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
               'corrections_enabled': False, 'noise_enabled': False,
               'epoch_matching_tolerance_s': 0.001,
               'max_epoch_offset_s': max(abs(row['epoch_offset_s']) for row in rows),
               'statistics': {}}
    for station in ['all', 'DSS-14', 'DSS-43', 'DSS-63']:
        summary['statistics'][station] = {}
        for kind in ('Range', 'DSN_TCP'):
            errors = [float(row['error_m_or_hz']) for row in rows if row['type'] == kind and (station == 'all' or row['station'] == station)]
            summary['statistics'][station][kind] = {'count': len(errors), 'unit': 'm' if kind == 'Range' else 'Hz',
                'mean': statistics.fmean(errors), 'rms': math.sqrt(statistics.fmean(x*x for x in errors)),
                'max_abs': max(abs(x) for x in errors)}
    (out / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print('Compared all 260 requested scalar observations. Statistics:', flush=True)
    print(json.dumps(summary['statistics']['all'], indent=2), flush=True)
    print('Corrections and noise disabled; residuals are diagnostic, not a tolerance-based validation pass.')


if __name__ == '__main__':
    main()
