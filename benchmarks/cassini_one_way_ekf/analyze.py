"""Validate and plot C++ one-way onboard EKF exports; never estimate states."""
from pathlib import Path
import argparse
import hashlib
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import chi2

ROOT = Path(__file__).resolve().parents[2]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output', type=Path, default=ROOT/'results/cassini_one_way_ekf')
args = parser.parse_args()
output = args.output.resolve()
colors = {'ideal': '#c46b35', 'dsac': '#187c9d'}
labels = {'ideal': 'Known ideal clock', 'dsac': 'Estimated DSAC-inspired clock'}
results = {}

def save(fig, folder, name):
    fig.savefig(folder/f'{name}.png', dpi=180)
    fig.savefig(folder/f'{name}.pdf')
    plt.close(fig)

def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()

for kind in colors:
    folder = output/kind
    if not (folder/'summary.json').exists():
        continue
    summary = json.loads((folder/'summary.json').read_text())
    state = pd.read_csv(folder/'trajectory.csv')
    innovations = pd.read_csv(folder/'innovations.csv')
    observations = pd.read_csv(folder/'observations.csv', float_precision='round_trip')
    covariance = pd.read_csv(folder/'covariance.csv')
    assert len(state) == len(covariance) == summary['diagnostic_samples']
    assert len(innovations) == len(observations) == summary['doppler_observations']
    assert int(observations.range_active.sum()) == summary['range_observations']
    assert state.elapsed_s.is_monotonic_increasing and state.elapsed_s.is_unique
    assert innovations.elapsed_s.is_monotonic_increasing and innovations.elapsed_s.is_unique
    assert (observations.minimum_elevation_deg >= 10).all()
    assert innovations.accepted.eq(1).all()
    assert np.isfinite(state.to_numpy()).all()
    np.testing.assert_array_equal(state.elapsed_s, covariance.elapsed_s)
    p = covariance.iloc[:, 1:].to_numpy().reshape(-1, 8, 8)
    scales = np.array([2, 2, 2, 1e-4, 1e-4, 1e-4, 1e-5, 1e-10])
    normalized = p/scales[None, :, None]/scales[None, None, :]
    np.testing.assert_allclose(normalized, normalized.transpose(0, 2, 1), atol=1e-10, rtol=1e-9)
    assert np.linalg.eigvalsh(normalized).min() >= -1e-10
    np.testing.assert_allclose(summary['position_rms_m'], np.sqrt(np.mean(state.position_error_m**2)), rtol=1e-12)
    np.testing.assert_allclose(summary['velocity_rms_mm_s'], np.sqrt(np.mean(state.velocity_error_mm_s**2)), rtol=1e-12)
    np.testing.assert_allclose(summary['mean_nis'], innovations.nis.mean(), rtol=1e-12)
    np.testing.assert_allclose(np.sqrt(np.trace(p[:, :3, :3], axis1=1, axis2=2))*1000, state.formal_position_sigma_m, rtol=1e-12)
    np.testing.assert_allclose(np.sqrt(np.trace(p[:, 3:6, 3:6], axis1=1, axis2=2))*1e6, state.formal_velocity_sigma_mm_s, rtol=1e-12)
    np.testing.assert_allclose(observations.observed_range_km-observations.truth_range_km, observations.range_noise_m/1000, rtol=0, atol=4e-7)
    np.testing.assert_allclose(observations.observed_doppler_hz-observations.truth_doppler_hz, observations.doppler_noise_hz, rtol=0, atol=2e-9)
    assert summary['maximum_half_stencil_noise_fraction'] < .01
    # In ideal mode the clock is known exactly rather than artificially regularized.
    if kind == 'ideal':
        np.testing.assert_array_equal(p[:, 6:, :], 0)
        np.testing.assert_allclose(state.clock_bias_error_s, 0, atol=0)
        np.testing.assert_allclose(state.clock_frequency_error, 0, atol=0)
    else:
        np.testing.assert_allclose(summary['clock_q_bias_s'], 2.25e-26, rtol=1e-15, atol=0)
        assert summary['clock_q_frequency_per_s'] == 0
    days = state.elapsed_s/86400
    analysis = pd.DataFrame({'elapsed_s': state.elapsed_s})
    fig, axes = plt.subplots(2, 1, figsize=(11, 7), sharex=True, layout='constrained')
    for ax, error, sigma, unit, title in zip(axes, ['position_error_m', 'velocity_error_mm_s'],
            ['formal_position_sigma_m', 'formal_velocity_sigma_mm_s'], ['m', 'mm/s'], ['Position', 'Velocity']):
        rms = np.sqrt(np.cumsum(state[error].to_numpy()**2)/np.arange(1, len(state)+1))
        bound = 3*np.sqrt(np.cumsum(state[sigma].to_numpy()**2)/np.arange(1, len(state)+1))
        analysis[f'{title.lower()}_cumulative_rms'] = rms
        analysis[f'{title.lower()}_minus_3sigma'] = -bound
        analysis[f'{title.lower()}_plus_3sigma'] = bound
        ax.fill_between(days, -bound, bound, color=colors[kind], alpha=.16, label='±3σ covariance reference')
        ax.plot(days, rms, color=colors[kind], label='Cumulative 3D RMS error')
        ax.set(ylabel=f'{title} RMS ({unit})', title=title)
        ax.grid(alpha=.2)
        ax.legend(fontsize=9)
    axes[-1].set_xlabel('Elapsed TDB days')
    fig.suptitle(f'Cassini one-way EKF · {labels[kind]}')
    save(fig, folder, 'rms_3sigma')
    analysis.to_csv(folder/'rms_analysis.csv', index=False)

    fig, axes = plt.subplots(2, 1, figsize=(11, 7), sharex=True, layout='constrained')
    for ax, prefix, suffix, unit in zip(axes, ['', '_velocity',], ['_error_m', '_error_mm_s'], ['m', 'mm/s']):
        for component, color in zip(['radial', 'along', 'cross'], ['#2878b5', '#dd8452', '#8a63ad']):
            column = f'{component}{prefix}{suffix}'
            ax.plot(days, state[column], label=component.capitalize(), color=color, lw=1)
        ax.set(ylabel=f'RTN error ({unit})')
        ax.legend()
        ax.grid(alpha=.2)
    axes[-1].set_xlabel('Elapsed TDB days')
    fig.suptitle(f'Forward EKF error against simulated truth · {labels[kind]}')
    save(fig, folder, 'orbit_errors')

    fig, axes = plt.subplots(2, 1, figsize=(11, 7), sharex=True, layout='constrained')
    for ax, value, sigma, factor, label in zip(axes,
            ['clock_bias_error_s', 'clock_frequency_error'], ['clock_bias_sigma_s', 'clock_frequency_sigma'],
            [1e6, 1e12], ['Clock bias error (µs)', 'Fractional frequency error (×10⁻¹²)']):
        ax.plot(days, state[value]*factor, color=colors[kind], label='Estimate − truth')
        ax.fill_between(days, -3*state[sigma]*factor, 3*state[sigma]*factor, color=colors[kind], alpha=.16, label='±3σ')
        ax.set_ylabel(label)
        ax.grid(alpha=.2)
        ax.legend()
    axes[-1].set_xlabel('Elapsed TDB days')
    fig.suptitle(f'Clock estimation · {labels[kind]}')
    save(fig, folder, 'clock_errors')

    fig, axes = plt.subplots(2, 1, figsize=(11, 7), sharex=True, layout='constrained')
    frequency_to_rate = -299792.458/7.175e9*1e6
    for station, color in zip(['DSS-14', 'DSS-43', 'DSS-63'], ['#2878b5', '#dd8452', '#8a63ad']):
        data = innovations[innovations.station == station]
        axes[0].scatter(data.elapsed_s/86400, data.range_innovation_m, s=7, label=station, color=color)
        axes[1].scatter(data.elapsed_s/86400, data.doppler_innovation_hz*frequency_to_rate, s=4, label=station, color=color)
    axes[0].set_ylabel('Range innovation (m)')
    axes[1].set(ylabel='Equivalent counted range-rate innovation (mm/s)', xlabel='Elapsed TDB days')
    for ax in axes:
        ax.grid(alpha=.2)
        ax.legend(fontsize=8)
    fig.suptitle(f'Unshaded pre-update innovations · {labels[kind]}')
    save(fig, folder, 'innovations')
    # Preserve the full view above; add an explicitly limited noise-scale view.
    fig, axes = plt.subplots(2, 1, figsize=(11, 7), sharex=True, layout='constrained')
    for station, color in zip(['DSS-14', 'DSS-43', 'DSS-63'], ['#2878b5', '#dd8452', '#8a63ad']):
        data = innovations[innovations.station == station]
        axes[0].scatter(data.elapsed_s/86400, data.range_innovation_m, s=7, label=station, color=color)
        axes[1].scatter(data.elapsed_s/86400, data.doppler_innovation_hz*frequency_to_rate, s=4, label=station, color=color)
    limits = [6*summary['range_sigma_m'], 6*summary['doppler_sigma_hz']*abs(frequency_to_rate)]
    values = [innovations.range_innovation_m, innovations.doppler_innovation_hz*frequency_to_rate]
    for ax, limit, value in zip(axes, limits, values):
        ax.set_ylim(-limit, limit)
        ax.text(.01, .97, f'{int((value.abs()>limit).sum())} initial/other point(s) outside view; retained in full plot and CSV', transform=ax.transAxes, va='top', fontsize=9)
        ax.grid(alpha=.2)
        ax.legend(fontsize=8)
    axes[0].set_ylabel('Range innovation (m)')
    axes[1].set(ylabel='Equivalent counted range-rate innovation (mm/s)', xlabel='Elapsed TDB days')
    fig.suptitle(f'Innovation detail · view limited to ±6× instrument noise · {kind}')
    save(fig, folder, 'innovations_detail')

    fig, axes = plt.subplots(2, 1, figsize=(11, 7), sharex=True, layout='constrained')
    axes[0].plot(innovations.elapsed_s/86400, innovations.nis/innovations.dimension, '.', ms=2, label='NIS / measurement dimension')
    axes[0].axhline(1, color='black', ls='--', label='Expected mean')
    for dimension in (1, 2):
        axes[0].axhline(chi2.ppf(.99, dimension)/dimension, ls=':', label=f'99% point threshold, {dimension} observable(s)')
    axes[0].set(ylabel='Normalized innovation statistic', yscale='log')
    nees_dimension = 6 if kind == 'ideal' else 8
    axes[1].plot(days, state.joint_nees/nees_dimension, color=colors[kind], label=f'NEES / {nees_dimension} active states')
    axes[1].axhline(1, color='black', ls='--', label='Expected ensemble mean')
    axes[1].set(ylabel='Normalized state-error statistic', xlabel='Elapsed TDB days', yscale='log')
    for ax in axes:
        ax.grid(alpha=.2)
        ax.legend(fontsize=8)
    fig.suptitle('Single-run consistency diagnostics · samples are temporally correlated')
    save(fig, folder, 'consistency')

    source_metadata = json.loads((ROOT/'results/cassini_two_way_wls/pass_schedule/truth/metadata.json').read_text())
    assert source_metadata['mass_kg'] == 3220 and source_metadata['area_m2'] == 12 and source_metadata['cr'] == 1.42 and source_metadata['gravity_mask'] == 2
    provenance_inputs = [ROOT/'results/cassini_two_way_wls/pass_schedule/truth'/name for name in ('trajectory.csv', 'metadata.json', 'summary.json')]
    provenance_inputs += [ROOT/'results/cassini_two_way_wls/pass_schedule/with_vlbi/passes.csv']
    sources = list((ROOT/'benchmarks/cassini_one_way_ekf').glob('*.cpp')) + list((ROOT/'benchmarks/cassini_one_way_ekf').glob('*.hpp'))
    sources += [ROOT/name for name in ('benchmarks/cassini_trajectory/SavedSimulation.hpp', 'benchmarks/cassini_trajectory/CassiniForceModel.hpp', 'benchmarks/cassini_trajectory/CassiniAttitude.hpp', 'benchmarks/cassini_trajectory/Sat427Gravity.hpp', 'src/filters/EKF.cpp', 'src/filters/CartesianPropagator.cpp', 'src/observations/radiometric/OneWayLink.cpp', 'src/observations/radiometric/OneWayMeasurement.cpp', 'src/observations/radiometric/EndpointClock.cpp', 'src/Clocks/DSAC.cpp', 'src/Clocks/ClockModel.cpp', 'src/Clocks/ClockTruthSimulator.cpp', 'src/RKF45Integrator.cpp')]
    validation = {'input_sha256': {str(p.relative_to(ROOT)): digest(p) for p in provenance_inputs},
        'source_sha256': {str(p.relative_to(ROOT)): digest(p) for p in sources},
        'output_sha256': {name: digest(folder/name) for name in ('observations.csv', 'trajectory.csv', 'innovations.csv', 'covariance.csv', 'initial_prior.csv', 'clock_truth.csv', 'summary.json')},
        'mean_nis_per_dimension': float((innovations.nis/innovations.dimension).mean()),
        'fraction_nis_above_pointwise_99_percent': float((innovations.nis > chi2.ppf(.99, innovations.dimension)).mean()),
        'mean_nees_per_active_state': float((state.joint_nees/nees_dimension).mean())}
    last = state.iloc[-1]
    r = last[['truth_x_km', 'truth_y_km', 'truth_z_km']].to_numpy(float)
    v = last[['truth_vx_km_s', 'truth_vy_km_s', 'truth_vz_km_s']].to_numpy(float)
    radial = r/np.linalg.norm(r)
    cross = np.cross(r, v); cross /= np.linalg.norm(cross)
    directions = {'radial': radial, 'along': np.cross(cross, radial), 'cross': cross}
    validation['endpoint_orbit_clock_correlations'] = {}
    if kind != 'ideal':
        for component, direction in directions.items():
            for clock_name, index in [('bias', 6), ('frequency', 7)]:
                denominator = np.sqrt((direction@p[-1, :3, :3]@direction)*p[-1, index, index])
                validation['endpoint_orbit_clock_correlations'][f'{component}_{clock_name}'] = float((direction@p[-1, :3, index])/denominator)
    (folder/'validation.json').write_text(json.dumps(validation, indent=2)+'\n')
    results[kind] = (summary, state, innovations, validation)
    print(kind, json.dumps({key: summary[key] for key in ('position_rms_m', 'endpoint_position_error_m', 'velocity_rms_mm_s', 'maximum_half_stencil_noise_fraction')}))

if len(results) == 2:
    ideal_obs = pd.read_csv(output/'ideal/observations.csv')
    dsac_obs = pd.read_csv(output/'dsac/observations.csv')
    for column in ['truth_elapsed_s', 'station', 'pass_id', 'range_active', 'range_noise_m', 'doppler_noise_hz']:
        np.testing.assert_array_equal(ideal_obs[column], dsac_obs[column])

if not results:
    raise RuntimeError('No completed C++ cases found')
fig, axes = plt.subplots(2, 1, figsize=(11, 7), sharex=True, layout='constrained')
for kind, (summary, state, innovations, validation) in results.items():
    axes[0].plot(state.elapsed_s/86400, state.position_error_m, color=colors[kind], label=labels[kind])
    axes[1].plot(state.elapsed_s/86400, state.velocity_error_mm_s, color=colors[kind], label=labels[kind])
axes[0].set_ylabel('3D position error (m)')
axes[1].set(ylabel='3D velocity error (mm/s)', xlabel='Elapsed TDB days')
for ax in axes:
    ax.grid(alpha=.2)
    ax.legend()
fig.suptitle('One-way forward EKF · identical orbital truth, separate clock assumptions')
save(fig, output, 'clock_case_comparison')
rows = [{'case': kind, **{key: s[key] for key in ('position_rms_m', 'velocity_rms_mm_s', 'endpoint_position_error_m', 'endpoint_velocity_error_mm_s', 'doppler_observations', 'range_observations')}, **{key: v[key] for key in ('mean_nis_per_dimension', 'mean_nees_per_active_state')}} for kind, (s, _, _, v) in results.items()]
pd.DataFrame(rows).to_csv(output/'comparison.csv', index=False)
lines = ['# Cassini one-way onboard EKF', '', 'C++ generates the one-way uplink observations and performs the causal EKF. Python validates and plots exports. Orbital truth is the unchanged retained simulation; the filter gets raw onboard tags and assigns coordinate epochs using its current estimated clock. No spacecraft SPK states, WLS estimates or VLBI measurements are supplied to the EKF.', '', '| Case | Position RMS (m) | Velocity RMS (mm/s) | Endpoint position error (m) | Doppler / range |', '|---|---:|---:|---:|---:|']
for row in rows:
    lines.append(f"| {labels[row['case']]} | {row['position_rms_m']:.3f} | {row['velocity_rms_mm_s']:.4f} | {row['endpoint_position_error_m']:.3f} | {row['doppler_observations']} / {row['range_observations']} |")
lines += ['', 'RMS uses the common nominal 600 s diagnostic grid plus the last 60 s endpoint, including tracking gaps and the initial prior. The final point is one minute before the saved truth boundary. Formal uncertainty is conditional on fixed forces and the assumed independent measurement noise; a single run does not establish statistical consistency.', '', 'The ideal-clock case has zero clock uncertainty. The DSAC-inspired case estimates bias/frequency with the existing white-FM model, known deterministic aging and matching process covariance. Its 14-day white-FM continuation is illustrative, not a fitted flight oscillator. The nominal clock reference is TT-equivalent with the existing smooth TDB−TT offset; complete spacecraft proper-time transfer is outside this baseline.', '', 'Tracking windows come from the retained realistic passes, but every new uplink is visibility-checked at its actual ground emission epochs, at both count endpoints. Independent 10 m range and 0.01 Hz counted Doppler errors are assumed; the earlier two-way shared pass biases are not reused. There is no orbital process noise or fitted force scale.', '', 'RMS plots have symmetric ±3σ covariance shading and no separate uncertainty curve. Orbit error and innovation plots are unshaded. Clock plots show component ±3σ bounds. NIS and NEES are diagnostics, not Monte Carlo ensemble tests.', '', 'The one-day half-step check changed solutions by at most 16.6 m (ideal clock) and 40.1 m (DSAC-inspired), below 0.32% of formal position uncertainty. Metre-level conclusions require tighter numerical validation. See `numerical_checks.json`.', '', 'The existing full-arc two-way WLS results are reference context only: those estimates use future observations and different observables/covariance assumptions, so this is not an isolated one-way versus two-way performance comparison.', '', '[DSAC mission results, NASA](https://ntrs.nasa.gov/citations/20230006982) · [Experiment assumptions and reproduction](../../benchmarks/cassini_one_way_ekf/README.md)', '', '![Clock-case comparison](clock_case_comparison.png)', '']
if 'dsac' in results:
    summary, state, innovations, validation = results['dsac']
    along_frequency = validation['endpoint_orbit_clock_correlations']['along_frequency']
    lines += ['The estimated-clock endpoint position uncertainty is '
              f"{state.formal_position_sigma_m.iloc[-1]/1000:.2f} km (root trace of Cartesian covariance). "
              f"Along-track position and clock frequency have correlation {along_frequency:.4f}. "
              'This strong coupling supports orbit–clock ambiguity as the main limitation of this setup; '
              'it does not demonstrate excessive oscillator instability. The declared initial clock calibration priors matter.', '',
              f"Mean NIS divided by observation dimension is {validation['mean_nis_per_dimension']:.4f}; "
              'this describes innovation fit for one realization, not proof of ensemble consistency.', '']
for kind in results:
    lines += [f'### {labels[kind]}', '', f'![Cumulative RMS with ±3σ]({kind}/rms_3sigma.png)', '', f'![Unshaded innovations]({kind}/innovations.png)', '', f'![Innovation noise-scale detail; full data retained above]({kind}/innovations_detail.png)', '', f'![Clock errors]({kind}/clock_errors.png)', '', f'![Consistency diagnostics]({kind}/consistency.png)', '']
(output/'REPORT.md').write_text('\n'.join(lines)+'\n')
