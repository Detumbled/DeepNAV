"""Compare completed C++ frequency-prior runs with identical simulated observations."""
from pathlib import Path
import argparse
import hashlib
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output', type=Path, default=ROOT/'results/cassini_one_way_ekf/frequency_prior_sweep')
output = parser.parse_args().output.resolve()


def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def save(fig, name):
    fig.savefig(output/f'{name}.png', dpi=180)
    fig.savefig(output/f'{name}.pdf')
    plt.close(fig)


cases = []
reference_observations = None
reference_clock = None
reference_prior = None
reference_summary = None
rows = []
for label, sigma, color in [('1e-10', 1e-10, '#8a63ad'), ('1e-11', 1e-11, '#dd8452'), ('1e-12', 1e-12, '#187c9d')]:
    folder = output/label/'dsac'
    summary = json.loads((folder/'summary.json').read_text())
    prior = pd.read_csv(folder/'initial_prior.csv', index_col='parameter', float_precision='round_trip')
    observations = pd.read_csv(folder/'observations.csv', float_precision='round_trip')
    clock = pd.read_csv(folder/'clock_truth.csv')
    state = pd.read_csv(folder/'trajectory.csv')
    innovations = pd.read_csv(folder/'innovations.csv')
    covariance = pd.read_csv(folder/'covariance.csv')
    assert prior.loc['fractional_frequency', 'standard_deviation'] == sigma
    assert prior.loc['fractional_frequency', 'mean'] == 0
    assert prior.loc['fractional_frequency', 'truth'] == 2e-12
    assert prior.loc['clock_bias_s', 'standard_deviation'] == 1e-5
    assert summary['duration_days'] == 14 and summary['estimated_parameters'] == 8
    assert innovations.accepted.eq(1).all()
    assert summary['maximum_half_stencil_noise_fraction'] < .01
    if reference_observations is None:
        reference_observations, reference_clock = observations, clock
        reference_prior, reference_summary = prior, summary
    else:
        pd.testing.assert_frame_equal(reference_observations, observations, check_exact=True)
        pd.testing.assert_frame_equal(reference_clock, clock, check_exact=True)
        pd.testing.assert_frame_equal(reference_prior.drop(index='fractional_frequency'), prior.drop(index='fractional_frequency'), check_exact=True)
        for key in ['clock_case', 'duration_days', 'seed', 'noise_enabled', 'truth_source', 'clock_configuration', 'clock_drift_per_s', 'clock_q_bias_s', 'clock_q_frequency_per_s', 'range_sigma_m', 'doppler_sigma_hz', 'initial_guess_scale', 'max_step_s', 'doppler_observations', 'range_observations', 'diagnostic_samples']:
            assert summary[key] == reference_summary[key], key
    np.testing.assert_allclose(summary['position_rms_m'], np.sqrt(np.mean(state.position_error_m**2)), rtol=1e-12)
    np.testing.assert_allclose(summary['velocity_rms_mm_s'], np.sqrt(np.mean(state.velocity_error_mm_s**2)), rtol=1e-12)
    np.testing.assert_array_equal(state.elapsed_s, covariance.elapsed_s)
    p = covariance.iloc[:, 1:].to_numpy().reshape(-1, 8, 8)
    scales = np.array([2, 2, 2, 1e-4, 1e-4, 1e-4, 1e-5, sigma])
    normalized = p/scales[None, :, None]/scales[None, None, :]
    np.testing.assert_allclose(normalized, normalized.transpose(0, 2, 1), atol=1e-9, rtol=1e-9)
    assert np.linalg.eigvalsh(normalized).min() >= -1e-9
    last = state.iloc[-1]
    r = last[['truth_x_km', 'truth_y_km', 'truth_z_km']].to_numpy(float)
    v = last[['truth_vx_km_s', 'truth_vy_km_s', 'truth_vz_km_s']].to_numpy(float)
    radial = r/np.linalg.norm(r)
    cross = np.cross(r, v); cross /= np.linalg.norm(cross)
    along = np.cross(cross, radial)
    correlation = (along@p[-1, :3, 7])/np.sqrt((along@p[-1, :3, :3]@along)*p[-1, 7, 7])
    row = {'frequency_prior_sigma': sigma,
        'position_rms_m': summary['position_rms_m'], 'velocity_rms_mm_s': summary['velocity_rms_mm_s'],
        'endpoint_position_error_m': summary['endpoint_position_error_m'],
        'endpoint_velocity_error_mm_s': summary['endpoint_velocity_error_mm_s'],
        'endpoint_formal_position_sigma_m': float(last.formal_position_sigma_m),
        'endpoint_frequency_error': float(last.clock_frequency_error),
        'endpoint_frequency_sigma': float(last.clock_frequency_sigma),
        'endpoint_frequency_error_over_sigma': float(last.clock_frequency_error/last.clock_frequency_sigma),
        'endpoint_bias_error_s': float(last.clock_bias_error_s),
        'endpoint_bias_sigma_s': float(last.clock_bias_sigma_s),
        'endpoint_along_frequency_correlation': float(correlation),
        'mean_nis_per_dimension': float((innovations.nis/innovations.dimension).mean()),
        'mean_nees_per_state': float(state.joint_nees.mean()/8),
        'endpoint_nees': float(last.joint_nees),
        'doppler_observations': len(observations), 'range_observations': int(observations.range_active.sum())}
    rows.append(row)
    cases.append((label, color, state, innovations))

pd.DataFrame(rows).to_csv(output/'comparison.csv', index=False)
provenance = {'identical_observations': True, 'identical_clock_truth': True,
    'changed_input': 'Initial fractional-frequency prior standard deviation only',
    'frequency_prior_mean': 0, 'true_initial_frequency': 2e-12,
    'clock_bias_prior_sigma_s': 1e-5,
    'baseline_reused': 'results/cassini_one_way_ekf/dsac',
    'input_sha256': {label: {name: digest(output/label/'dsac'/name) for name in ['observations.csv', 'clock_truth.csv', 'initial_prior.csv']} for label, _, _, _ in cases}}
(output/'validation.json').write_text(json.dumps(provenance, indent=2)+'\n')

fig, axes = plt.subplots(2, 1, figsize=(11, 7), sharex=True, layout='constrained')
for label, color, state, _ in cases:
    axes[0].plot(state.elapsed_s/86400, state.position_error_m/1000, color=color, label=f'σ_y = {label}')
    axes[1].plot(state.elapsed_s/86400, state.velocity_error_mm_s, color=color, label=f'σ_y = {label}')
axes[0].set_ylabel('3D position error (km)')
axes[1].set(ylabel='3D velocity error (mm/s)', xlabel='Elapsed TDB days')
for ax in axes:
    ax.grid(alpha=.2)
    ax.legend()
fig.suptitle('Clock-frequency prior sweep · identical one-way observations and clock truth')
save(fig, 'orbit_errors')

fig, axes = plt.subplots(3, 2, figsize=(13, 11), sharex=True, layout='constrained')
for row_axes, (label, color, state, _) in zip(axes, cases):
    for ax, error, sigma, factor, title, unit in zip(row_axes,
            ['position_error_m', 'velocity_error_mm_s'], ['formal_position_sigma_m', 'formal_velocity_sigma_mm_s'],
            [1/1000, 1], ['Position', 'Velocity'], ['km', 'mm/s']):
        n = np.arange(1, len(state)+1)
        rms = np.sqrt(np.cumsum(state[error].to_numpy()**2)/n)*factor
        bound = 3*np.sqrt(np.cumsum(state[sigma].to_numpy()**2)/n)*factor
        ax.fill_between(state.elapsed_s/86400, -bound, bound, alpha=.16, color=color, label='±3σ covariance reference')
        ax.plot(state.elapsed_s/86400, rms, color=color, label='Cumulative 3D RMS error')
        ax.set(ylabel=f'RMS ({unit})', title=f'{title} · initial σ_y = {label}')
        ax.grid(alpha=.2)
        ax.legend(fontsize=8)
for ax in axes[-1]:
    ax.set_xlabel('Elapsed TDB days')
fig.suptitle('Frequency-prior sensitivity · RMS includes initial convergence and tracking gaps')
save(fig, 'rms_3sigma')

fig, axes = plt.subplots(2, 1, figsize=(11, 7), sharex=True, layout='constrained')
for label, color, state, _ in cases:
    axes[0].plot(state.elapsed_s/86400, state.clock_bias_error_s*1e6, color=color, label=f'σ_y = {label}')
    axes[1].plot(state.elapsed_s/86400, state.clock_frequency_error*1e12, color=color, label=f'σ_y = {label}')
axes[0].set_ylabel('Clock bias error (µs)')
axes[1].set(ylabel='Fractional frequency error (×10⁻¹²)', xlabel='Elapsed TDB days')
for ax in axes:
    ax.grid(alpha=.2)
    ax.legend()
fig.suptitle('Clock recovery · unchanged zero prior mean and simulated initial offset 2×10⁻¹²')
save(fig, 'clock_errors')

lines = ['# One-way EKF clock-frequency prior sweep', '',
    'Only the initial fractional-frequency prior standard deviation changes. Observations, clock truth, force model, orbit prior, clock-bias prior (10 µs), frequency prior mean (zero), process noise and noise seeds are identical. The original 10⁻¹⁰ result is reused with its numerical outputs and provenance preserved.', '',
    '| Initial frequency σ | Position RMS (km) | Velocity RMS (mm/s) | Endpoint error (km) | Endpoint formal position σ (km) | Along/frequency correlation | Mean NIS/dimension |',
    '|---|---:|---:|---:|---:|---:|---:|']
for row in rows:
    lines.append(f"| {row['frequency_prior_sigma']:.0e} | {row['position_rms_m']/1000:.3f} | {row['velocity_rms_mm_s']:.3f} | {row['endpoint_position_error_m']/1000:.3f} | {row['endpoint_formal_position_sigma_m']/1000:.3f} | {row['endpoint_along_frequency_correlation']:.4f} | {row['mean_nis_per_dimension']:.4f} |")
tightest = rows[-1]
best = min(rows, key=lambda r: r['position_rms_m'])
reduction = 100*(1-best['position_rms_m']/rows[0]['position_rms_m'])
lines += ['', f"The lowest position RMS in this realization uses σ_y={best['frequency_prior_sigma']:.0e}, reducing RMS by {reduction:.1f}% relative to the unchanged baseline. This is a sensitivity experiment, not evidence that this prior describes an independently calibrated flight clock.", '',
    f"The tightest prior reduces the endpoint along-track/frequency correlation to {tightest['endpoint_along_frequency_correlation']:.4f}, but its endpoint frequency error remains {tightest['endpoint_frequency_error']:.3e} ({tightest['endpoint_frequency_error_over_sigma']:.2f} posterior σ). It largely holds frequency near the zero prior mean rather than recovering the injected offset.", '',
    'The true initial frequency is 2×10⁻¹² while every prior mean remains zero: the tightest prior is initially two standard deviations from truth. Smaller covariance alone is not success; inspect clock recovery, orbital errors and consistency statistics together. Results from one shared-noise realization do not establish average performance or statistical calibration.', '',
    'No SRP or thermal coefficients were added or estimated. RMS is evaluated on the same nominal 600 s diagnostic grid including gaps, initial convergence and the final minute-before-boundary point. Formal RMS shading is a symmetric covariance reference, not a confidence interval on the nonnegative RMS. Error plots are unshaded.', '',
    '[Comparison CSV](comparison.csv) · [Input identity checks](validation.json) · [Benchmark instructions](../../../benchmarks/cassini_one_way_ekf/README.md)', '',
    '![Orbit errors](orbit_errors.png)', '', '![Cumulative RMS and ±3σ](rms_3sigma.png)', '', '![Clock errors](clock_errors.png)', '']
(output/'REPORT.md').write_text('\n'.join(lines)+'\n')
print(pd.DataFrame(rows).to_string(index=False))
