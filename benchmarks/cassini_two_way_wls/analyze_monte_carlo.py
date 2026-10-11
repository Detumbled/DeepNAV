"""Analyze completed C++ Cassini Monte Carlo exports; no simulation or estimation."""
from pathlib import Path
import argparse
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def analyze(study):
    study = study.resolve()
    expected = int((study/'manifest.txt').read_text().splitlines()[1].split()[0])
    metrics = pd.read_csv(study/'metrics.csv', float_precision='round_trip')
    status = metrics[['trial_id', 'variant', 'status', 'message']].drop_duplicates()
    assert len(status) == 2*expected
    assert set(status.variant) == {'no_vlbi', 'with_vlbi'}
    assert all(set(g.trial_id) == set(range(expected)) for _, g in status.groupby('variant'))
    output = study/'analysis'
    output.mkdir(exist_ok=True)
    status.to_csv(output/'fit_status.csv', index=False)
    variants = ['no_vlbi', 'with_vlbi']
    labels = ['Without VLBI', 'With VLBI']
    colors = ['#c46b35', '#187c9d']
    parameter_rows, epoch_rows, trial_rows = [], [], []
    ensembles = {}
    for variant in variants:
        good = status[(status.variant == variant) & (status.status == 'success')]
        parameters, covariance, epochs = [], [], []
        for trial in good.trial_id:
            folder = study/f'trial_{trial}'/variant
            s = json.loads((folder/'summary.json').read_text())
            assert s['converged'] and s['truth'] == 'simulation' and s['estimated_parameters'] == 8
            assert s['gravity_mask'] == 2 and s['mass_kg'] == 3220
            p = pd.read_csv(folder/'initial_estimate.csv', float_precision='round_trip')
            cov = np.loadtxt(folder/'initial_covariance.csv', delimiter=',')
            np.linalg.cholesky(cov/np.sqrt(np.outer(np.diag(cov), np.diag(cov))))
            e = pd.read_csv(folder/'epoch_error_covariance.csv', float_precision='round_trip')
            if epochs:
                np.testing.assert_array_equal(e.elapsed_s, epochs[0].elapsed_s)
                np.testing.assert_array_equal(p.component, parameters[0].component)
                np.testing.assert_array_equal(p.truth, parameters[0].truth)
            parameters.append(p)
            covariance.append(cov)
            epochs.append(e)
            trial_rows.append(dict(trial_id=int(trial), variant=variant,
                position_rms_m=s['fitted_position_rms_m'], endpoint_error_m=s['fitted_endpoint_error_m'],
                iterations=s['iterations'], range_residual_rms_m=s['postfit_range_rms_m'],
                doppler_residual_rms_hz=s['postfit_doppler_rms_hz'], vlbi_residual_rms_m=s['postfit_vlbi_rms_m']))
        if not epochs:
            continue
        errors = np.stack([(p.estimate-p.truth).to_numpy() for p in parameters])
        covariances = np.stack(covariance)
        for j, name in enumerate(parameters[0].component):
            sigma = np.sqrt(covariances[:, j, j])
            parameter_rows.append(dict(variant=variant, component=name, unit=parameters[0].unit.iloc[j],
                successful_trials=len(epochs), mean_estimate=float(np.mean([p.estimate.iloc[j] for p in parameters])),
                truth=float(parameters[0].truth.iloc[j]), bias=float(errors[:, j].mean()),
                empirical_sigma=float(errors[:, j].std(ddof=1)) if len(epochs)>1 else np.nan,
                formal_sigma_rms=float(np.sqrt(np.mean(sigma**2))),
                coverage_1sigma=float(np.mean(abs(errors[:, j]) <= sigma)),
                coverage_2sigma=float(np.mean(abs(errors[:, j]) <= 2*sigma)),
                coverage_3sigma=float(np.mean(abs(errors[:, j]) <= 3*sigma))))
        error_columns = ['dx_m', 'dy_m', 'dz_m', 'dvx_mm_s', 'dvy_mm_s', 'dvz_mm_s']
        state_errors = np.stack([e[error_columns].to_numpy() for e in epochs])
        matrices = np.stack([e[[f'cov_{i}{j}' for i in range(6) for j in range(6)]].to_numpy().reshape(-1, 6, 6) for e in epochs])
        rows = []
        for k, elapsed in enumerate(epochs[0].elapsed_s):
            error = state_errors[:, k]
            cov = matrices[:, k]
            sigma = np.sqrt(np.diagonal(cov, axis1=1, axis2=2))
            normalized = cov/(sigma[:, :, None]*sigma[:, None, :])
            residual = error/sigma
            root = np.linalg.cholesky(normalized)
            solved = np.linalg.solve(root, residual[..., None])[..., 0]
            mean = error.mean(axis=0)
            rows.append(dict(variant=variant, elapsed_s=float(elapsed), successful_trials=len(epochs),
                position_rms_m=float(np.sqrt(np.mean(np.sum(error[:, :3]**2, axis=1)))),
                velocity_rms_mm_s=float(np.sqrt(np.mean(np.sum(error[:, 3:]**2, axis=1)))),
                position_bias_m=float(np.linalg.norm(mean[:3])), velocity_bias_mm_s=float(np.linalg.norm(mean[3:])),
                formal_position_rms_m=float(np.sqrt(np.trace(cov[:, :3, :3], axis1=1, axis2=2).mean())),
                formal_velocity_rms_mm_s=float(np.sqrt(np.trace(cov[:, 3:, 3:], axis1=1, axis2=2).mean())),
                mean_state_nees=float(np.mean(np.sum(solved**2, axis=1)))))
        epoch_rows.extend(rows)
        ensembles[variant] = pd.DataFrame(rows)
    trials = pd.DataFrame(trial_rows)
    pd.DataFrame(parameter_rows).to_csv(output/'parameter_statistics.csv', index=False)
    pd.DataFrame(epoch_rows).to_csv(output/'epoch_statistics.csv', index=False)
    trials.to_csv(output/'trial_statistics.csv', index=False)
    summary = dict(requested_trials=expected, variants={}, paired_successful_trials=0)
    for variant in variants:
        counts = status[status.variant == variant].status.value_counts().to_dict()
        entry = dict(success=int(counts.get('success', 0)), not_converged=int(counts.get('not_converged', 0)),
                     failed=int(counts.get('failed', 0)))
        if not trials.empty:
            t = trials[trials.variant == variant]
            if len(t):
                entry.update(position_rms_m_mean=float(t.position_rms_m.mean()),
                    position_rms_m_median=float(t.position_rms_m.median()),
                    position_rms_m_p05=float(t.position_rms_m.quantile(.05)),
                    position_rms_m_p95=float(t.position_rms_m.quantile(.95)))
        summary['variants'][variant] = entry
    def save(fig, name):
        fig.savefig(output/(name+'.png'), dpi=180)
        fig.savefig(output/(name+'.pdf'))
        plt.close(fig)
    figures = []
    if ensembles:
        fig, axes = plt.subplots(2, 2, figsize=(13, 8), layout='constrained')
        for variant, label, color in zip(variants, labels, colors):
            if variant not in ensembles:
                continue
            e = ensembles[variant]
            days = e.elapsed_s/86400
            for ax, column, formal, ylabel in [
                (axes[0, 0], 'position_rms_m', 'formal_position_rms_m', '3D position RMS (m)'),
                (axes[1, 0], 'velocity_rms_mm_s', 'formal_velocity_rms_mm_s', '3D velocity RMS (mm/s)')]:
                ax.plot(days, e[column], color=color, label=label)
                ax.plot(days, e[formal], '--', color=color, alpha=.65, label=label+' formal')
                ax.set_ylabel(ylabel)
            axes[0, 1].plot(days, e.position_bias_m, color=color, label=label)
            axes[1, 1].plot(days, e.velocity_bias_mm_s, color=color, label=label)
        axes[0, 1].set_ylabel('Magnitude of ensemble position bias (m)')
        axes[1, 1].set_ylabel('Magnitude of ensemble velocity bias (mm/s)')
        for ax in axes.flat:
            ax.grid(alpha=.2)
            ax.set_xlabel('Elapsed TDB days from 9 October 2004')
            ax.legend(fontsize=8)
        fig.suptitle('Cassini two-week Monte Carlo · independent noise realizations · fixed simulated truth')
        save(fig, 'ensemble_errors')
        figures.append('ensemble_errors')
    if not trials.empty:
        paired = trials.pivot(index='trial_id', columns='variant', values='position_rms_m').dropna()
        if set(variants).issubset(paired.columns):
            paired['improvement_factor'] = paired.no_vlbi/paired.with_vlbi
            paired.to_csv(output/'paired_improvement.csv')
            summary.update(paired_successful_trials=len(paired),
                median_improvement_factor=float(paired.improvement_factor.median()),
                fraction_improved_with_vlbi=float((paired.with_vlbi < paired.no_vlbi).mean()))
            fig, axes = plt.subplots(1, 2, figsize=(12, 5), layout='constrained')
            for variant, label, color in zip(variants, labels, colors):
                data = trials[trials.variant == variant].position_rms_m
                axes[0].hist(data, bins=max(5, int(np.sqrt(len(data)))), alpha=.55, color=color, label=label)
            axes[0].set(xlabel='Full-arc position RMS (m)', ylabel='Successful fits')
            axes[0].legend()
            axes[1].scatter(paired.no_vlbi, paired.with_vlbi, s=18, color=colors[1])
            maximum = max(paired.no_vlbi.max(), paired.with_vlbi.max())
            axes[1].plot([0, maximum], [0, maximum], 'k--', label='Equal RMS')
            axes[1].set(xlabel='Without VLBI RMS (m)', ylabel='With VLBI RMS (m)')
            axes[1].legend()
            for ax in axes:
                ax.grid(alpha=.2)
            fig.suptitle('Paired fits · identical range and Doppler observations in each trial')
            save(fig, 'paired_improvement')
            figures.append('paired_improvement')
    (output/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    lines = ['# Cassini two-week Monte Carlo', '', f'{expected} requested independent noise realizations; two estimator variants.', '',
        '| Variant | Successful | Not converged | Failed | Median arc position RMS (m) | 5–95% range (m) |',
        '|---|---:|---:|---:|---:|---:|']
    for variant, label in zip(variants, labels):
        e = summary['variants'][variant]
        median = f"{e['position_rms_m_median']:.2f}" if e['success'] else '—'
        spread = f"{e['position_rms_m_p05']:.2f}–{e['position_rms_m_p95']:.2f}" if e['success'] else '—'
        lines.append(f"| {label} | {e['success']} | {e['not_converged']} | {e['failed']} | {median} | {spread} |")
    if summary['paired_successful_trials']:
        lines += ['', f"Paired successful trials: {summary['paired_successful_trials']}. Median per-trial VLBI improvement factor: {summary['median_improvement_factor']:.3f}. Fraction with lower RMS using VLBI: {summary['fraction_improved_with_vlbi']:.1%}."]
    lines += ['', 'Truth, geometry, dynamical model, starting guess and priors remain fixed; only synthetic measurement errors change. Range errors include independent sample noise and a shared pass offset. Both variants use identical radio observations. All generation, propagation, estimation and numerical exports are C++; Python performs analysis and plotting.', '',
        'Statistics above condition on successful fits; all failure and nonconvergence counts are retained separately. Epoch RMS uses errors across realizations at each saved epoch, rather than treating epochs as independent trials. Bias is the magnitude of the ensemble mean error. Formal RMS is the square root of the ensemble-average covariance trace.', '',
        'The parameter table reports empirical scatter, bias, formal uncertainty and component coverage. NEES and coverage are descriptive for this fixed-truth experiment with an offset prior; they do not establish unconditional Bayesian chi-square consistency. The existing shaded single-fit plots are unchanged.', '',
        '[Fit status](fit_status.csv) · [Trial statistics](trial_statistics.csv) · [Parameter statistics](parameter_statistics.csv) · [Epoch statistics](epoch_statistics.csv)']
    for name in figures:
        lines += ['', f'![{name.replace("_", " ")}]({name}.png)']
    (output/'REPORT.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--study', type=Path, required=True)
    analyze(parser.parse_args().study)
