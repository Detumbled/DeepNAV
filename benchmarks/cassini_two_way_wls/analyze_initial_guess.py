"""Compare C++ WLS outputs from original and larger numerical starting guesses."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
root=Path(__file__).resolve().parents[2]/'results/cassini_two_way_wls'
output=root/'relaxed_initial_guess'
lines=['# Larger numerical starting guess','', 'Only the numerical starting position/velocity offsets are multiplied by ten: [10, −8, 5] km and [0.10, −0.08, 0.06] m/s. Force scales start at 1. Prior means and uncertainties, simulated truth, force model, observations and noise remain identical to the original runs. The convergence test remains an undamped joint correction below 0.1 posterior σ, with at least two iterations.', '', '| Fit | Original iterations | Larger-guess iterations | Original RMS (m) | Larger-guess RMS (m) |','|---|---:|---:|---:|---:|']
fig,axes=plt.subplots(1,2,figsize=(12,4.5),layout='constrained')
metrics={}
notes=[]
for ax,name,label in zip(axes,['no_vlbi','with_vlbi'],['Range + Doppler','Range + Doppler + VLBI']):
    folders=[root/name,output/name]
    original=json.loads((folders[0]/'summary.json').read_text())
    if not (folders[1]/'summary.json').exists():
        ax.text(.5,.5,'No converged output',ha='center',va='center',transform=ax.transAxes);ax.set_title(label)
        lines.append(f"| {label} | {original['iterations']} | No converged output | {original['fitted_position_rms_m']:.3f} | — |")
        metrics[name]={'converged':False};continue
    modified=json.loads((folders[1]/'summary.json').read_text())
    assert modified['converged'] and modified['initial_guess_scale']==10
    np.testing.assert_array_equal(original['initial_prior_mean'],modified['initial_prior_mean'])
    np.testing.assert_array_equal(original['initial_prior_sigma'],modified['initial_prior_sigma'])
    for filename in ['observations.csv','vlbi_observations.csv']:
        if (folders[0]/filename).exists():pd.testing.assert_frame_equal(pd.read_csv(folders[0]/filename),pd.read_csv(folders[1]/filename),check_exact=True)
    for folder,text,color in zip(folders,['Original guess','10× larger starting error'],['#c46b35','#187c9d']):
        orbit=pd.read_csv(folder/'trajectory.csv');ax.plot(orbit.elapsed_s/86400,orbit.fitted_position_error_m,label=text,color=color)
    ax.set(title=label,xlabel='Elapsed TDB days',ylabel='WLS − simulation 3D error (m)');ax.grid(alpha=.2);ax.legend(fontsize=9)
    sig=np.asarray(modified['initial_prior_sigma']);delta=(np.asarray(modified['initial_estimated_state'])-original['initial_estimated_state'])/sig
    covariance=np.loadtxt(folders[1]/'initial_covariance.csv',delimiter=',')/np.outer(sig,sig)
    distance=float(np.sqrt(delta@np.linalg.solve(covariance,delta)))
    metrics[name]={'converged':True,'original_rms_m':original['fitted_position_rms_m'],'larger_guess_rms_m':modified['fitted_position_rms_m'],'original_iterations':original['iterations'],'larger_guess_iterations':modified['iterations'],'final_undamped_step_sigma':modified['undamped_posterior_step_sigma'],'estimate_difference_posterior_sigma':distance}
    lines.append(f"| {label} | {original['iterations']} | {modified['iterations']} | {original['fitted_position_rms_m']:.3f} | {modified['fitted_position_rms_m']:.3f} |")
    notes+=['',f"{label}: final undamped correction {modified['undamped_posterior_step_sigma']:.6g} posterior σ. Difference between the two final eight-parameter estimates: {distance:.6g} posterior σ. Finite convergence tolerance and counted-Doppler numerical precision permit a small dependence on the starting point.",'']
fig.suptitle('Cassini starting-guess sensitivity · fixed prior and identical observations')
fig.savefig(output/'initial_guess_comparison.png',dpi=180);fig.savefig(output/'initial_guess_comparison.pdf');plt.close(fig)
lines+=notes
lines+=['','[Original comparison](../COMPARISON.md) · [Larger-guess comparison](COMPARISON.md)','', '![Starting-guess comparison](initial_guess_comparison.png)','']
(output/'INITIAL_GUESS.md').write_text('\n'.join(lines)+'\n');(output/'initial_guess_summary.json').write_text(json.dumps(metrics,indent=2)+'\n');print(json.dumps(metrics,indent=2))
