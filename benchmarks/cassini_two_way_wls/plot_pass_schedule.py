"""Plot completed C++ pass-schedule fits and analyze shared range errors."""
from pathlib import Path
import json,argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,default=Path('results/cassini_two_way_wls/pass_schedule'));root=parser.parse_args().output.resolve()
folders=[root/'no_vlbi',root/'with_vlbi'];summaries=[json.loads((p/'summary.json').read_text()) for p in folders]
passes=[pd.read_csv(p/'passes.csv') for p in folders];pd.testing.assert_frame_equal(*passes,check_exact=True)
p=passes[0];colors={'DSS-14':'#2878b5','DSS-43':'#dd8452','DSS-63':'#8a63ad'}
fig,axes=plt.subplots(2,2,figsize=(13,8),layout='constrained')
for station,color in colors.items():
 d=p[p.station==station];axes[0,0].barh(d.pass_id,(d.end_elapsed_s-d.start_elapsed_s)/86400,left=d.start_elapsed_s/86400,color=color,label=station,height=.65)
axes[0,0].set(xlabel='Elapsed TDB days from 9 October 2004',ylabel='Pass ID',title='Visibility-qualified non-overlapping tracking passes');axes[0,0].legend(fontsize=9)
axes[0,1].bar(p.pass_id,p.minimum_elevation_deg,color=[colors[a] for a in p.station]);axes[0,1].axhline(10,color='black',ls='--',label='10° cutoff');axes[0,1].set(xlabel='Pass ID',ylabel='Minimum elevation (deg)',title='Both legs and both count endpoints');axes[0,1].legend()
for folder,label,color in zip(folders,['Without VLBI','With VLBI'],['#c46b35','#187c9d']):
 d=pd.read_csv(folder/'pass_range_analysis.csv');axes[1,0].plot(d.pass_id,d.postfit_range_mean_m,'o-',label=label,color=color)
axes[1,0].plot(p.pass_id,p.true_range_bias_m,'kx',label='Injected shared offset (diagnostic only)');axes[1,0].set(xlabel='Pass ID',ylabel='Pass-mean range residual (m)',title='Shared offsets persist across range samples');axes[1,0].legend(fontsize=8)
s=summaries[0];n=np.arange(1,101);axes[1,1].plot(n,np.sqrt(s['range_pass_bias_sigma_m']**2+s['range_independent_sigma_m']**2/n),label='Correlated pass mean',color='#187c9d');axes[1,1].plot(n,s['range_sigma_m']/np.sqrt(n),'--',label='Treating marginal errors as independent',color='#c46b35');axes[1,1].set(xlabel='Range samples in a pass',ylabel='Pass-mean uncertainty (m)',title='More samples cannot remove a shared bias');axes[1,1].legend(fontsize=8)
for ax in axes.flat:ax.grid(alpha=.2)
fig.suptitle('Cassini pass experiment · 1 m independent range noise + 10 m shared offset per pass')
fig.savefig(root/'pass_schedule.png',dpi=180);fig.savefig(root/'pass_schedule.pdf');plt.close(fig)
rows=[]
previous=root.parent/'relaxed_initial_guess'
for name,label,folder,new in zip(['no_vlbi','with_vlbi'],['Without VLBI','With VLBI'],folders,summaries):
 row=dict(fit=label,position_rms_m=new['fitted_position_rms_m'],endpoint_error_m=new['fitted_endpoint_error_m'],formal_position_rms_m=float(np.sqrt(np.mean(pd.read_csv(folder/'trajectory.csv').formal_position_sigma_m**2))),range_residual_rms_m=new['postfit_range_rms_m'],doppler_residual_rms_hz=new['postfit_doppler_rms_hz'],Cs=new['initial_estimated_state'][6],Cs_sigma=new['initial_posterior_sigma'][6],thermal_scale=new['initial_estimated_state'][7],thermal_sigma=new['initial_posterior_sigma'][7],iterations=new['iterations'])
 if (previous/name/'summary.json').exists():row['previous_schedule_position_rms_m']=json.loads((previous/name/'summary.json').read_text())['fitted_position_rms_m']
 rows.append(row)
pd.DataFrame(rows).to_csv(root/'pass_experiment_metrics.csv',index=False)
lines=['# Scheduled-pass Cassini experiment','',f"{len(p)} passes, {int(p.doppler_points.sum())} Doppler observations and {int(p.range_points.sum())} range observations. Every retained pass has its midpoint in a different elapsed day, can cross midnight, and does not overlap adjacent passes. Station pass counts: {p.station.value_counts().to_dict()}.",'', 'All numerical generation, covariance whitening, propagation and WLS estimation are C++. Python verifies and plots the completed exports. Truth remains the saved simulation. The 10× starting error, original priors and 0.1 posterior-σ convergence criterion remain fixed.','', 'Range covariance for each pass is R=σ_ind²I+σ_shared²11ᵀ, with assumed σ_ind=1 m and σ_shared=10 m. The shared offset is independently drawn once per pass. Whitening marginalizes its Gaussian uncertainty without estimating extra states; WLS is not told the injected offsets. The pass-mean uncertainty for 18 samples is approximately 10.003 m. Doppler remains independent 0.01 Hz noise and VLBI remains independent 1 cm path noise. These values are synthetic assumptions rather than a fitted flight error budget.','', '| Result | Without VLBI | With 60 VLBI |','|---|---:|---:|']
for title,key in [('3D position RMS (m)','position_rms_m'),('Endpoint error (m)','endpoint_error_m'),('Formal position RMS (m)','formal_position_rms_m'),('Range residual RMS (m)','range_residual_rms_m'),('Doppler residual RMS (Hz)','doppler_residual_rms_hz'),('SRP scale C_s','Cs'),('SRP scale formal σ','Cs_sigma'),('Thermal scale λ','thermal_scale'),('Thermal scale formal σ','thermal_sigma'),('WLS iterations','iterations')]:lines.append(f"| {title} | {rows[0][key]:.6g} | {rows[1][key]:.6g} |")
lines+=['', 'Actual DSN allocation logs, antenna pointing restrictions and correlated Doppler media errors are not reconstructed.','', '[Paired trajectory comparison](COMPARISON.md) · [Observation CSV](with_vlbi/observations.csv) · [Tracking passes](with_vlbi/passes.csv) · [Eight-parameter estimate](with_vlbi/initial_estimate.csv) · [Interactive 3D comparison](vlbi_trajectory_3d.html)','', '![Pass schedule and common-range-error diagnostics](pass_schedule.png)','', '![Trajectory error comparison](vlbi_comparison.png)','']
(root/'PASS_EXPERIMENT.md').write_text('\n'.join(lines)+'\n');print(json.dumps(rows,indent=2))
