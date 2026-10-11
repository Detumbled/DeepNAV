"""Plot completed VLBI-fit errors and residuals; all covariance propagation is C++."""
from pathlib import Path
import argparse,json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output',type=Path,default=Path('results/cassini_two_way_wls/pass_schedule/with_vlbi'))
folder=parser.parse_args().output.resolve()
s=json.loads((folder/'summary.json').read_text())
orbit=pd.read_csv(folder/'trajectory.csv');u=pd.read_csv(folder/'trajectory_uncertainty.csv')
r=pd.read_csv(folder/'residuals.csv');v=pd.read_csv(folder/'vlbi_residuals.csv')
assert s['converged'] and len(v)==s['vlbi_observations']
np.testing.assert_array_equal(orbit.elapsed_s,u.elapsed_s)
np.testing.assert_allclose(orbit.formal_position_sigma_m,u.formal_position_sigma_m,rtol=1e-10)
np.testing.assert_allclose(orbit.fitted_position_error_m,u.position_error_m,rtol=1e-10,atol=1e-5)
np.testing.assert_allclose(orbit.fitted_velocity_error_mm_s,u.velocity_error_mm_s,rtol=1e-10,atol=1e-5)
np.testing.assert_allclose(u.radial_sigma_m**2+u.along_sigma_m**2+u.cross_sigma_m**2,u.formal_position_sigma_m**2,rtol=1e-12)
assert np.isfinite(u.to_numpy()).all() and (u.filter(like='sigma_m')>0).all().all()
cov=np.loadtxt(folder/'initial_covariance.csv',delimiter=',')
np.testing.assert_allclose(u.formal_position_sigma_m.iloc[0],np.sqrt(np.trace(cov[:3,:3]))*1000,rtol=1e-12)
np.testing.assert_allclose(u.formal_velocity_sigma_mm_s.iloc[0],np.sqrt(np.trace(cov[3:6,3:6]))*1e6,rtol=1e-12)
days=orbit.elapsed_s.to_numpy()/86400
error=orbit.fitted_position_error_m.to_numpy()
cumulative=np.sqrt(np.cumsum(error**2)/np.arange(1,len(error)+1))
np.testing.assert_allclose(cumulative[-1],s['fitted_position_rms_m'],rtol=1e-12)
color='#187c9d';shade='#80b6ca'
def save(fig,name):
    fig.savefig(folder/(name+'.png'),dpi=180);fig.savefig(folder/(name+'.pdf'));plt.close(fig)

fig,axes=plt.subplots(2,1,figsize=(12,8),layout='constrained')
metrics={};analysis=u.copy()
for ax,name,unit,errcol,sigcol in zip(axes,['position','velocity'],['m','mm/s'],['position_error_m','velocity_error_mm_s'],['formal_position_sigma_m','formal_velocity_sigma_mm_s']):
    e=u[errcol].to_numpy();sigma=u[sigcol].to_numpy()
    rms=np.sqrt(np.cumsum(e**2)/np.arange(1,len(e)+1))
    formal_rms=np.sqrt(np.cumsum(sigma**2)/np.arange(1,len(e)+1))
    ax.fill_between(days,-3*formal_rms,3*formal_rms,color=shade,alpha=.3,label='±3σ covariance reference')
    ax.plot(days,rms,color=color,lw=2,label=f'Cumulative 3D error RMS; final {rms[-1]:.3f} {unit}')
    ax.axhline(0,color='black',lw=.7)
    ax.set(title=f'{name.title()} RMS against saved simulation',ylabel=f'RMS ({unit})',xlabel='Elapsed TDB days from 9 October 2004')
    ax.legend(fontsize=9);ax.grid(alpha=.2)
    metrics[name+'_rms']=float(rms[-1]);metrics[name+'_formal_rms']=float(formal_rms[-1])
    analysis['cumulative_'+errcol]=rms;analysis['cumulative_'+sigcol]=formal_rms
fig.suptitle('Cassini · VLBI-augmented WLS · position and velocity RMS\nSymmetric ±3σ covariance reference; RMS remains nonnegative')
save(fig,'vlbi_rms_3sigma')

fig,axes=plt.subplots(3,1,figsize=(13,9),layout='constrained')
# First-order two-way Doppler-to-range-rate conversion at the constant reference frequency.
reference_hz=7.175e9*880/749
rate_scale=-299792.458e6/(2*reference_hz) # mm/s per Hz; receding gives negative Doppler
for ax,data,column,sigma,label,unit,scale in zip(axes,[r,r,v],['postfit_range_m','postfit_doppler_hz','postfit_m'],[s['range_sigma_m'],s['doppler_sigma_hz']*abs(rate_scale),s['vlbi_sigma_m']*1000],['Two-way range','Equivalent two-way range rate (60 s count)','VLBI differential path'],['m','mm/s','mm'],[1,rate_scale,1000]):
    y=data[column].to_numpy()*scale;valid=np.isfinite(y)
    rms=float(np.sqrt(np.mean(y[valid]**2)))
    if column=='postfit_m':
        groups=data.station_one+'/'+data.station_two
    else:groups=data.station
    for group in groups.unique():
        mask=(groups==group).to_numpy()&valid
        ax.scatter(data.elapsed_s[mask]/86400,y[mask],s=12 if unit=='mm' else 5,alpha=.7,label=group)
    ax.axhline(0,color='black',lw=.7)
    ax.set(title=f'{label} · postfit RMS {rms:.4g} {unit}',ylabel=f'Observed − predicted ({unit})',xlabel='Elapsed TDB days from 9 October 2004')
    ax.legend(fontsize=8,ncol=4);ax.grid(alpha=.2)
fig.suptitle('VLBI-augmented fit · postfit observation residuals')
save(fig,'vlbi_residuals_3sigma')
analysis.to_csv(folder/'vlbi_uncertainty_analysis.csv',index=False)
pd.DataFrame(dict(elapsed_s=r.elapsed_s,station=r.station,postfit_range_m=r.postfit_range_m,postfit_equivalent_range_rate_mm_s=r.postfit_doppler_hz*rate_scale)).to_csv(folder/'range_rate_residuals.csv',index=False)
(folder/'VLBI_UNCERTAINTY.md').write_text('''# VLBI fit: residuals and formal uncertainty

The existing eight-parameter WLS solution is unchanged. C++ repropagates its saved estimate and initial covariance using the same force model and STM, exporting both position and velocity covariance traces.

The cumulative position and velocity RMS are sqrt(mean(||error||²)) from the start through each displayed epoch. Their final values are full-arc RMS. Formal RMS uncertainty is sqrt(mean(trace(P))) over the same epochs. Shading is displayed symmetrically from minus to plus three times this formal uncertainty. RMS remains nonnegative; the negative half is a visual covariance reference, not negative RMS. This is not a confidence interval on the measured RMS or a calibrated three-dimensional probability region.

Residual panels contain no uncertainty shading. Equivalent range-rate residual = -c × Doppler residual / (2 × reference frequency), where the constant reference is (880/749) × 7.175 GHz; this first-order conversion preserves the 60-second Doppler count, rather than replacing it with instantaneous geometric range rate. Range samples remain correlated within a pass. The RMS panels retain the shaded covariance reference without a separate formal-uncertainty curve.

[Propagated covariance CSV](trajectory_uncertainty.csv) · [Analysis CSV](vlbi_uncertainty_analysis.csv) · [Range-rate residual CSV](range_rate_residuals.csv)

![Position and velocity RMS](vlbi_rms_3sigma.png)

![Postfit residuals](vlbi_residuals_3sigma.png)
''')
print(json.dumps(metrics,indent=2))
