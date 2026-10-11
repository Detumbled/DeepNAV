"""Analyze paired C++ eight-parameter WLS CSVs; generate scientific comparison figures."""
from pathlib import Path
import argparse,json,hashlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=Path('results/cassini_two_way_wls/pass_schedule'))
    root=parser.parse_args().output.resolve()
    folders=[root/'no_vlbi',root/'with_vlbi']
    summaries=[json.loads((p/'summary.json').read_text()) for p in folders]
    orbits=[pd.read_csv(p/'trajectory.csv') for p in folders]
    observations=[pd.read_csv(p/'observations.csv') for p in folders]
    pd.testing.assert_frame_equal(*observations,check_exact=True)
    for key in ['start_utc','duration_days','mass_kg','area_m2','cr','gravity_mask','thermal_reference_vector_km_s2','initial_truth_state','initial_seed_state','initial_prior_mean','initial_prior_sigma','range_sigma_m','doppler_sigma_hz','noise_seed']:
        assert summaries[0][key]==summaries[1][key],key
    if summaries[0].get('pass_schedule',False):
        for key in ['range_independent_sigma_m','range_pass_bias_sigma_m','range_cadence_s','range_observations','doppler_observations','pass_count']:
            assert summaries[0][key]==summaries[1][key]
        pd.testing.assert_frame_equal(pd.read_csv(folders[0]/'passes.csv'),pd.read_csv(folders[1]/'passes.csv'),check_exact=True)
    assert all(s['converged'] for s in summaries)
    assert summaries[0]['vlbi_observations']==0 and summaries[1]['vlbi_observations']==60
    truthcols=[f'truth_{a}_km' for a in 'xyz']
    np.testing.assert_array_equal(orbits[0][truthcols],orbits[1][truthcols])
    np.testing.assert_array_equal(orbits[0].elapsed_s,orbits[1].elapsed_s)
    labels=['Range + Doppler','Range + Doppler + 60 VLBI'];colors=['#c46b35','#187c9d']
    days=orbits[0].elapsed_s/86400
    def save(fig,name):
        fig.savefig(root/(name+'.png'),dpi=180);fig.savefig(root/(name+'.pdf'));plt.close(fig)
    fig,axes=plt.subplots(2,2,figsize=(13,8),layout='constrained')
    columns=['fitted_position_error_m','radial_error_m','along_error_m','cross_error_m']
    for ax,column,title in zip(axes.flat,columns,['3D position error','Radial error','Along-track error','Cross-track error']):
        for data,label,color in zip(orbits,labels,colors):
            ax.plot(days,data[column],label=label,color=color)
            if column==columns[0]:ax.plot(days,data.formal_position_sigma_m,'--',color=color,alpha=.65,label=label+' formal √trace(Prr)')
        ax.set(xlabel='Elapsed TDB days from 9 October 2004',ylabel='Meters',title=title);ax.grid(alpha=.2)
    axes[0,0].legend(fontsize=8)
    fig.suptitle('Cassini WLS − saved simulated trajectory · identical dynamics, priors and two-way data')
    save(fig,'vlbi_comparison')
    residuals=[pd.read_csv(p/'residuals.csv') for p in folders]
    vlbi=pd.read_csv(folders[1]/'vlbi_observations.csv');vr=pd.read_csv(folders[1]/'vlbi_residuals.csv')
    fig,axes=plt.subplots(3,1,figsize=(12,9),layout='constrained')
    for d,label,color in zip(residuals,labels,colors):
        axes[0].plot(d.elapsed_s/86400,d.postfit_range_m,'.',ms=2,alpha=.45,color=color,label=label)
        axes[1].plot(d.elapsed_s/86400,d.postfit_doppler_hz,'.',ms=2,alpha=.45,color=color)
    for pair,d in vr.groupby(['station_one','station_two']):axes[2].plot(d.elapsed_s/86400,d.postfit_m*1000,'o',ms=4,label='/'.join(pair))
    axes[2].axhline(10,color='gray',ls='--',label='Assumed ±1σ = 10 mm');axes[2].axhline(-10,color='gray',ls='--')
    for ax,y in zip(axes,['Range residual (m)','Doppler residual (Hz)','VLBI residual (mm)']):ax.set(xlabel='Elapsed TDB days',ylabel=y);ax.grid(alpha=.2)
    axes[0].legend(fontsize=9);axes[2].legend(fontsize=9,ncol=2)
    fig.suptitle('Postfit observation − prediction residuals · synthetic observation errors')
    save(fig,'vlbi_residual_comparison')
    fig,axes=plt.subplots(1,3,figsize=(14,4),layout='constrained')
    for j,name in enumerate(['C_s','λ']):
        ax=axes[j]
        for i,(s,label,color) in enumerate(zip(summaries,labels,colors)):
            ax.errorbar(i,s['initial_estimated_state'][6+j],yerr=s['initial_posterior_sigma'][6+j],fmt='o',capsize=5,color=color)
        ax.axhline(1,color='black',ls='--',label='Truth = 1');ax.set_xticks([0,1],['Two-way','Two-way + VLBI']);ax.set(ylabel='Estimate ± formal 1σ',title=name);ax.grid(alpha=.2);ax.legend()
    for d,label,color in zip(orbits,labels,colors):axes[2].plot(days,d.formal_position_sigma_m,label=label,color=color)
    axes[2].set(xlabel='Elapsed TDB days',ylabel='Formal √trace(Prr) (m)',title='Position uncertainty');axes[2].legend(fontsize=8);axes[2].grid(alpha=.2)
    save(fig,'vlbi_parameter_comparison')
    truth=orbits[0][truthcols].to_numpy();fits=[d[[f'fitted_{a}_km' for a in 'xyz']].to_numpy() for d in orbits]
    fig=plt.figure(figsize=(13,6),layout='constrained');a=fig.add_subplot(121,projection='3d');b=fig.add_subplot(122,projection='3d')
    a.plot(*(truth/1e6).T,color='#303c45',label='Saved simulation')
    for fit,label,color in zip(fits,labels,colors):a.plot(*(fit/1e6).T,label=label,color=color,ls='--');b.plot(*((fit-truth)*1000).T,label=label,color=color)
    a.set(xlabel='X (10⁶ km)',ylabel='Y (10⁶ km)',zlabel='Z (10⁶ km)',title='Saturn-centered J2000 trajectory');b.set(xlabel='ΔX (m)',ylabel='ΔY (m)',zlabel='ΔZ (m)',title='WLS − simulated truth');a.legend(fontsize=8);b.legend(fontsize=8)
    save(fig,'vlbi_trajectory_3d')
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    fig=make_subplots(rows=1,cols=2,specs=[[{'type':'scene'},{'type':'scene'}]],subplot_titles=['Saturn-centered J2000 trajectory','WLS − simulated truth (m)'])
    fig.add_trace(go.Scatter3d(x=truth[:,0]/1e6,y=truth[:,1]/1e6,z=truth[:,2]/1e6,mode='lines',name='Saved simulation',line=dict(color='#303c45')),row=1,col=1)
    for fit,label,color in zip(fits,labels,colors):
        for col,points in [(1,fit/1e6),(2,(fit-truth)*1000)]:fig.add_trace(go.Scatter3d(x=points[:,0],y=points[:,1],z=points[:,2],mode='lines',name=label,legendgroup=label,showlegend=col==1,line=dict(color=color)),row=1,col=col)
    fig.update_layout(height=700,title='Cassini eight-parameter WLS comparison',scene=dict(xaxis_title='X (10⁶ km)',yaxis_title='Y (10⁶ km)',zaxis_title='Z (10⁶ km)',aspectmode='data'),scene2=dict(xaxis_title='ΔX (m)',yaxis_title='ΔY (m)',zaxis_title='ΔZ (m)',aspectmode='data'))
    fig.write_html(root/'vlbi_trajectory_3d.html',include_plotlyjs=True)
    merged=orbits[0][['elapsed_s']+truthcols].copy()
    for prefix,data in zip(['no_vlbi','with_vlbi'],orbits):
        for col in ['fitted_x_km','fitted_y_km','fitted_z_km']+columns+['formal_position_sigma_m']:merged[prefix+'_'+col]=data[col]
    merged.to_csv(root/'comparison.csv',index=False,float_format='%.17g')
    metrics={}
    for label,s,d in zip(['no_vlbi','with_vlbi'],summaries,orbits):
        metrics[label]={key:s[key] for key in ['fitted_position_rms_m','fitted_endpoint_error_m','postfit_range_rms_m','postfit_doppler_rms_hz','iterations','weighted_postfit_rms']}
        metrics[label].update({axis+'_rms_m':float(np.sqrt(np.mean(d[col]**2))) for axis,col in zip(['radial','along','cross'],columns[1:])})
        metrics[label]['formal_position_rms_m']=float(np.sqrt(np.mean(d.formal_position_sigma_m**2)))
    metrics['position_rms_improvement_factor']=metrics['no_vlbi']['fitted_position_rms_m']/metrics['with_vlbi']['fitted_position_rms_m']
    assert set((vlbi.elapsed_s/86400).astype(int))==set(range(14))
    metrics['vlbi_first_day']=float(vlbi.elapsed_s.iloc[0]/86400);metrics['vlbi_last_day']=float(vlbi.elapsed_s.iloc[-1]/86400);metrics['vlbi_maximum_gap_hours']=float(vlbi.elapsed_s.diff().max()/3600);metrics['vlbi_minimum_elevation_deg']=float(vlbi[['elevation_one_deg','elevation_two_deg']].min().min())
    metrics['vlbi_count']=len(vlbi);metrics['vlbi_pair_counts']=summaries[1]['vlbi_pair_counts'];metrics['postfit_vlbi_rms_m']=summaries[1]['postfit_vlbi_rms_m'];metrics['identical_two_way_csv_sha256']=hashlib.sha256((folders[0]/'observations.csv').read_bytes()).hexdigest()
    (root/'comparison_summary.json').write_text(json.dumps(metrics,indent=2)+'\n')
    lines=['# Cassini WLS with and without VLBI','', f"Initial-guess position/velocity offset multiplier: {summaries[0].get('initial_guess_scale',1):g}. The two fits retain identical prior means and uncertainties.", '', 'Both fits use the saved simulation as truth over 9–23 October 2004. C++ generates all observations, estimates the same eight parameters [r₀, v₀, C_s, λ], and exports the numerical results. Python validates paired inputs and analyzes/plots those exports. Two-way data, noise samples, priors, initial guess, forces and spacecraft mass are identical.', '', '| Metric | Range + Doppler | Range + Doppler + 60 VLBI |','|---|---:|---:|']
    for name,key,unit in [('3D position RMS','fitted_position_rms_m','m'),('Endpoint error','fitted_endpoint_error_m','m'),('Radial RMS','radial_rms_m','m'),('Along-track RMS','along_rms_m','m'),('Cross-track RMS','cross_rms_m','m'),('Formal position RMS √mean(trace(Prr))','formal_position_rms_m','m'),('Range residual RMS','postfit_range_rms_m','m'),('Doppler residual RMS','postfit_doppler_rms_hz','Hz')]:lines.append(f"| {name} ({unit}) | {metrics['no_vlbi'][key]:.6g} | {metrics['with_vlbi'][key]:.6g} |")
    lines+=['',f"Position RMS improvement factor: {metrics['position_rms_improvement_factor']:.3f}. Added VLBI postfit residual RMS: {metrics['postfit_vlbi_rms_m']*1000:.3f} mm.",'', 'The 60 scalar measurements target evenly spaced epochs across the arc, selecting nearby jointly visible times without duplicate epochs. Both stations must exceed 10° elevation at their respective receptions; baseline preferences rotate among DSS-14/43, DSS-14/63 and DSS-43/63. Pair counts: '+str(metrics['vlbi_pair_counts'])+'.','', 'VLBI uses a common spacecraft emission event and solves the second station’s reception time, including station motion, solar Shapiro and troposphere. The observable is c times the signed second-minus-first coordinate TDB arrival delay. Differential path is evaluated directly to avoid subtraction of large absolute times and ranges.','', '**Synthetic assumption:** independent 1 cm differential-path noise (about 33 ps), known station clocks and matched propagation corrections. This is a near-field delay experiment; quasar calibration, real clock biases and media mismodeling are not included. The demonstrated improvement is conditional on these assumptions and this noise realization.','', f"The largest gains are along-track and cross-track. Radial RMS changes from {metrics['no_vlbi']['radial_rms_m']:.3f} to {metrics['with_vlbi']['radial_rms_m']:.3f} m; a joint noisy estimate need not improve every component in one realization. Thermal-scale formal uncertainty changes from {summaries[0]['initial_posterior_sigma'][7]:.6g} to {summaries[1]['initial_posterior_sigma'][7]:.6g}, while the SRP-scale uncertainty remains close to its 0.2 prior.", '', 'The original simulated force model and truth trajectory remain unchanged. Formal covariance comes from the final undamped eight-parameter information matrix; both fits must meet the same posterior-stationarity criterion. VLBI derivatives and half-step prediction differences are checked against 1% of the assumed noise.','', '[Without VLBI numerical results](no_vlbi/BASELINE.md) · [With VLBI numerical results](with_vlbi/BASELINE.md) · [Comparison CSV](comparison.csv) · [60 generated VLBI observations](with_vlbi/vlbi_observations.csv) · [Interactive 3D plot](vlbi_trajectory_3d.html)','', '![Position error comparison](vlbi_comparison.png)','', '![3D trajectory and position error](vlbi_trajectory_3d.png)','', '![Observation residuals](vlbi_residual_comparison.png)','', '![Force scales and uncertainty](vlbi_parameter_comparison.png)','']
    if summaries[0].get('pass_schedule',False):
        s=summaries[0]
        lines+=['', '## Pass-based experiment', '', f"{s['pass_count']} daily contiguous passes: {s['doppler_observations']} Doppler observations with 60 s counts/cadence and {s['range_observations']} range observations every 20 minutes. Range adds independent {s['range_independent_sigma_m']:g} m noise and a shared {s['range_pass_bias_sigma_m']:g} m offset per pass. The covariance R_pass=σ_ind²I+σ_bias²11ᵀ is used by C++ whitening, marginalizing offsets without expanding the eight-parameter state. The injected offsets are not revealed to WLS. The 10× larger starting error and original priors are retained. Doppler remains independent 0.01 Hz Gaussian noise; this is a controlled synthetic experiment rather than an actual DSN tracking-log reconstruction.", '', '[Pass schedule](with_vlbi/passes.csv) · [Pass diagnostics](with_vlbi/pass_range_analysis.csv)']
    (root/'COMPARISON.md').write_text('\n'.join(lines));print(json.dumps(metrics,indent=2))
if __name__=='__main__':main()
