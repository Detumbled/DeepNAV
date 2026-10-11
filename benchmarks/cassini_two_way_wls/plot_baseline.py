"""Validate and plot a completed Cassini coherent two-way batch WLS baseline."""
from pathlib import Path
import json,sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

folder=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path('results/cassini_two_way_wls').resolve()
s=json.loads((folder/'summary.json').read_text());orbit=pd.read_csv(folder/'trajectory.csv');r=pd.read_csv(folder/'residuals.csv');observations=pd.read_csv(folder/'observations.csv');it=pd.read_csv(folder/'iterations.csv');schedule=pd.read_csv(folder/'schedule.csv');cov=np.loadtxt(folder/'initial_covariance.csv',delimiter=',')
assert s['converged'] and len(r)==s['tracking_epochs']==len(observations)
assert s['scalar_observations']==2*len(r)
assert np.isfinite(orbit.to_numpy()).all() and np.isfinite(cov).all()
np.linalg.cholesky(cov);np.testing.assert_allclose(cov,cov.T,rtol=1e-10,atol=1e-15)
np.testing.assert_array_equal(r.elapsed_s,observations.elapsed_s)
assert (r.station==observations.station).all()
assert int(schedule.tracking.sum())==len(r)
for key,column in [('fitted_position_rms_m','fitted_position_error_m'),('nominal_position_rms_m','nominal_position_error_m'),('seed_position_rms_m','seed_position_error_m')]:np.testing.assert_allclose(s[key],np.sqrt(np.mean(orbit[column]**2)),rtol=1e-12)
np.testing.assert_allclose(s['postfit_range_rms_m'],np.sqrt(np.mean(r.postfit_range_m**2)),rtol=1e-12)
np.testing.assert_allclose(s['postfit_doppler_rms_hz'],np.sqrt(np.mean(r.postfit_doppler_hz**2)),rtol=1e-12)
np.testing.assert_allclose(s['initial_posterior_sigma'],np.sqrt(np.diag(cov)),rtol=1e-12)
initial_csv=pd.read_csv(folder/'initial_estimate.csv')
assert initial_csv.component.tolist()==['x','y','z','vx','vy','vz','C_s','lambda']
np.testing.assert_allclose(initial_csv.estimate,s['initial_estimated_state'],rtol=1e-12)
np.testing.assert_allclose(initial_csv.posterior_sigma,s['initial_posterior_sigma'],rtol=1e-12)
colors={'DSS-14':'#2878b5','DSS-43':'#dd8452','DSS-63':'#55a868'}
def save(fig,name):
    fig.savefig(folder/(name+'.png'),dpi=180);fig.savefig(folder/(name+'.pdf'));plt.close(fig)
def clean(axes):
    for ax in np.asarray(axes).ravel():ax.grid(alpha=.22)

days=orbit.elapsed_s/86400;rdays=r.elapsed_s/86400
fig,axes=plt.subplots(2,2,figsize=(13,8),layout='constrained');clean(axes)
axes[0,0].semilogy(days,orbit.seed_position_error_m,label='Biased initial guess',color='#a8a8a8')
axes[0,0].semilogy(days,np.maximum(orbit.nominal_position_error_m,1e-6),label='Reintegrated nominal model',color='#dd8452')
axes[0,0].semilogy(days,np.maximum(orbit.fitted_position_error_m,1e-6),label='Two-way WLS',color='#008f95')
axes[0,0].semilogy(days,orbit.formal_position_sigma_m,'--',label='Formal √trace(P position)',color='#8a63ad')
axes[0,0].set_ylabel('3D position error (m)');axes[0,0].legend(fontsize=9)
for label,key,color in [('Radial','radial_error_m','#dd8452'),('Along-track','along_error_m','#2878b5'),('Cross-track','cross_error_m','#8a63ad')]:axes[0,1].plot(days,orbit[key],label=label,color=color)
axes[0,1].set_ylabel('WLS − simulation RTN error (m)');axes[0,1].legend(fontsize=9)
for station,color in colors.items():
    mask=r.station==station;axes[1,0].scatter(rdays[mask],r.postfit_range_m[mask],s=4,color=color,label=station,alpha=.65);axes[1,1].scatter(rdays[mask],r.postfit_doppler_hz[mask],s=4,color=color,alpha=.65)
axes[1,0].set_ylabel('Postfit range residual (m)');axes[1,1].set_ylabel('Postfit counted Doppler residual (Hz)');axes[1,0].legend(fontsize=9)
for ax in axes.ravel():ax.set_xlabel('Elapsed TDB days from initial state')
fig.suptitle(f"Cassini · {s['duration_days']:g}-day coherent two-way WLS · {s['mass_kg']:g} kg\nEight parameters: six initial states + SRP and recoil scales; truth={s['truth']}")
save(fig,'diagnostic_overview')

fig,axes=plt.subplots(2,2,figsize=(13,8),layout='constrained');clean(axes)
for j,(name,unit,sigma) in enumerate([('range','m',s['range_sigma_m']),('doppler','hz',s['doppler_sigma_hz'])]):
    for phase,color in [('prefit','#a8a8a8'),('postfit','#008f95')]:axes[j,0].plot(rdays,r[f'{phase}_{name}_{unit}']/sigma,'.',ms=2,label=phase,color=color,alpha=.65)
    axes[j,0].set_ylabel(f'{name.title()} residual / assumed σ');axes[j,0].legend();axes[j,0].set_xlabel('Elapsed TDB days')
    axes[j,1].hist(r[f'postfit_{name}_{unit}']/sigma,bins=45,density=True,color='#008f95',alpha=.7)
    x=np.linspace(-4,4,200);axes[j,1].plot(x,np.exp(-x*x/2)/np.sqrt(2*np.pi),'k--',label='N(0,1) reference');axes[j,1].set_xlabel('Postfit residual / assumed σ');axes[j,1].set_ylabel('Density');axes[j,1].legend()
fig.suptitle('Synthetic two-way measurement residuals · observation minus prediction')
save(fig,'residuals')

fig,axes=plt.subplots(1,2,figsize=(12,4),layout='constrained');clean(axes)
axes[0].plot(it.iteration,it.weighted_postfit_rms,'o-',color='#008f95');axes[0].set_xlabel('WLS iteration');axes[0].set_ylabel('Nonlinear weighted residual RMS')
axes[1].semilogy(it.iteration,np.maximum(it.normalized_correction_norm,1e-12),'o-',color='#2878b5');axes[1].set_xlabel('WLS iteration');axes[1].set_ylabel('Prior-scaled correction norm')
fig.suptitle('Batch WLS convergence · same initial state and same observation set at every iteration')
save(fig,'convergence')

fig=plt.figure(figsize=(12,6),layout='constrained');a=fig.add_subplot(121,projection='3d');b=fig.add_subplot(122,projection='3d')
truth_prefix='truth'
truth_label='Saved simulation'
truth=orbit[[f'{truth_prefix}_{axis}_km' for axis in 'xyz']].to_numpy();fit=orbit[['fitted_x_km','fitted_y_km','fitted_z_km']].to_numpy()
for points,label,color in [(truth,truth_label,'#2d3947'),(fit,'WLS','#008f95')]:a.plot(*(points/1e6).T,label=label,color=color)
a.set(xlabel='X (10⁶ km)',ylabel='Y (10⁶ km)',zlabel='Z (10⁶ km)',title='Saturn-centered J2000 trajectory');a.legend()
delta=(fit-truth)*1000;b.plot(*delta.T,color='#008f95');b.scatter(*delta[-1],color='#dd8452');b.set(xlabel='ΔX (m)',ylabel='ΔY (m)',zlabel='ΔZ (m)',title=f'WLS minus {truth_label} position')
save(fig,'trajectory_3d')
try:
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    figure=make_subplots(rows=1,cols=2,specs=[[{'type':'scene'},{'type':'scene'}]],subplot_titles=['Saturn-centered trajectory',f'WLS − {truth_label} position'])
    for points,label,color in [(truth/1e6,truth_label,'#2d3947'),(fit/1e6,'WLS','#008f95')]:figure.add_trace(go.Scatter3d(x=points[:,0],y=points[:,1],z=points[:,2],mode='lines',name=label,line=dict(color=color)),row=1,col=1)
    figure.add_trace(go.Scatter3d(x=delta[:,0],y=delta[:,1],z=delta[:,2],mode='lines',name='Position residual (m)',line=dict(color='#008f95')),row=1,col=2)
    figure.update_layout(title='Cassini · coherent two-way WLS baseline',scene=dict(xaxis_title='X (10⁶ km)',yaxis_title='Y (10⁶ km)',zaxis_title='Z (10⁶ km)',aspectmode='data'),scene2=dict(xaxis_title='ΔX (m)',yaxis_title='ΔY (m)',zaxis_title='ΔZ (m)',aspectmode='data'),height=700)
    figure.write_html(folder/'trajectory_3d.html',include_plotlyjs=True)
except ImportError:pass

# Parameter uncertainty and correlation are part of the eight-state baseline.
correlation=cov/np.sqrt(np.outer(np.diag(cov),np.diag(cov)))
fig,axes=plt.subplots(1,2,figsize=(12,5),layout='constrained')
for j,label in enumerate(['C_s (SRP scale)','λ (thermal scale)']):
    axes[0].errorbar(j-.08,s['initial_prior_mean'][6+j],yerr=s['initial_prior_sigma'][6+j],fmt='o',color='#a8a8a8',label='Prior 1σ' if j==0 else None)
    axes[0].errorbar(j+.08,s['initial_estimated_state'][6+j],yerr=s['initial_posterior_sigma'][6+j],fmt='o',color='#008f95',label='Posterior 1σ' if j==0 else None)
axes[0].axhline(1,color='black',ls='--',label='Truth');axes[0].set_xticks([0,1],['C_s','λ']);axes[0].set_ylabel('Dimensionless force scale');axes[0].legend();axes[0].grid(alpha=.2)
labels=['X','Y','Z','VX','VY','VZ','C_s','λ']
mesh=axes[1].imshow(correlation,vmin=-1,vmax=1,cmap='RdBu_r');axes[1].set_xticks(range(8),labels);axes[1].set_yticks(range(8),labels);axes[1].set_title('Posterior parameter correlation');fig.colorbar(mesh,ax=axes[1])
save(fig,'parameters')
lines=['# Cassini simulated-truth two-way WLS baseline','',
    f"**{s['duration_days']:g} TDB days from {s['start_utc']}; {s['mass_kg']:g} kg; eight estimated parameters.**",'',
    '## Ground truth','',
    'Observations are generated in C++ from the saved `multibody_j2_srp_thermal` rows of `results/cassini_trajectory/trajectory.csv`. Cubic Hermite interpolation uses saved positions and velocities, with acceleration from the retained model for velocity interpolation. Spacecraft SPK states are not used as observation truth or as the accuracy reference. Planetary, station and attitude kernels remain inputs to the dynamical and radio models.', '',
    '| Result | Value |','|---|---:|',
    f"| Tracking epochs / scalar observations | {s['tracking_epochs']} / {s['scalar_observations']} |",
    f"| Iterations / converged | {s['iterations']} / {s['converged']} |",
    f"| WLS 3D position RMS against simulation | {s['fitted_position_rms_m']:.3f} m |",
    f"| Endpoint position error | {s['fitted_endpoint_error_m']:.3f} m |",
    f"| Prefit → postfit range RMS | {s['prefit_range_rms_m']:.3f} → {s['postfit_range_rms_m']:.3f} m |",
    f"| Prefit → postfit Doppler RMS | {s['prefit_doppler_rms_hz']:.6f} → {s['postfit_doppler_rms_hz']:.6f} Hz |",
    f"| Weighted postfit RMS / reduced χ² (N−8) | {s['weighted_postfit_rms']:.6f} / {s['reduced_chi_square']:.6f} |",'',
    '## Augmented estimator','',
    'The estimated vector is [X₀, Y₀, Z₀, VX₀, VY₀, VZ₀, C_s, λ]. The scales are constant over the arc: a_SRP = C_s a_SRP,nominal and a_thermal = λ a_thermal,nominal. Both truth scales are 1. C_s scales the entire nominal attitude-dependent SRP acceleration, rather than changing the reflection law through Cr. λ scales the retained signed XYZ recoil vector, preserving its attitude dependence and radioactive decay.', '',
    'Nominal forces are multibody tides, direct Saturn J2, uniform rings, attitude-dependent HGA SRP and thermal recoil. Gravity mask 2 retains the previously selected model. Mass, area, Cr and all remaining coefficients stay fixed.', '',
    'Initial position/velocity offsets are [1, −0.8, 0.5] km and [0.01, −0.008, 0.006] m/s. Initial force scales are 1. Gaussian priors have 2 km position and 0.1 m/s velocity standard deviations, centered on the perturbed initial state. Force-scale priors are centered at 1, with assumed 20% SRP and 50% thermal uncertainties. These illustrative priors are part of the solution, not uncertainty estimates from a Cassini calibration.', '',
    'C++ integrates the eight-state vector and 8×8 variational matrix. The two parameter sensitivity columns are driven by the nominal SRP and recoil accelerations. Full counted-link central differences use the propagated sensitivities; accepted corrections use fresh nonlinear propagation. The existing WLS and batch driver minimize weighted radio residuals plus the prior penalty. No Cartesian truth samples are fitted, and no process noise, state resets or maneuver impulses are estimated.', '',
    '| Parameter | Estimate | Truth | Posterior 1σ | Prior 1σ |','|---|---:|---:|---:|---:|']
for j,name in enumerate(['X₀ (km)','Y₀ (km)','Z₀ (km)','VX₀ (km/s)','VY₀ (km/s)','VZ₀ (km/s)','C_s','λ']):
    lines.append(f"| {name} | {s['initial_estimated_state'][j]:.12g} | {s['initial_truth_state'][j]:.12g} | {s['initial_posterior_sigma'][j]:.6g} | {s['initial_prior_sigma'][j]:.6g} |")
lines+=['', 'Uncertainties are conditional on the remaining fixed model and assumed independent noise. Coefficients close to 1 do not demonstrate recovery if posterior uncertainties remain close to their priors. The force-scale uncertainty plot and full correlation matrix expose this limitation. This is a retrospective batch estimate rather than a sequential filter.', '',
    '## Radio setup','',
    f"DSS-14/43/63 counts: {s['station_counts']}. One station per epoch, highest minimum elevation over both legs at both count endpoints, cutoff 10°. Cadence {s['cadence_s']:g} s; end-tagged count {s['count_s']:g} s. Reception begins 6000 s after the state epoch so retarded spacecraft events lie inside the arc. This is an idealized visibility schedule rather than a flight tracking schedule.", '',
    f"Independent Gaussian noise: {s['range_sigma_m']:g} m range and {s['doppler_sigma_hz']:g} Hz Doppler; seed {s['noise_seed']}; noise enabled {s['noise_enabled']}. Range is c·RTLT/2. Coherent X/X ratio 880/749, 7.175 GHz constant uplink, known 1 µs delay, received-minus-reference Doppler. Matching solar Shapiro, 2.3 m zenith troposphere and smooth SPICE TDB−TT ground clock are used for generation and prediction. No independent spacecraft oscillator, dispersive plasma, observation biases or clock parameters are estimated.", '',
    '## Verification','',
    'Force position derivatives and the augmented STM are checked against finite differences/nonlinear DP853 propagation. Coefficient radio columns use large affine sensitivity stencils to resolve the small signals above the counted-Doppler numerical floor; this does not propagate a physical trajectory with those large coefficients. A half-stencil check and a final half-step nonlinear integration check are performed.',
    f"Final undamped posterior correction: {s['undamped_posterior_step_sigma']:.6g} σ, required below 0.1 σ. A trial step below 0.1 posterior σ allows 0.1% objective-resolution tolerance for counted-Doppler arithmetic; larger steps must decrease the objective.",
    f"Radio half-stencil relative difference: {s['radio_stencil_relative_difference']:.6g}; discrepancy over posterior uncertainty is {s['radio_stencil_posterior_noise_fraction']:.6g} of measurement noise (required <0.01). Weighted design condition number: {s['measurement_design_condition_number']:.6g}.",
    f"Half-step checks: position {s['nominal_refined_maximum_m']*1000:.3f} mm, range {s['nominal_refined_range_maximum_m']*1000:.3f} mm, Doppler {s['nominal_refined_doppler_maximum_hz']:.6g} Hz.",
    f"Reintegrated nominal model differs from the saved truth by at most {s.get('retained_orbit_reproduction_maximum_m',0)*1000:.3f} mm at sampled epochs.", '',
    f"Radio differentiation offsets (km, km/s, dimensionless scales): {s['radio_differentiation_steps']}.", '',
    'All observation generation, propagation, estimation, covariance and numerical CSV exports are C++. Python only validates/analyzes those files and generates plots. The summary stores source, kernel, input and observation hashes.', '',
    '## Outputs','',
    '[Generated observations](observations.csv) · [Initial eight-parameter estimate](initial_estimate.csv) · [Estimated trajectory](trajectory.csv) · [Residuals](residuals.csv) · [8×8 covariance](initial_covariance.csv)', '',
    '[Rotatable 3D trajectory](trajectory_3d.html)', '',
    '![Trajectory and radio residuals](diagnostic_overview.png)', '',
    '![Force scales and correlations](parameters.png)', '',
    '![Radio residual distributions](residuals.png)', '',
    '![WLS convergence](convergence.png)', '']
(folder/'BASELINE.md').write_text('\n'.join(lines))
print(json.dumps({key:s[key] for key in ['converged','iterations','tracking_epochs','fitted_position_rms_m','postfit_range_rms_m','postfit_doppler_rms_hz','nominal_refined_maximum_m']},indent=2))
