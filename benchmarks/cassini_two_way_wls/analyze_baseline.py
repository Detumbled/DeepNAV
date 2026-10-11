"""Analyze completed C++ coherent two-way WLS CSV outputs and plot results."""
from pathlib import Path
import argparse, hashlib, json, runpy, sys
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]

def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=ROOT/'results/cassini_two_way_wls/pass_schedule/no_vlbi')
    args = parser.parse_args()
    source = ROOT/'results/cassini_two_way_wls/pass_schedule/truth'
    metadata = json.loads((source/'metadata.json').read_text())
    source_summary = json.loads((source/'summary.json').read_text())
    assert metadata['gravity_mask'] == 2 and metadata['numerical_validation']
    assert metadata['thermal_reference_mass_kg'] == 2781 and metadata['thermal_reference_utc'] == '2005-08-13T00:00:00'
    assert metadata['thermal_half_life_years'] == 87.74
    output=args.output.resolve()
    summary_path=output/'summary.json';summary=json.loads(summary_path.read_text())
    assert summary['converged'] and summary['estimated_parameters']==8 and summary['truth']=='simulation'
    assert summary['mass_kg']==metadata['mass_kg'] and summary['gravity_mask']==metadata['gravity_mask']
    assert summary['area_m2']==metadata['area_m2'] and summary['cr']==metadata['cr']
    np.testing.assert_array_equal(summary['thermal_reference_vector_km_s2'],[metadata['thermal_reference_'+axis+'_km_s2'] for axis in 'xyz'])
    trajectory=pd.read_csv(output/'trajectory.csv')
    assert np.isfinite(trajectory.to_numpy()).all()
    if summary['duration_days']==metadata['duration_days'] and summary['start_utc']==metadata['start_utc']:
        old=pd.read_csv(source/'trajectory.csv');old=old[old.model=='multibody_j2_srp_thermal'].reset_index(drop=True)
        sampled=trajectory[trajectory.elapsed_s.isin(old.elapsed_s)].reset_index(drop=True)
        np.testing.assert_array_equal(old.elapsed_s,sampled.elapsed_s)
        np.testing.assert_allclose(old[['x_km','y_km','z_km']].to_numpy(),sampled[['truth_x_km','truth_y_km','truth_z_km']].to_numpy(),rtol=0,atol=1e-8)
        difference=np.linalg.norm(old[['x_km','y_km','z_km']].to_numpy()-sampled[['nominal_x_km','nominal_y_km','nominal_z_km']].to_numpy(),axis=1)*1000
        assert trajectory.nominal_position_error_m.max()<.1, 'Interpolated simulation differs from reintegrated truth'
        summary['retained_orbit_reproduction_maximum_m']=float(difference.max())
        assert difference.max()<.1, 'WLS force model does not reproduce retained orbit'
    source_paths=[ROOT/name for name in ['benchmarks/cassini_two_way_wls/baseline.cpp','benchmarks/cassini_trajectory/SavedSimulation.hpp','benchmarks/cassini_trajectory/CassiniForceModel.hpp','benchmarks/cassini_trajectory/CassiniAttitude.hpp','benchmarks/cassini_trajectory/Sat427Gravity.hpp','include/RKF45Integrator.hpp','src/RKF45Integrator.cpp','src/dynamics/EphemerisInterpolator.cpp','src/observations/radiometric/TwoWayLink.cpp','src/observations/radiometric/TwoWayMeasurement.cpp','src/observations/radiometric/VlbiLink.cpp','src/observations/radiometric/VlbiMeasurement.cpp','include/observations/radiometric/VlbiLink.hpp','include/observations/radiometric/VlbiMeasurement.hpp','include/observations/radiometric/PassRangeWhitening.hpp','src/filters/WLS.cpp','src/filters/BatchLeastSquaresDriver.cpp','include/perturbations/Gravitational.hpp','include/perturbations/J2.hpp','src/perturbations/J2.cpp','include/perturbations/AnnularGravity.hpp','include/perturbations/FlatPlateSRP.hpp','include/perturbations/ThermalRecoil.hpp']]
    summary['source_sha256']={str(path.relative_to(ROOT)):digest(path) for path in source_paths}
    summary['source_arc_metadata']=metadata
    summary['source_recoil_study']=source_summary.get('recoil_study')
    summary['input_sha256']={str(path.relative_to(ROOT)):digest(path) for path in [source/'metadata.json',source/'trajectory.csv',source/'summary.json']}
    radio_kernels=['earthstns_itrf93_201023.bsp','earth_1962_250826_2125_combined.bpc','earth_assoc_itrf93.tf','Cassini/fk/earth_topo_201023.tf']
    kernel_paths=[Path(path) for path in metadata['kernel_paths']]+[ROOT/'Kernels'/name for name in radio_kernels]
    summary['kernel_sha256']={str(path.relative_to(ROOT)):digest(path) for path in kernel_paths}
    summary['observation_sha256']=digest(output/'observations.csv')
    summary['setup']={'elevation_cutoff_deg':10,'reception_start_margin_s':6000,'station_names':['DSS-14','DSS-43','DSS-63'],'station_selection':'Highest minimum elevation across uplink/downlink at both count endpoints; one station per epoch','range_convention':'c * round-trip light time / 2, km','doppler_convention':'End-tagged counted received-minus-reference frequency, Hz','turnaround_ratio':'880/749','uplink_hz':7.175e9,'transponder_delay_s':1e-6,'zenith_troposphere_m':2.3,'solar_shapiro':True,'ground_clock':'SPICE smooth TDB-TT','spacecraft_clock':'Ideal coherent transponder; no independent spacecraft oscillator state','fit':'Six initial Saturn-relative Cartesian components plus Cs and lambda; Gaussian priors; fixed remaining force coefficients; no state resets or process noise','covariance':'Formal eight-parameter covariance conditional on remaining fixed forces and assumed independent observation errors'}
    if summary.get('vlbi_observations',0):
        vlbi=pd.read_csv(output/'vlbi_observations.csv');residual=pd.read_csv(output/'vlbi_residuals.csv')
        assert len(vlbi)==len(residual)==summary['vlbi_observations']
        assert (vlbi[['elevation_one_deg','elevation_two_deg']]>=10).all().all()
        assert vlbi.elapsed_s.is_monotonic_increasing and vlbi.elapsed_s.is_unique
        assert (vlbi.elapsed_s>=6000).all() and (vlbi.elapsed_s<=summary['duration_days']*86400).all()
        np.testing.assert_array_equal(vlbi.elapsed_s,residual.elapsed_s)
        np.testing.assert_allclose(vlbi.observed_delay_km-vlbi.truth_delay_km,vlbi.noise_m/1000,rtol=0,atol=2e-12)
        np.testing.assert_allclose(summary['postfit_vlbi_rms_m'],np.sqrt(np.mean(residual.postfit_m**2)),rtol=1e-12)
        assert summary['refined_vlbi_maximum_m']<.01*summary['vlbi_sigma_m']
        assert summary['vlbi_stencil_posterior_noise_fraction']<.01
        summary['vlbi_observation_sha256']=digest(output/'vlbi_observations.csv')
        summary['setup']['vlbi']='Near-field common-emission differential TDB reception delay times c; known station clocks; solar Shapiro and troposphere; independent synthetic 1 cm path noise; not quasar-calibrated delta-DOR'
    if summary.get('pass_schedule',False):
        passes=pd.read_csv(output/'passes.csv');obs=pd.read_csv(output/'observations.csv',float_precision='round_trip');res=pd.read_csv(output/'residuals.csv')
        assert len(passes)==summary['pass_count']==int(summary['duration_days'])
        assert obs.range_active.sum()==summary['range_observations'] and len(obs)==summary['doppler_observations']
        assert summary['scalar_observations']==len(obs)+int(obs.range_active.sum())+summary['vlbi_observations']
        assert (passes.minimum_elevation_deg>=10).all() and (passes.usable_hours<=summary['pass_target_hours']).all()
        assert (passes.end_elapsed_s>passes.start_elapsed_s).all()
        assert (passes.start_elapsed_s.iloc[1:].to_numpy()>=passes.end_elapsed_s.iloc[:-1].to_numpy()).all()
        midpoint=(passes.start_elapsed_s+passes.end_elapsed_s)/2
        assert ((midpoint>=passes.pass_id*86400)&(midpoint<(passes.pass_id+1)*86400)).all()
        np.testing.assert_allclose(passes.usable_hours,(passes.end_elapsed_s-passes.start_elapsed_s)/3600,rtol=0,atol=1e-12)
        np.testing.assert_array_equal(obs.elapsed_s,res.elapsed_s)
        chi=float(np.sum((res.postfit_doppler_hz/summary['doppler_sigma_hz'])**2))
        means=[]
        for pass_id,d in obs.groupby('pass_id'):
            p=passes[passes.pass_id==pass_id].iloc[0];assert len(d)==p.doppler_points and int(d.range_active.sum())==p.range_points
            assert d.station.nunique()==1 and d.station.iloc[0]==p.station
            np.testing.assert_allclose(np.diff(d.elapsed_s),60,rtol=0,atol=0)
            ranges=d[d.range_active==1];np.testing.assert_allclose(np.diff(ranges.elapsed_s),summary['range_cadence_s'],rtol=0,atol=0)
            np.testing.assert_allclose(ranges.range_noise_m,ranges.range_independent_noise_m+ranges.range_bias_m,rtol=0,atol=1e-12)
            np.testing.assert_allclose(ranges.range_bias_m,p.true_range_bias_m,rtol=0,atol=1e-12)
            np.testing.assert_allclose(ranges.observed_range_km-ranges.truth_range_km,ranges.range_noise_m/1000,rtol=0,atol=5e-7)
            assert d[d.range_active==0].observed_range_km.isna().all()
            r=res.loc[ranges.index,'postfit_range_m'].to_numpy();n=len(r)
            sigma=summary['range_independent_sigma_m'];bias=summary['range_pass_bias_sigma_m']
            w=r/sigma;w=w+(sigma/np.sqrt(sigma*sigma+n*bias*bias)-1)*w.mean();chi+=float(w@w)
            means.append(dict(pass_id=int(pass_id),station=p.station,range_points=n,injected_bias_m=float(p.true_range_bias_m),postfit_range_mean_m=float(r.mean()),within_pass_range_rms_m=float(np.sqrt(np.mean((r-r.mean())**2))),pass_mean_sigma_m=float(np.sqrt(bias*bias+sigma*sigma/n))))
        if summary['vlbi_observations']:
            vr=pd.read_csv(output/'vlbi_residuals.csv');chi+=float(np.sum((vr.postfit_m/summary['vlbi_sigma_m'])**2))
        np.testing.assert_allclose(summary['weighted_postfit_rms'],np.sqrt(chi/summary['scalar_observations']),rtol=1e-9)
        pd.DataFrame(means).to_csv(output/'pass_range_analysis.csv',index=False)
        summary['setup']['station_selection']='One contiguous visibility-qualified pass per day, aiming for six hours; rotating station preference'
        summary['setup']['range_covariance']='Per pass: sigma_independent^2 I + sigma_shared_bias^2 11^T; symmetric whitening; Gaussian biases marginalized; eight estimated parameters'
        summary['setup']['covariance']='Formal eight-parameter covariance conditional on fixed forces, correlated range blocks and independent Doppler/VLBI errors'
        summary['setup']['range_noise']='Independent 1 m samples plus independent 10 m shared Gaussian offset per pass (synthetic assumptions)'
        summary['setup']['doppler_noise']='Independent 0.01 Hz Gaussian samples; media and antenna errors beyond the shared range offset are not simulated'
        summary['pass_csv_sha256']=digest(output/'passes.csv')
    summary_path.write_text(json.dumps(summary,indent=2)+'\n')
    sys.argv=[str(Path(__file__).with_name('plot_baseline.py')),str(output)]
    runpy.run_path(sys.argv[0],run_name='__main__')

if __name__=='__main__':main()
