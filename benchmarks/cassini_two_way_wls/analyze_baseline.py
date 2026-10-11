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
    parser.add_argument('--output', type=Path, default=ROOT/'results/cassini_two_way_wls')
    args = parser.parse_args()
    source = ROOT/'results/cassini_trajectory'
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
        np.testing.assert_array_equal(old.elapsed_s,trajectory.elapsed_s)
        np.testing.assert_allclose(old[['x_km','y_km','z_km']].to_numpy(),trajectory[['truth_x_km','truth_y_km','truth_z_km']].to_numpy(),rtol=0,atol=1e-8)
        difference=np.linalg.norm(old[['x_km','y_km','z_km']].to_numpy()-trajectory[['nominal_x_km','nominal_y_km','nominal_z_km']].to_numpy(),axis=1)*1000
        summary['retained_orbit_reproduction_maximum_m']=float(difference.max())
        assert difference.max()<.1, 'WLS force model does not reproduce retained orbit'
    source_paths=[ROOT/name for name in ['benchmarks/cassini_two_way_wls/baseline.cpp','benchmarks/cassini_trajectory/CassiniForceModel.hpp','benchmarks/cassini_trajectory/CassiniAttitude.hpp','benchmarks/cassini_trajectory/Sat427Gravity.hpp','include/RKF45Integrator.hpp','src/RKF45Integrator.cpp','src/dynamics/EphemerisInterpolator.cpp','src/observations/radiometric/TwoWayLink.cpp','src/observations/radiometric/TwoWayMeasurement.cpp','src/filters/WLS.cpp','src/filters/BatchLeastSquaresDriver.cpp','include/perturbations/Gravitational.hpp','include/perturbations/J2.hpp','src/perturbations/J2.cpp','include/perturbations/AnnularGravity.hpp','include/perturbations/FlatPlateSRP.hpp','include/perturbations/ThermalRecoil.hpp']]
    summary['source_sha256']={str(path.relative_to(ROOT)):digest(path) for path in source_paths}
    summary['source_arc_metadata']=metadata
    summary['source_recoil_study']=source_summary.get('recoil_study')
    summary['input_sha256']={str(path.relative_to(ROOT)):digest(path) for path in [source/'metadata.json',source/'trajectory.csv',source/'summary.json']}
    radio_kernels=['earthstns_itrf93_201023.bsp','earth_1962_250826_2125_combined.bpc','earth_assoc_itrf93.tf','Cassini/fk/earth_topo_201023.tf']
    kernel_paths=[Path(path) for path in metadata['kernel_paths']]+[ROOT/'Kernels'/name for name in radio_kernels]
    summary['kernel_sha256']={str(path.relative_to(ROOT)):digest(path) for path in kernel_paths}
    summary['observation_sha256']=digest(output/'observations.csv')
    summary['setup']={'elevation_cutoff_deg':10,'reception_start_margin_s':6000,'station_names':['DSS-14','DSS-43','DSS-63'],'station_selection':'Highest minimum elevation across uplink/downlink at both count endpoints; one station per epoch','range_convention':'c * round-trip light time / 2, km','doppler_convention':'End-tagged counted received-minus-reference frequency, Hz','turnaround_ratio':'880/749','uplink_hz':7.175e9,'transponder_delay_s':1e-6,'zenith_troposphere_m':2.3,'solar_shapiro':True,'ground_clock':'SPICE smooth TDB-TT','spacecraft_clock':'Ideal coherent transponder; no independent spacecraft oscillator state','fit':'Six initial Saturn-relative Cartesian components plus Cs and lambda; Gaussian priors; fixed remaining force coefficients; no state resets or process noise','covariance':'Formal eight-parameter covariance conditional on remaining fixed forces and assumed independent observation errors'}
    summary_path.write_text(json.dumps(summary,indent=2)+'\n')
    sys.argv=[str(Path(__file__).with_name('plot_baseline.py')),str(output)]
    runpy.run_path(sys.argv[0],run_name='__main__')

if __name__=='__main__':main()
