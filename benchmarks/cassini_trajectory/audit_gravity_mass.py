"""Check SAT427 pole against CSPICE and verify October mass from NASA PDS."""
import csv
import hashlib
import json
import re
import subprocess
import tempfile
from pathlib import Path
import numpy as np
import spiceypy as spice

ROOT = Path(__file__).resolve().parents[2]
DATA = Path(__file__).resolve().parent / 'data'
LABEL = ROOT / 'Kernels/Cassini/spk/200128RU_SCPSE_04251_04303.bsp.lbl'
BASE = 'https://pds-atmospheres.nmsu.edu/data_and_services/atmospheres_data/Cassini/files/'

def main():
    text = LABEL.read_text()
    constants = {k: float(v) for k, v in re.findall(r'\b([A-Z][A-Za-z0-9_]*)\s+([+-]?\d+\.\d+E[+-]\d+)', text)}
    masses = {name.upper(): float(value) for name, value in re.findall(
        r'\b(Mimas|Enceladus|Tethys|Dione|Rhea|Titan|Hyperion|Iapetus|Phoebe|Helene|Saturn)\s+\d+\s+(\d+\.\d+E[+-]\d+)', text)}
    assert constants['POLTIM'] == 2451545.0
    ring_gm=sum(constants[f'RINGM{i}'] for i in (1,2,4))
    np.testing.assert_allclose(constants['System'],sum(masses.values())+ring_gm,rtol=0,atol=1e-6)
    # CSPICE independently evaluates the same text-PCK convention. Include ±100
    # years to catch days/years/centuries mistakes, not only the current arc.
    for name in ['naif0012.tls', 'pck00011.tpc', 'gm_de440.tpc']:
        spice.furnsh(str(ROOT / 'Kernels' / name))
    spice.pdpool('BODY699_POLE_RA', [constants['ZACPL6'], constants['DACPL6'], 0])
    spice.pdpool('BODY699_POLE_DEC', [constants['ZDEPL6'], constants['DDEPL6'], 0])
    spice.pdpool('BODY699_NUT_PREC_RA', [constants[f'RA_AMPL{i}'] for i in range(1,4)])
    spice.pdpool('BODY699_NUT_PREC_DEC', [constants[f'DE_AMPL{i}'] for i in range(1,4)])
    spice.pdpool('BODY699_NUT_PREC_PM', [0, 0, 0])
    assert all(constants[f'RA_FAZE{i}'] == constants[f'DE_FAZE{i}'] and
               constants[f'RA_FREQ{i}'] == constants[f'DE_FREQ{i}'] for i in range(1,4))
    spice.pdpool('BODY6_NUT_PREC_ANGLES', [constants[f'{key}{i}'] for i in range(1,4) for key in ('RA_FAZE','RA_FREQ')])
    with tempfile.TemporaryDirectory(prefix='deepnav-gravity-audit-') as tmp:
        subprocess.run([str(ROOT/'build-clang/compare_cassini_spk'), '--gravity-audit', '--output', tmp], check=True)
        rows = list(csv.DictReader((Path(tmp)/'gravity_audit.csv').open()))
    max_difference = 0
    for row in rows:
        pole = np.array([float(row[key]) for key in ('px','py','pz')])
        reference = spice.pxform('IAU_SATURN','J2000',float(row['epoch_tdb']))[:,2]
        np.testing.assert_allclose(pole, reference, rtol=0, atol=5e-15)
        np.testing.assert_allclose(np.linalg.norm(pole), 1, rtol=0, atol=3e-16)
        max_difference=max(max_difference,float(np.linalg.norm(pole-reference)))
    mass_rows=list(csv.DictReader((DATA/'spacecraft_mass.csv').open()))
    relevant=[row for row in mass_rows if row['Maneuver'] in ('OTM-003','OTM-004','OTM-005')]
    assert len(relevant)==3
    assert all(float(row['Spacecraft Mass (kg)']) == 3220 for row in relevant)
    assert '"pre-burn"' in (DATA/'spacecraft_mass_README.txt').read_text()
    manifest=dict(gravity_release='SAT427.23',
        gravity_label_url='https://naif.jpl.nasa.gov/pub/naif/CASSINI/kernels/spk/'+LABEL.name,
        gravity_label_sha256=hashlib.sha256(LABEL.read_bytes()).hexdigest(),
        gravity_constants=constants, body_gm_km3_s2=masses,
        ring_total_gm_km3_s2=ring_gm,
        ring_mass_check='System GM equals sum of Saturn and modeled nonzero moons plus RINGM1/2/4 (absolute tolerance 1e-6 km3/s2).',
        pole_formula='RA=RA0+rate*T+sum(A*sin(phase+frequency*T)); DEC=DEC0+rate*T+sum(A*cos(phase+frequency*T)); degrees; T=TDB Julian centuries since POLTIM',
        pole_reference='https://doi.org/10.3847/1538-3881/ac90c9',
        max_cspice_pole_difference=max_difference,
        solar_gm_convention='SAT427 GM10 includes inner planets; subtract explicit barycenter GMs 1–4 from gm_de440.tpc to avoid double counting. This is a resolved approximation, not SAT427’s lumped solar model.',
        mass_kg=3220, arc_start_utc='2004-10-09T00:00:00', arc_end_utc='2004-10-23T00:00:00',
        mass_policy='Constant estimate from rounded pre-burn OTM-004 mass; the table provides maneuver snapshots, not a continuous mass history. No interpolation across burns. Huygens still attached.',
        mass_rows=relevant,
        mass_sources=[dict(path=str((DATA/name).relative_to(ROOT)), url=BASE+url,
            sha256=hashlib.sha256((DATA/name).read_bytes()).hexdigest()) for name,url in (
            ('spacecraft_mass.csv','Cassini%20Spacecraft%20Mass%20Table.csv'),
            ('spacecraft_mass_README.txt','README%20for%20Cassini%20Spacecraft%20Mass%20table.txt'))],
        maneuver_time_note='Mass table lists OTM-004 2004-10-23 06:16 UTC; weekly NASA report says 07:29:28 UTC. Both are after the arc; benchmark uses earlier time conservatively.')
    (DATA/'gravity_mass_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(f'PASS: six pole epochs match CSPICE (max vector difference {max_difference:.3g}); NASA PDS October mass = 3220 kg.')

if __name__=='__main__': main()
