KPL/MK

Cassini ISS NAC geometry for N1476124698_2.IMG and N1475349073_2.IMG.
Load from the DeepNAV repository root. Paths are relative to the process
working directory, not to this meta-kernel's directory.

The bundled reconstruction is loaded after DE442 so its corrected DE437
planetary states and SAT427 satellite states take precedence on its arc.

\begindata

PATH_VALUES  = ( 'Kernels/Cassini' )
PATH_SYMBOLS = ( 'CAS' )

KERNELS_TO_LOAD = (
    '$CAS/common/naif0012.tls'
    '$CAS/common/pck00011.tpc'
    '$CAS/common/gm_de440.tpc'
    '$CAS/common/de442.bsp'
    '$CAS/common/earthstns_itrf93_201023.bsp'
    '$CAS/common/earth_1962_250826_2125_combined.bpc'
    '$CAS/common/earth_assoc_itrf93.tf'
    '$CAS/fk/earth_topo_201023.tf'
    '$CAS/spk/200128RU_SCPSE_04251_04303.bsp'
    '$CAS/sclk/cas00172.tsc'
    '$CAS/fk/cas_v40.tf'
    '$CAS/ik/cas_iss_v10.ti'
    '$CAS/ck/04271_04276ra.bc'
    '$CAS/ck/04281_04286ra.bc'
)

\begintext
