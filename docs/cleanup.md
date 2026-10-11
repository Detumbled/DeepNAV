# Repository cleanup

Removed prior generated results from Output EKF, Output clocks, Output optical,
Output_montecarlo, tests, montecarlo/data, build benchmark result folders, and
the GMAT result reports. Removed the Cassini one-way
C++ demo, its plotting script, study runner, documentation and CMake targets.
Reusable one-way radio models and tests remain; the batch helper now lives in
include/filters/OneWayBatchBaseline.hpp. New trajectory results go under results/.

A recoverable copy of removed files is temporarily stored outside the repository at:

/var/folders/ny/9gsmq_q12v31spm3nfzgjzz80000gn/T/deepnav-previous-results-oeec0joo

Retained only the nominal SAT427/PDS-3220-kg Cassini trajectory. Removed the
alternate-mass/recoil outputs and study driver from the repository; numerical
validation of the retained physical model remains. Recoverable archive:

/tmp/deepnav-sensitivity-archive-8wnawsoh
