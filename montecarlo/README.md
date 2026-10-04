# Monte Carlo studies

This directory is reserved for ensemble runners, study configurations, results,
and plotting scripts. No Monte Carlo runner is implemented yet.

Reuse the `deepnav_ekf` library and observation/dynamics models rather than
duplicating the estimator. Draw initial errors from the declared prior and vary
measurement/process noise seeds. Future validation should compare ensemble state
errors with reported covariance, evaluate NEES/NIS confidence intervals, and
measure the fraction of runs remaining below the navigation requirement over the
entire evaluation interval. Keep results in an ignored `results/` directory.

The single-run EKF demo and plots are separate: `ekf_demo.cpp` and
`tests/plot_ekf.py`.
