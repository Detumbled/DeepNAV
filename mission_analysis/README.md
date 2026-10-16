# Mission analysis

Reserved for trajectory assessment, maneuver planning and mission-level trade
studies. No mission-analysis algorithms or build targets are implemented yet.

- `include/mission_analysis/`: future public C++ interfaces (`fd::mission_analysis`).
- `src/`: implementations.
- `tests/`: module tests.

Reuse DeepNAV's dynamics, propagation and observation libraries. Keep
mission-specific inputs and studies in `benchmarks/`, and generated outputs in
`results/`. Add this module to the root CMake project when it has source files.
