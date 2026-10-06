# Radiometric modularization performance check

Measured locally on 2026-10-06. The pre-refactor executable was not retained;
the baseline reconstructs its standalone `log1p` Shapiro formula and original
graphics-library linkage. Both demo variants use the same compiled demo object,
explicit Cassini parameters, event algorithms and Debug flags. The comparison
therefore measures the most recent modularization, not Voyager or the initial
implementation before later demo refinements.

The three-count demo ran 20 alternating trials per variant, after warm-up, with
console output redirected. Median complete execution was **117.335 ms before** and
**113.320 ms after** (-3.42%). This includes process startup and local kernel loading;
it is not cold-disk performance.

The kernel microbenchmark compiled both variants with the same AppleClang `-O2`
flags, loaded kernels once, and excluded console output and kernel loading from
timing. Each trial measured nine batches of 10,000 analytic counts or 1,000
SPICE-backed counts per case; four alternating trials per variant were performed.
Each count uses the start/end two-way event chains and 60-second count interval.
The table reports the median of the four trial medians, in microseconds per count.

| Case | Before, us/count | After, us/count | Change |
|---|---:|---:|---:|
| numerical_count | 0.445 | 0.445 | -0.05% |
| spice_CN_count | 129.105 | 132.815 | +2.87% |
| spice_library_count | 137.766 | 140.584 | +2.05% |
| spice_corrected_count | 236.102 | 240.447 | +1.84% |

The analytic numerical model has effectively unchanged timing. SPICE-backed
counts measured approximately 2–3% slower in this run; the CN case also showed
this difference despite unchanged CN calculation code. These small timing
changes should not be attributed solely to formula reuse or modularity without
additional profiling. Complete demo execution was approximately 3% faster.

This is a small local workload check, not a throughput guarantee, a full
orbit-determination benchmark, or a comparison against Voyager's one-way model.
No interpolation cache or parallel SPICE queries were introduced. Real mission
performance depends on observation count, providers, corrections and kernel use.

## State reuse and LSK snapshot optimization

The subsequent source optimization was compared against the retained pre-change
microbenchmark executable (not the earlier reconstructed linkage baseline).
Both benchmark binaries used AppleClang `-O2`, the same kernels, 60-second counts,
and no console output inside timed loops. Two alternating trials per binary each
measured nine batches. Warm kernels were used. These are kernel timings, separate
from the Release build now used for the complete application.

| Case | Before, us/count | Optimized, us/count | Time reduction |
|---|---:|---:|---:|
| CSPICE CN | 137.343 | 123.702 | 9.9% |
| Library with SPICE states | 151.323 | 96.488 | 36.2% |
| Library with corrections | 248.830 | 161.978 | 34.9% |

The receiver is evaluated once per library leg and the converged geometric
emitter is reused. A `SpiceTdbClock` instance reads the three LSK constant groups
once, then evaluates clock offsets without kernel-pool reads. The one-shot
`spiceTdbMinusTt` helper remains available but is not used by the optimized demos.
Clock offsets reused within a count are also evaluated only once.

The very cheap numerical-only case varied substantially: optimized trial
medians ranged from 0.385 to 0.628 us/count, compared with 0.448–0.489 before.
This run does not establish a numerical-only speedup; the robust gains above
apply to the SPICE-backed paths.

`tests/benchmark_two_way.cpp` preserves the benchmark workloads. Build/run with:

```sh
cmake --build build-clang --target benchmark_two_way -j 4
./build-clang/benchmark_two_way
```

The current application uses `CMAKE_BUILD_TYPE=Release` in `build-clang/`; Debug
can be selected by reconfiguring that same directory. No LTO or fast-math flags
were introduced. Receiver-state query counts and clock snapshot/reload semantics
are covered by tests, as well as numerical agreement with CSPICE and analytic
range/Doppler cases.
