#!/usr/bin/env python3
"""Self-contained clock-model checks. Requires only NumPy; produces no files/plots.

Run: python3 test_clocks.py
All coefficients below are synthetic validation values, not hardware specifications.
"""
import numpy as np

C = 299_792_458.0


def transition(dt, q_b, q_y):
    if not np.isfinite([dt, q_b, q_y]).all() or dt <= 0 or min(q_b, q_y) < 0:
        raise ValueError("Require finite dt > 0 and q_b, q_y >= 0")
    f = np.array([[1.0, dt], [0.0, 1.0]])
    q = np.array([[q_b * dt + q_y * dt**3 / 3, q_y * dt**2 / 2],
                  [q_y * dt**2 / 2, q_y * dt]])
    return f, q


def draw_increments(count, dt, q_b, q_y, rng):
    """Exact jointly distributed bias/frequency increments, including zero noise."""
    transition(dt, q_b, q_y)  # Validate inputs.
    z_b, z_y, z_i = rng.standard_normal((3, count))
    dy = np.sqrt(q_y * dt) * z_y
    db = (np.sqrt(q_b * dt) * z_b
          + np.sqrt(q_y * dt**3) * (z_y / 2 + z_i / np.sqrt(12)))
    return db, dy


def simulate(steps, dt, q_b, q_y, b0=0.0, y0=0.0, drift=0.0, seed=2026):
    """Return ideal elapsed time, bias [s], fractional frequency [dimensionless]."""
    if not isinstance(steps, (int, np.integer)) or steps < 1:
        raise ValueError("steps must be a positive integer")
    if not np.isfinite([b0, y0, drift]).all():
        raise ValueError("Initial state and drift must be finite")
    db, dy = draw_increments(steps, dt, q_b, q_y, np.random.default_rng(seed))
    y = np.r_[y0, y0 + np.cumsum(drift * dt + dy)]
    b = np.r_[b0, b0 + np.cumsum(y[:-1] * dt + 0.5 * drift * dt**2 + db)]
    t = np.arange(steps + 1) * dt
    return t, b, y


def overlapping_adev(bias, dt, m):
    """Overlapping Allan deviation from uniformly sampled time error [s]."""
    bias = np.asarray(bias, dtype=float)
    if (bias.ndim != 1 or not np.isfinite(bias).all() or not np.isfinite(dt)
            or dt <= 0 or not isinstance(m, (int, np.integer))
            or m < 1 or 2 * m >= len(bias)):
        raise ValueError("Need finite 1D phase data, dt > 0, integer m >= 1, N > 2m")
    second_difference = bias[2*m:] - 2 * bias[m:-m] + bias[:-2*m]
    return np.sqrt(np.mean(second_difference**2) / (2 * (m * dt)**2))


def theoretical_adev(tau, q_b, q_y, drift=0.0):
    # Drift contributes only when NOT removed from the phase series.
    return np.sqrt(q_b / tau + q_y * tau / 3 + 0.5 * (drift * tau)**2)


def range_envelope(t, q_b, q_y, b0=0.0, y0=0.0, drift=0.0,
                   sigma_b0=0.0, sigma_y0=0.0, k=3.0):
    """Conservative pointwise bound; independent initial b/y uncertainties.

    Deterministic terms use absolute values to prevent cancellation in the bound.
    This is NOT a probability bound for the entire sample path.
    """
    deterministic = abs(b0) + abs(y0) * t + 0.5 * abs(drift) * t**2
    variance = sigma_b0**2 + (sigma_y0 * t)**2 + q_b * t + q_y * t**3 / 3
    return C * (deterministic + k * np.sqrt(variance))


def calibration_interval(threshold_m, horizon_s, **parameters):
    """First threshold contact for the monotonic conservative envelope.

    None means no crossing within horizon, not unlimited autonomy.
    """
    if (not np.isfinite([threshold_m, horizon_s]).all()
            or min(threshold_m, horizon_s) <= 0):
        raise ValueError("threshold and horizon must be finite and positive")
    for name, value in parameters.items():
        if not np.isfinite(value):
            raise ValueError(f"{name} must be finite")
        if name in {"q_b", "q_y", "sigma_b0", "sigma_y0", "k"} and value < 0:
            raise ValueError(f"{name} must be nonnegative")
    if range_envelope(0.0, **parameters) >= threshold_m:
        return 0.0
    if range_envelope(horizon_s, **parameters) < threshold_m:
        return None
    lo, hi = 0.0, horizon_s
    for _ in range(80):
        mid = (lo + hi) / 2
        if range_envelope(mid, **parameters) < threshold_m:
            lo = mid
        else:
            hi = mid
    return hi


def check(condition, message):
    if not condition:
        raise AssertionError(message)


def main():
    # 1. Exact deterministic evolution; no noise and no matrix factorization issues.
    t, b, y = simulate(1000, 2.0, 0.0, 0.0, b0=2e-9, y0=1e-11, drift=2e-15)
    np.testing.assert_allclose(b, 2e-9 + 1e-11*t + 0.5*2e-15*t**2,
                               rtol=1e-12, atol=1e-20)
    np.testing.assert_allclose(y, 1e-11 + 2e-15*t, rtol=1e-12, atol=1e-23)
    print("PASS: deterministic bias, frequency and drift")

    # 2. Independent empirical verification of the exact process-noise covariance.
    dt, qb, qy = 10.0, 4e-24, 3e-26
    _, expected_q = transition(dt, qb, qy)
    db, dy = draw_increments(200_000, dt, qb, qy, np.random.default_rng(42))
    measured_q = np.cov(np.vstack((db, dy)), ddof=1)
    np.testing.assert_allclose(measured_q, expected_q, rtol=0.03, atol=0.0)
    print("PASS: empirical Q, including bias/frequency cross-covariance")

    # 3. Finite-record Allan checks at tau = 1..256 s. Each case uses its own series.
    # 20% is a finite-sample tolerance, not a specification of clock accuracy.
    cases = [("white frequency", 1e-24, 0.0),
             ("random walk frequency", 0.0, 3e-26),
             ("combined", 1e-24, 3e-26)]
    for index, (name, qb, qy) in enumerate(cases):
        t, b, _ = simulate(2**18, 1.0, qb, qy, seed=100 + index)
        ratios = []
        for m in (1, 4, 16, 64, 256):
            measured = overlapping_adev(b, 1.0, m)
            expected = theoretical_adev(float(m), qb, qy)
            ratio = measured / expected
            check(abs(ratio - 1) < 0.20, f"{name}, tau={m}: ratio={ratio}")
            ratios.append(f"{m}s:{ratio:.3f}")
        print(f"PASS: {name}; simulated/theoretical ADEV: " + ", ".join(ratios))

    # 4. Allan invariance to offset + constant frequency, and effect of drift.
    linear = 1e-8 + 2e-11*t
    check(overlapping_adev(linear, 1.0, 64) < 1e-22,
          "Constant bias and frequency should not contribute to ADEV")
    drift = 2e-15
    quadratic = 0.5 * drift * t**2
    np.testing.assert_allclose(overlapping_adev(quadratic, 1.0, 64),
                               abs(drift)*64/np.sqrt(2), rtol=1e-7)
    print("PASS: Allan offset invariance and deterministic drift contribution")

    # 5. Covariance propagation against the continuous-time analytic result.
    qb, qy, dt = 1e-24, 3e-26, 7.0
    f, q = transition(dt, qb, qy)
    p = np.zeros((2, 2))
    for _ in range(100):
        p = f @ p @ f.T + q
    duration = 100*dt
    np.testing.assert_allclose(p[0, 0], qb*duration + qy*duration**3/3, rtol=1e-12)
    print("PASS: propagated timing variance")

    # 6. Calibrating the calibration-interval calculation itself.
    interval = calibration_interval(1.0, 30*86400.0, q_b=0.0, q_y=0.0, y0=3e-15)
    np.testing.assert_allclose(interval, 1/(C*3e-15), rtol=1e-12)
    check(calibration_interval(1.0, 86400.0, q_b=0.0, q_y=0.0) is None,
          "Ideal clock must not cross the threshold")
    check(calibration_interval(1.0, 86400.0, q_b=0.0, q_y=0.0, b0=1e-8) == 0,
          "Initial bias already exceeds budget")
    print("PASS: threshold calculation, zero-noise and exceeded-initial-budget cases")
    print(f"Illustration ONLY: constant residual y=3e-15 -> 1 m in {interval/86400:.2f} days")
    print("All checks passed. No plots or output files generated.")


if __name__ == "__main__":
    main()
