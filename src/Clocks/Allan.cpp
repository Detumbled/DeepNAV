#include "Clocks/Allan.hpp"
#include "Validation.hpp"

#include <Eigen/QR>

#include <cmath>
#include <limits>
#include <stdexcept>

namespace fd::clocks {

using namespace detail;

ClockParameters fitClockAllanData(
    std::span<const AllanDatum> data, const AllanFitAssumptions& a) {
    positive(a.valid_tau_min_s);
    positive(a.valid_tau_max_s);
    nonnegative(a.max_relative_variance_residual);
    finite(a.frequency_drift_per_s);
    if (!a.drift_removed || a.valid_tau_min_s > a.valid_tau_max_s || data.empty())
        throw std::invalid_argument("Allan fit needs data, a valid tau interval and explicit drift removal.");
    if (a.noise != ClockNoiseAssumption::Both
        && a.noise != ClockNoiseAssumption::WhiteFrequencyOnly
        && a.noise != ClockNoiseAssumption::RandomWalkFrequencyOnly)
        throw std::invalid_argument("Unknown clock noise assumption.");
    const bool both = a.noise == ClockNoiseAssumption::Both;
    if (both && data.size() < 2)
        throw std::invalid_argument("One Allan point cannot identify two noise intensities.");

    // Scale columns and observations before QR: physical q_b and q_y may
    // differ by many orders of magnitude. No optimization dependency needed.
    Eigen::MatrixXd design(data.size(), both ? 2 : 1);
    Eigen::VectorXd variance(data.size());
    for (std::size_t i = 0; i < data.size(); ++i) {
        positive(data[i].tau_s);
        nonnegative(data[i].adev);
        if (data[i].tau_s < a.valid_tau_min_s || data[i].tau_s > a.valid_tau_max_s)
            throw std::invalid_argument("Allan datum is outside the declared fit interval.");
        variance[i] = checked(data[i].adev * data[i].adev);
        design(i, 0) = checked(a.noise == ClockNoiseAssumption::RandomWalkFrequencyOnly
            ? data[i].tau_s / 3.0 : 1.0 / data[i].tau_s);
        if (both) design(i, 1) = data[i].tau_s / 3.0;
    }
    Eigen::VectorXd scales(design.cols());
    for (Eigen::Index j = 0; j < design.cols(); ++j) {
        scales[j] = design.col(j).maxCoeff();
        design.col(j) /= scales[j];
    }
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(design);
    qr.setThreshold(1.0e-8);
    if (qr.rank() != design.cols())
        throw std::invalid_argument("Allan fit is ill-conditioned; use distinct, well-separated tau values.");
    const double variance_scale = variance.maxCoeff();
    Eigen::VectorXd coefficients = Eigen::VectorXd::Zero(design.cols());
    if (variance_scale > 0.0) coefficients = qr.solve(variance / variance_scale);
    if (!coefficients.allFinite() || (coefficients.array() < -1.0e-12).any())
        throw std::invalid_argument("Allan data require negative noise intensities; incompatible clock model.");
    coefficients = coefficients.cwiseMax(0.0); // Only QR roundoff at the zero boundary.
    const Eigen::VectorXd predicted = design * coefficients * variance_scale;
    for (Eigen::Index i = 0; i < variance.size(); ++i) {
        const double allowance = a.max_relative_variance_residual * variance[i]
            + 64.0 * std::numeric_limits<double>::epsilon() * variance_scale;
        if (std::abs(predicted[i] - variance[i]) > allowance)
            throw std::invalid_argument("Allan variance residual exceeds tolerance; incompatible clock model.");
    }
    ClockParameters p{a.frequency_drift_per_s, 0.0, 0.0};
    if (a.noise == ClockNoiseAssumption::RandomWalkFrequencyOnly)
        p.q_frequency_per_s = checked(coefficients[0] * variance_scale / scales[0]);
    else {
        p.q_bias_s = checked(coefficients[0] * variance_scale / scales[0]);
        if (both) p.q_frequency_per_s = checked(coefficients[1] * variance_scale / scales[1]);
    }
    validate(p);
    return p;
}

double overlappingAllanDeviation(std::span<const double> bias, double dt, std::size_t m) {
    positive(dt);
    if (m == 0 || bias.empty() || m > (bias.size() - 1) / 2)
        throw std::invalid_argument("Allan deviation requires m >= 1 and N > 2m.");
    for (double b : bias) finite(b);
    const long double tau = static_cast<long double>(m) * dt;
    long double sum = 0.0;
    for (std::size_t i = 0; i < bias.size() - 2 * m; ++i) {
        const long double difference = static_cast<long double>(bias[i + 2 * m])
            - 2.0L * bias[i + m] + bias[i];
        sum += difference * difference;
    }
    return checked(static_cast<double>(std::sqrt(sum / (2.0L * tau * tau * (bias.size() - 2 * m)))));
}

double theoreticalAllanDeviation(double tau, const ClockParameters& p, bool include_drift) {
    positive(tau);
    validate(p);
    const double drift = include_drift ? p.frequency_drift_per_s * tau : 0.0;
    return checked(std::sqrt(p.q_bias_s / tau + p.q_frequency_per_s * tau / 3.0 + 0.5 * drift * drift));
}

} // namespace fd::clocks
