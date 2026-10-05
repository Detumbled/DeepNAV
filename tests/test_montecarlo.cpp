#include "Statistics.hpp"
#include <Eigen/Core>
#include <cmath>
#include <iostream>

namespace {
void require(bool condition, const char* message) {
    if (!condition)
        throw std::runtime_error(message);
}
} // namespace

int main() {
    try {
        Eigen::Matrix2d covariance;
        covariance << 4, 1, 1, 9;
        const Eigen::Vector2d error(1, 2);
        require(std::abs(fd::montecarlo::nees(error, covariance) - .6) < 1e-12,
                "NEES lost the covariance cross terms.");
        const Eigen::Vector2d units(1000, 1e-10);
        require(
            std::abs(fd::montecarlo::nees(units.asDiagonal() * error,
                                          units.asDiagonal() * covariance * units.asDiagonal()) -
                     .6) < 1e-12,
            "NEES depends on physical units.");
        std::mt19937_64 generator(1234), repeated(1234);
        require((fd::montecarlo::drawError(covariance, generator) -
                 fd::montecarlo::drawError(covariance, repeated))
                        .norm() == 0,
                "Prior sampling is not reproducible.");
        Eigen::Vector2d sum = Eigen::Vector2d::Zero();
        Eigen::Matrix2d products = Eigen::Matrix2d::Zero();
        double meanNees = 0;
        constexpr int count = 20000;
        for (int j = 0; j < count; ++j) {
            const Eigen::Vector2d sample = fd::montecarlo::drawError(covariance, generator);
            sum += sample;
            products += sample * sample.transpose();
            meanNees += fd::montecarlo::nees(sample, covariance);
        }
        const Eigen::Vector2d mean = sum / count;
        const Eigen::Matrix2d empirical = products / count - mean * mean.transpose();
        require(mean.norm() < .05 && (empirical - covariance).norm() / covariance.norm() < .04 &&
                    std::abs(meanNees / count - 2) < .05,
                "Correlated prior sampling does not reproduce covariance/expected NEES.");
        bool rejected = false;
        try {
            (void)fd::montecarlo::nees(error, Eigen::Matrix2d::Zero());
        } catch (const std::invalid_argument&) {
            rejected = true;
        }
        require(rejected, "Singular covariance must not silently produce a NEES.");
        std::cout
            << "Correlated sampling, reproducibility, mixed-unit NEES and validation passed.\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
