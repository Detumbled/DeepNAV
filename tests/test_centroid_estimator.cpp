#include "dynamics/LinearStateProvider.hpp"
#include "opnav/CameraModel.hpp"
#include "opnav/core/ApparentDirection.hpp"
#include "opnav/image/PhotocenterCorrection.hpp"
#include <Eigen/Cholesky>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numbers>

using namespace fd::opnav;
namespace img = fd::opnav::image;

namespace {
void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}
void near(double actual, double expected, double tolerance = 1e-7) {
    require(std::isfinite(actual) && std::abs(actual-expected) <= tolerance, "Numerical comparison failed.");
}
template<class Exception = std::invalid_argument, class F>
void rejects(F f) {
    try { f(); } catch (const Exception&) { return; }
    throw std::runtime_error("Expected exception was not thrown.");
}
void success(const img::GaussianFitResult& fit) {
    require(fit.status == img::FitStatus::Converged && fit.measurement.has_value()
        && fit.parameterCovariance.has_value(), "Gaussian fit did not produce a measurement.");
    require(fit.measurement->covariance.determinant() > 0, "Centroid covariance is not positive definite.");
}

void pixelIntegrationAndPartials() {
    const img::CircularGaussian model{{4.23, 5.41}, 120, 0.65, 8};
    const PixelCoordinates pixel{4, 5};
    const auto analytic = img::evaluateGaussianPixel(model, pixel);
    // Independent midpoint quadrature of continuous intensity across one pixel.
    constexpr int n = 400;
    double integral = 0;
    for (int y = 0; y < n; ++y)
        for (int x = 0; x < n; ++x) {
            const double dx = pixel.sample-0.5+(x+0.5)/n-model.center.sample;
            const double dy = pixel.line-0.5+(y+0.5)/n-model.center.line;
            integral += model.heightDn*std::exp(-(dx*dx+dy*dy)/(2*model.sigmaPixels*model.sigmaPixels));
        }
    near(analytic.dn, model.backgroundDn+integral/(n*n), 1e-4);
    for (int i = 0; i < 5; ++i) {
        auto plus = model, minus = model;
        auto change = [i](img::CircularGaussian& m, double d) {
            switch (i) {
                case 0: m.center.sample += d; break;
                case 1: m.center.line += d; break;
                case 2: m.heightDn += d; break;
                case 3: m.sigmaPixels += d; break;
                case 4: m.backgroundDn += d; break;
            }
        };
        const double step = 1e-5;
        change(plus, step); change(minus, -step);
        const double numeric = (img::evaluateGaussianPixel(plus, pixel).dn
            - img::evaluateGaussianPixel(minus, pixel).dn)/(2*step);
        near(analytic.partials[i], numeric, 1e-6);
    }
    const auto image = img::renderCircularGaussian({31, 31}, {{15.2, 14.6}, 120, 0.65, 0});
    near(image.sum(), 2*std::numbers::pi*120*0.65*0.65, 1e-9);
    std::cout << "PASS: pixel integration, analytic partials, and flux conservation\n";
}

void noiselessRecoveryAndMasks() {
    const PixelCoordinates origin{100, 200};
    for (double sigma : {0.4, 0.7, 1.3}) {
        for (double phase : {0.0, 0.23, 0.49}) {
            const img::CircularGaussian truth{{105+phase, 205.37}, 1800, sigma, 12};
            const auto dn = img::renderCircularGaussian({11, 11}, truth, origin);
            const img::Image variance = img::Image::Constant(11, 11, 9);
            const auto fit = img::fitCircularGaussian(dn, variance, origin);
            success(fit);
            near(fit.model.center.sample, truth.center.sample, 1e-6);
            near(fit.model.center.line, truth.center.line, 1e-6);
            near(fit.model.heightDn, truth.heightDn, 1e-3);
            near(fit.model.sigmaPixels, sigma, 1e-6);
            near(fit.model.backgroundDn, 12, 1e-5);
        }
    }
    const img::CircularGaussian truth{{5.24, 5.36}, 1800, 0.8, 12};
    auto dn = img::renderCircularGaussian({11, 11}, truth);
    auto variance = img::Image::Constant(11, 11, 9).eval();
    img::PixelMask mask = img::PixelMask::Zero(11, 11);
    dn(5, 5) = 1e8; mask(5, 5) = 1; // Saturated/bad central pixel is excluded.
    dn(0, 0) = std::numeric_limits<double>::quiet_NaN();
    variance(0, 0) = 0; mask(0, 0) = 1;
    const auto fit = img::fitCircularGaussian(dn, variance, {}, &mask);
    success(fit);
    near(fit.model.center.sample, truth.center.sample, 1e-6);
    near(fit.model.center.line, truth.center.line, 1e-6);
    require(fit.degreesOfFreedom == 114, "Masked degrees of freedom are wrong.");
    std::cout << "PASS: five-parameter recovery, subpixel phases, narrow PSF, and masks\n";
}

void covarianceMonteCarlo() {
    constexpr int trials = 200;
    const img::CircularGaussian truth{{5.21, 5.38}, 2000, 0.75, 20};
    Eigen::Vector2d sum = Eigen::Vector2d::Zero();
    Eigen::Matrix2d errorSquares = Eigen::Matrix2d::Zero();
    Eigen::Matrix2d predicted = Eigen::Matrix2d::Zero();
    double normalizedError = 0, reducedChiSquared = 0;
    for (int i = 0; i < trials; ++i) {
        const auto image = img::simulateCircularGaussian({11, 11}, truth, {2, 3}, 4000+i);
        const auto fit = img::fitCircularGaussian(image.dn, image.varianceDn2);
        success(fit);
        const Eigen::Vector2d error{fit.model.center.sample-truth.center.sample,
                                   fit.model.center.line-truth.center.line};
        sum += error;
        errorSquares.noalias() += error*error.transpose();
        predicted += fit.measurement->covariance;
        normalizedError += error.dot(fit.measurement->covariance.ldlt().solve(error));
        reducedChiSquared += fit.chiSquared/fit.degreesOfFreedom;
    }
    const Eigen::Vector2d bias = sum/trials;
    const Eigen::Matrix2d scatter = (errorSquares-trials*bias*bias.transpose())/(trials-1);
    predicted /= trials;
    for (int i = 0; i < 2; ++i) {
        require(std::abs(bias[i]) < 4*std::sqrt(predicted(i, i)/trials), "Significant centroid bias.");
        const double ratio = scatter(i, i)/predicted(i, i);
        require(ratio > 0.65 && ratio < 1.45, "Empirical scatter disagrees with covariance.");
    }
    require(normalizedError/trials > 1.4 && normalizedError/trials < 2.6,
            "Two-dimensional normalized errors disagree with covariance.");
    require(reducedChiSquared/trials > 0.85 && reducedChiSquared/trials < 1.15,
            "Residuals disagree with supplied detector noise.");
    std::cout << "PASS: 200 seeded noisy images; mean normalized squared error="
              << normalizedError/trials << ", reduced chi^2=" << reducedChiSquared/trials << '\n';
}

void photocenterCorrection() {
    img::Image brightness = img::Image::Zero(3, 3);
    brightness(1, 1) = 3; brightness(1, 2) = 1;
    auto offset = img::brightnessPhotocenterOffset(brightness, {1, 1});
    near(offset.pixels.x(), 0.25); near(offset.pixels.y(), 0);
    const auto scaled = img::brightnessPhotocenterOffset((brightness*100).eval(), {1, 1});
    near((scaled.pixels-offset.pixels).norm(), 0);
    const auto translated = img::brightnessPhotocenterOffset(brightness, {101, 201}, {100, 200});
    near((translated.pixels-offset.pixels).norm(), 0);
    const img::CircularGaussian source{{5.25, 5.4}, 2000, 0.8, 12};
    const auto fit = img::fitCircularGaussian(img::renderCircularGaussian({11, 11}, source),
                                             img::Image::Constant(11, 11, 4));
    success(fit);
    offset.covariance << 0.0004, 0.0001, 0.0001, 0.0009;
    const auto center = img::correctPhotocenter(*fit.measurement, offset);
    near(center.center.sample, 5, 1e-6); near(center.center.line, 5.4, 1e-6);
    near((center.covariance-fit.measurement->covariance-offset.covariance).norm(), 0);
    near(fit.measurement->center.sample, 5.25, 1e-6); // Raw light center is preserved.
    std::cout << "PASS: brightness first moments, correction sign, and offset uncertainty\n";
}

void failureCases() {
    const img::Image variance = img::Image::Constant(9, 9, 25);
    require(!img::fitCircularGaussian(img::Image::Constant(9, 9, 10), variance).measurement,
            "Blank image produced a measurement.");
    const auto faint = img::fitCircularGaussian(
        img::renderCircularGaussian({9, 9}, {{4.2, 4.3}, 0.1, 0.8, 10}), variance);
    require(faint.status == img::FitStatus::NoSignal && !faint.measurement,
            "Low-SNR fit produced a measurement.");
    const img::CircularGaussian truth{{4.2, 4.3}, 1000, 0.8, 10};
    const auto dn = img::renderCircularGaussian({9, 9}, truth);
    const auto incomplete = img::fitCircularGaussian(dn, variance, {}, nullptr,
        img::CircularGaussian{{3, 3}, 800, 1.2, 5}, {1, 1e-12, 5});
    require(incomplete.status == img::FitStatus::NonConverged && !incomplete.measurement,
            "Nonconverged fit produced a measurement.");
    img::PixelMask mask = img::PixelMask::Ones(9, 9);
    mask.row(4).setZero();
    const auto singular = img::fitCircularGaussian(dn, variance, {}, &mask, truth);
    require(singular.status == img::FitStatus::Singular && !singular.measurement,
            "Unobservable five-parameter fit produced a measurement.");
    rejects([&] { (void)img::fitCircularGaussian(dn, img::Image::Zero(9, 9)); });
    rejects([&] { (void)img::fitCircularGaussian(dn, img::Image::Ones(8, 9)); });
    rejects([&] { (void)img::fitCircularGaussian(dn, variance, {}, nullptr, std::nullopt, {0}); });
    auto bad = dn; bad(0, 0) = std::numeric_limits<double>::quiet_NaN();
    rejects([&] { (void)img::fitCircularGaussian(bad, variance); });
    rejects([&] { (void)img::renderCircularGaussian({9, 9}, {{4, 4}, 1, 0, 0}); });
    rejects([&] { (void)img::simulateCircularGaussian({9, 9}, truth, {0, 1}, 1); });
    rejects([&] { (void)img::brightnessPhotocenterOffset(img::Image::Zero(3, 3), {1, 1}); });
    auto invalidOffset = img::PhotocenterOffset{};
    invalidOffset.covariance << 1, 2, 2, 1;
    rejects([&] { (void)img::correctPhotocenter({{1, 1}, Eigen::Matrix2d::Identity()}, invalidOffset); });
    std::cout << "PASS: no signal, singular/nonconverged fits, and invalid inputs\n";
}

void stateToImageToObservation(const char* demoPath) {
    const TdbEpoch epoch(0);
    const fd::dynamics::CartesianState observer{{0, 0, 0}, {3, 10, 0}};
    const fd::dynamics::LinearStateProvider moon(epoch, {{0, 0, 299792.458}, {0, 0, 0}}, "J2000", "SSB");
    const CameraModel camera(CameraIntrinsics{50, {5, 5}, 20*Eigen::Matrix2d::Identity()});
    const auto apparent = core::computeApparentDirection(epoch, observer, moon);
    const auto pixel = camera.project(apparent.vectorKm, CameraAttitude(Eigen::Matrix3d::Identity()));
    // Unresolved symmetric moon: photocenter offset is zero in this first model.
    const img::CircularGaussian truth{pixel, 5000, 0.7, 30};
    const auto image = img::simulateCircularGaussian({11, 11}, truth, {2, 3}, 42);
    const auto fit = img::fitCircularGaussian(image.dn, image.varianceDn2);
    success(fit);
    const Eigen::Vector2d error{fit.model.center.sample-pixel.sample, fit.model.center.line-pixel.line};
    require(error.dot(fit.measurement->covariance.ldlt().solve(error)) < 25,
            "Synthetic moon centroid disagrees with geometry.");
    std::cout << "PASS: state -> apparent direction -> pixel -> noisy image -> centroid\n"
              << "  true pixel=(" << pixel.sample << ',' << pixel.line << "), fitted=("
              << fit.model.center.sample << ',' << fit.model.center.line << "), error=" << error.norm()
              << " px, sigma=(" << std::sqrt(fit.measurement->covariance(0, 0)) << ','
              << std::sqrt(fit.measurement->covariance(1, 1)) << ") px\n";
    if (demoPath) {
        std::ofstream csv(demoPath);
        csv.exceptions(std::ios::failbit | std::ios::badbit);
        csv << "sample,line,dn,variance_dn2,true_sample,true_line,fit_sample,fit_line,"
               "variance_sample,variance_line,covariance_sample_line\n" << std::setprecision(17);
        for (int l = 0; l < 11; ++l)
            for (int s = 0; s < 11; ++s)
                csv << s << ',' << l << ',' << image.dn(l, s) << ',' << image.varianceDn2(l, s)
                    << ',' << pixel.sample << ',' << pixel.line << ',' << fit.model.center.sample
                    << ',' << fit.model.center.line << ',' << fit.measurement->covariance(0, 0)
                    << ',' << fit.measurement->covariance(1, 1) << ',' << fit.measurement->covariance(0, 1) << '\n';
    }
}
} // namespace

int main(int argc, char** argv) {
    try {
        if (argc > 2) throw std::invalid_argument("Usage: test_centroid_estimator [demo.csv]");
        pixelIntegrationAndPartials();
        noiselessRecoveryAndMasks();
        covarianceMonteCarlo();
        photocenterCorrection();
        failureCases();
        stateToImageToObservation(argc == 2 ? argv[1] : nullptr);
        std::cout << "All image-estimator tests passed.\n";
    } catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n'; return 1;
    }
}
