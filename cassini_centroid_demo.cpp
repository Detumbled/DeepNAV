#include "opnav/image/CassiniImage.hpp"
#include "opnav/image/CentroidEstimator.hpp"
#include "opnav/image/EllipticalGaussian.hpp"
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

namespace img = fd::opnav::image;
namespace {
double median(std::vector<double> values) {
    if (values.empty()) throw std::runtime_error("No usable background pixels.");
    const auto mid = values.begin()+values.size()/2;
    std::nth_element(values.begin(), mid, values.end());
    return values.size()%2 ? *mid : 0.5*(*mid+*std::max_element(values.begin(), mid));
}
std::ofstream output(const std::filesystem::path& path) {
    std::ofstream stream(path);
    if (!stream) throw std::runtime_error("Cannot write " + path.string());
    stream << std::setprecision(17);
    stream.exceptions(std::ios::badbit | std::ios::failbit);
    return stream;
}
fd::opnav::PixelCoordinates brightnessCentroid(const img::Image& window, const img::PixelMask& mask,
    fd::opnav::PixelCoordinates origin, double background) {
    double sum=0, x=0, y=0;
    for (int l=0;l<window.rows();++l) for (int s=0;s<window.cols();++s) {
        if (mask(l,s)) continue;
        const double signal=std::max(0.0,window(l,s)-background);
        sum+=signal; x+=signal*s; y+=signal*l;
    }
    if (!(sum>0)) throw std::runtime_error("No positive signal for brightness centroid.");
    return {origin.sample+x/sum,origin.line+y/sum};
}
}

int main(int argc, char** argv) {
    try {
        const std::filesystem::path root = DEEPNAV_SOURCE_DIR;
        std::filesystem::path input = root/"Input Optical/N1476124698_2.IMG";
        std::filesystem::path directory = root/"Output optical/cassini";
        int sample = -1, line = -1, radius = 10;
        for (int i = 1; i < argc; ++i) {
            const std::string arg = argv[i];
            if (arg == "--help") {
                std::cout << "cassini_centroid_demo [image.IMG] [--sample X --line Y] [--radius R] [--output DIR]\n"
                          << "Coordinates are zero-based original image pixels; default seed is brightest usable pixel.\n";
                return 0;
            }
            if (arg == "--sample" || arg == "--line" || arg == "--radius" || arg == "--output") {
                if (++i == argc) throw std::invalid_argument("Missing option value.");
                if (arg == "--output") directory = argv[i];
                else {
                    std::size_t used = 0;
                    const int value = std::stoi(argv[i], &used);
                    if (used != std::string(argv[i]).size() || value < 0)
                        throw std::invalid_argument("Expected a nonnegative integer option.");
                    if (arg == "--sample") sample = value;
                    else if (arg == "--line") line = value;
                    else radius = value;
                }
            } else if (arg.starts_with("--")) throw std::invalid_argument("Unknown option: " + arg);
            else input = arg;
        }
        const auto image = img::loadCassiniImage(input);
        const auto& dn = image.dn;
        if ((sample < 0) != (line < 0)) throw std::invalid_argument("Provide both --sample and --line.");
        if (sample < 0) {
            double maximum = -1;
            for (int y = 0; y < dn.rows(); ++y)
                for (int x = 0; x < dn.cols(); ++x)
                    if (dn(y,x) >= 0 && dn(y,x) < 4095 && dn(y,x) > maximum) {
                        maximum = dn(y,x); sample = x; line = y;
                    }
        }
        if (sample < 0 || sample >= dn.cols() || line < 0 || line >= dn.rows() ||
            radius < 2 || radius > 50)
            throw std::invalid_argument("Seed must be inside the image; radius must be 2..50.");
        const int x0 = std::max(0, sample-radius), y0 = std::max(0, line-radius);
        const int width = std::min(int(dn.cols())-1, sample+radius)-x0+1;
        const int height = std::min(int(dn.rows())-1, line+radius)-y0+1;
        const img::Image window = dn.block(y0,x0,height,width);
        // Estimate local baseline variance outside the fit window. It includes
        // background fluctuations; add only source shot noise, avoiding double counting.
        std::vector<double> background;
        const int inner = radius+3, outer = radius+15;
        for (int y = std::max(0,line-outer); y <= std::min(int(dn.rows())-1,line+outer); ++y)
            for (int x = std::max(0,sample-outer); x <= std::min(int(dn.cols())-1,sample+outer); ++x)
                if (std::max(std::abs(x-sample),std::abs(y-line)) >= inner && dn(y,x)>=0 && dn(y,x)<4095)
                    background.push_back(dn(y,x));
        const double baseline = median(background);
        for (auto& value : background) value = std::abs(value-baseline);
        const double noise = 1.4826*median(background);
        const double baselineVariance = std::max(1.0/12.0,noise*noise);
        const img::Image variance = (window.array()-baseline).max(0)/image.electronsPerDn + baselineVariance;
        img::PixelMask mask(height,width);
        for (int y = 0; y < height; ++y)
            for (int x = 0; x < width; ++x) mask(y,x) = window(y,x)<0 || window(y,x)>=4095;
        const fd::opnav::PixelCoordinates origin{double(x0),double(y0)};
        const auto fit = img::fitCircularGaussian(window,variance,origin,&mask);
        if (!fit.measurement) throw std::runtime_error("Gaussian fit failed; status=" + std::to_string(int(fit.status)));
        const auto model = img::renderCircularGaussian({width,height},fit.model,origin);
        const auto ellipse = img::fitEllipticalGaussian(window,variance,origin,&mask);
        if (!ellipse.measurement) throw std::runtime_error("Elliptical fit failed; status=" + std::to_string(int(ellipse.status)));
        const auto ellipticalModel = img::renderEllipticalGaussian({width,height},ellipse.model,origin);
        const auto moment = brightnessCentroid(window,mask,origin,baseline);
        std::filesystem::create_directories(directory);
        auto full = output(directory/"image.csv");
        for (int y = 0; y < dn.rows(); ++y) {
            for (int x = 0; x < dn.cols(); ++x) full << (x ? "," : "") << dn(y,x);
            full << '\n';
        }
        auto roi = output(directory/"fit.csv");
        roi << "sample,line,dn,model_dn,residual_dn,variance_dn2,masked,elliptical_model_dn,elliptical_residual_dn\n";
        for (int y = 0; y < height; ++y)
            for (int x = 0; x < width; ++x)
                roi << x0+x << ',' << y0+y << ',' << window(y,x) << ',' << model(y,x) << ','
                    << window(y,x)-model(y,x) << ',' << variance(y,x) << ',' << int(mask(y,x)) << ','
                    << ellipticalModel(y,x) << ',' << window(y,x)-ellipticalModel(y,x) << '\n';
        const auto& cov = fit.measurement->covariance;
        const double reduced = fit.chiSquared/fit.degreesOfFreedom;
        auto summary = output(directory/"summary.json");
        summary << "{\n\"image\": " << std::quoted(input.filename().string())
                << ",\n\"mid_time_utc\": " << std::quoted(image.imageMidTimeUtc)
                << ",\n\"target_description\": " << std::quoted(image.targetDescription)
                << ",\n\"sample\": " << fit.model.center.sample << ",\n\"line\": " << fit.model.center.line
                << ",\n\"height_dn\": " << fit.model.heightDn << ",\n\"sigma_pixels\": " << fit.model.sigmaPixels
                << ",\n\"background_dn\": " << fit.model.backgroundDn
                << ",\n\"background_noise_dn\": " << std::sqrt(baselineVariance)
                << ",\n\"gain_electrons_per_dn\": " << image.electronsPerDn
                << ",\n\"exposure_ms\": " << image.exposureMilliseconds
                << ",\n\"reduced_chi_squared\": " << reduced
                << ",\n\"iterations\": " << fit.iterations
                << ",\n\"covariance\": [[" << cov(0,0) << ',' << cov(0,1) << "],[" << cov(1,0) << ',' << cov(1,1)
                << "]],\n\"brightness_centroid\": [" << moment.sample << ',' << moment.line << ']'
                << ",\n\"elliptical\": {\"sample\": " << ellipse.model.center.sample
                << ", \"line\": " << ellipse.model.center.line << ", \"height_dn\": " << ellipse.model.heightDn
                << ", \"sigma_major_pixels\": " << ellipse.model.sigmaMajorPixels
                << ", \"sigma_minor_pixels\": " << ellipse.model.sigmaMinorPixels
                << ", \"angle_radians\": " << ellipse.model.angleRadians
                << ", \"background_dn\": " << ellipse.model.backgroundDn
                << ", \"reduced_chi_squared\": " << ellipse.chiSquared/ellipse.degreesOfFreedom
                << ", \"iterations\": " << ellipse.iterations << ", \"covariance\": [["
                << ellipse.measurement->covariance(0,0) << ',' << ellipse.measurement->covariance(0,1) << "],["
                << ellipse.measurement->covariance(1,0) << ',' << ellipse.measurement->covariance(1,1) << "]]},"
                << "\n\"noise_model\": \"local background MAD plus source DN / label gain; uncalibrated diagnostic\"\n}\n";
        auto stability=output(directory/"window_comparison.csv");
        stability << "radius,circular_sample,circular_line,elliptical_sample,elliptical_line,moment_sample,moment_line,circular_reduced_chi_squared,elliptical_reduced_chi_squared,circular_status,elliptical_status\n";
        std::vector<int> radii{std::max(2,radius-2),radius,std::min(50,radius+2)};
        radii.erase(std::unique(radii.begin(),radii.end()),radii.end());
        for (const int r:radii) {
            const int sx=std::max(0,sample-r),sy=std::max(0,line-r);
            const int w=std::min(int(dn.cols())-1,sample+r)-sx+1;
            const int h=std::min(int(dn.rows())-1,line+r)-sy+1;
            const img::Image patch=dn.block(sy,sx,h,w);
            const img::Image weights=(patch.array()-baseline).max(0)/image.electronsPerDn+baselineVariance;
            img::PixelMask excluded(h,w);
            for (int y=0;y<h;++y) for (int x=0;x<w;++x) excluded(y,x)=patch(y,x)<0 || patch(y,x)>=4095;
            const fd::opnav::PixelCoordinates offset{double(sx),double(sy)};
            const auto circular = r==radius ? fit : img::fitCircularGaussian(patch,weights,offset,&excluded);
            const auto elliptical = r==radius ? ellipse : img::fitEllipticalGaussian(patch,weights,offset,&excluded);
            const auto centroid=brightnessCentroid(patch,excluded,offset,baseline);
            const double missing=std::numeric_limits<double>::quiet_NaN();
            stability << r << ',' << (circular.measurement ? circular.model.center.sample : missing) << ','
                << (circular.measurement ? circular.model.center.line : missing) << ','
                << (elliptical.measurement ? elliptical.model.center.sample : missing) << ','
                << (elliptical.measurement ? elliptical.model.center.line : missing) << ','
                << centroid.sample << ',' << centroid.line << ','
                << circular.chiSquared/circular.degreesOfFreedom << ',' << elliptical.chiSquared/elliptical.degreesOfFreedom << ','
                << int(circular.status) << ',' << int(elliptical.status) << '\n';
        }
        std::cout << std::setprecision(10) << "Loaded " << dn.cols() << " x " << dn.rows() << " original DN pixels\n"
                  << "Fitted centroid (sample,line): (" << fit.model.center.sample << ',' << fit.model.center.line << ")\n"
                  << "Sigma: " << fit.model.sigmaPixels << " px; reduced chi^2: " << reduced << '\n'
                  << "Elliptical centroid: (" << ellipse.model.center.sample << ',' << ellipse.model.center.line << ")\n"
                  << "Elliptical widths: (" << ellipse.model.sigmaMajorPixels << ',' << ellipse.model.sigmaMinorPixels
                  << ") px; angle: " << ellipse.model.angleRadians << " rad; reduced chi^2: "
                  << ellipse.chiSquared/ellipse.degreesOfFreedom << '\n'
                  << "Background-subtracted brightness centroid: (" << moment.sample << ',' << moment.line << ")\n"
                  << "Saved diagnostic files in " << directory << '\n'
                  << "Noise and covariance are estimates from uncalibrated EDR data.\n";
        if (reduced > 3)
            std::cout << "Large residuals: review source shape and noise model before using this centroid for navigation.\n";
    } catch (const std::exception& error) {
        std::cerr << "ERROR: " << error.what() << '\n';
        return 1;
    }
}
