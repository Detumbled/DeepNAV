#include "opnav/image/EllipticalGaussian.hpp"
#include <Eigen/Cholesky>
#include <cmath>
#include <iostream>
#include <limits>
#include <numbers>
#include <random>

namespace img=fd::opnav::image;
using fd::opnav::PixelCoordinates;
namespace {
void require(bool condition,const char* message) { if (!condition) throw std::runtime_error(message); }
void near(double actual,double expected,double tolerance=1e-6) {
    require(std::isfinite(actual) && std::abs(actual-expected)<=tolerance,"Numerical comparison failed.");
}
void success(const img::EllipticalGaussianFitResult& fit) {
    require(fit.status==img::FitStatus::Converged && fit.measurement.has_value(),"Ellipse fit failed.");
    require(fit.measurement->covariance.determinant()>0,"Invalid centroid covariance.");
}
template<class F> void rejects(F operation) {
    try { operation(); } catch (const std::invalid_argument&) { return; }
    throw std::runtime_error("Invalid input was accepted.");
}
void integrationAndPartials() {
    const img::EllipticalGaussian m{{5.27,5.43},120,1.3,0.37,0.61,8};
    const PixelCoordinates pixel{5,5};
    const auto analytic=img::evaluateEllipticalGaussianPixel(m,pixel);
    double sum=0;
    constexpr int n=500;
    // Independent midpoint integration in the rotated principal-axis coordinates.
    for (int y=0;y<n;++y) for (int x=0;x<n;++x) {
        const double dx=4.5+(x+.5)/n-m.center.sample,dy=4.5+(y+.5)/n-m.center.line;
        const double u=std::cos(m.angleRadians)*dx+std::sin(m.angleRadians)*dy;
        const double v=-std::sin(m.angleRadians)*dx+std::cos(m.angleRadians)*dy;
        sum += m.heightDn*std::exp(-.5*(u*u/(m.sigmaMajorPixels*m.sigmaMajorPixels)+v*v/(m.sigmaMinorPixels*m.sigmaMinorPixels)));
    }
    near(analytic.dn,m.backgroundDn+sum/(n*n),1e-4);
    for (int k=0;k<7;++k) {
        auto plus=m,minus=m;
        auto change=[k](img::EllipticalGaussian& model,double d) {
            switch(k) {
                case 0:model.center.sample+=d;break; case 1:model.center.line+=d;break;
                case 2:model.heightDn+=d;break; case 3:model.sigmaMajorPixels+=d;break;
                case 4:model.sigmaMinorPixels+=d;break; case 5:model.angleRadians+=d;break;
                case 6:model.backgroundDn+=d;break;
            }
        };
        change(plus,1e-5);change(minus,-1e-5);
        near(analytic.partials[k],(img::evaluateEllipticalGaussianPixel(plus,pixel).dn-
            img::evaluateEllipticalGaussianPixel(minus,pixel).dn)/2e-5,2e-6);
    }
    for (double sigma:{.2,.45,1.3})
        for (double phase:{0.,.49}) {
            const img::EllipticalGaussian round{{5+phase,5.3},120,sigma,sigma,.9,8};
            const auto elliptic=img::evaluateEllipticalGaussianPixel(round,pixel);
            const auto circular=img::evaluateGaussianPixel({round.center,120,sigma,8},pixel);
            near(elliptic.dn,circular.dn,1e-8);
        }
    auto flux=m;flux.center={13.27,13.43};flux.backgroundDn=0;
    const auto raster=img::renderEllipticalGaussian({27,27},flux);
    near(raster.sum(),2*std::numbers::pi*m.heightDn*m.sigmaMajorPixels*m.sigmaMinorPixels,1e-7);
    std::cout << "PASS: independent integration, seven analytic partials, circular limit, total flux\n";
}
void recoveryAndMasks() {
    const PixelCoordinates origin{100,200};
    for (double angle:{0.,.6,-1.1}) {
        const img::EllipticalGaussian truth{{107.23,207.38},1500,1.4,.45,angle,12};
        auto dn=img::renderEllipticalGaussian({15,15},truth,origin);
        auto variance=img::Image::Constant(15,15,4).eval();
        img::PixelMask mask=img::PixelMask::Zero(15,15);
        dn(0,0)=std::numeric_limits<double>::quiet_NaN();variance(0,0)=0;mask(0,0)=1;
        dn(7,8)=1e8;mask(7,8)=1;
        const auto fit=img::fitEllipticalGaussian(dn,variance,origin,&mask);
        success(fit);
        near(fit.model.center.sample,truth.center.sample);
        near(fit.model.center.line,truth.center.line);
        near(fit.model.sigmaMajorPixels,truth.sigmaMajorPixels);
        near(fit.model.sigmaMinorPixels,truth.sigmaMinorPixels);
        near(std::sin(fit.model.angleRadians-truth.angleRadians),0);
        near(fit.model.heightDn,truth.heightDn,1e-3);near(fit.model.backgroundDn,12);
        require(fit.parameterCovariance.has_value(),"Missing physical covariance for anisotropic model.");
        require(fit.degreesOfFreedom==216,"Incorrect masked degrees of freedom.");
    }
    const img::EllipticalGaussian round{{5.21,5.34},1500,.8,.8,0,12};
    const auto dn=img::renderEllipticalGaussian({11,11},round);
    const auto fit=img::fitEllipticalGaussian(dn,img::Image::Constant(11,11,4));
    success(fit);near(fit.model.center.sample,round.center.sample);
    require(!fit.parameterCovariance,"Exactly circular angle must not have a reported uncertainty.");
    // Search with a shifted image origin and a source close to an edge.
    const img::EllipticalGaussian edge{{102.25,203.3},1500,1.,.6,.4,12};
    const auto frame=img::renderEllipticalGaussian({31,25},edge,origin);
    const auto found=img::fitBrightestEllipticalGaussian(frame,img::Image::Constant(25,31,4),5,origin);
    success(found);near(found.model.center.sample,edge.center.sample);near(found.model.center.line,edge.center.line);
    std::cout << "PASS: seven-parameter recovery, rotated/narrow ellipses, masks, origins, edge search, circular degeneracy\n";
}
void noisyCovariance() {
    const img::EllipticalGaussian truth{{5.23,5.38},1500,1.1,.6,.55,12};
    const auto expected=img::renderEllipticalGaussian({11,11},truth);
    const img::Image variance=img::Image::Constant(11,11,4);
    std::mt19937_64 generator(7191);
    std::normal_distribution<double> noise(0,2);
    double normalized=0,fullNormalized=0,reduced=0;
    constexpr int trials=60;
    for (int i=0;i<trials;++i) {
        auto dn=expected;
        for (Eigen::Index j=0;j<dn.size();++j) dn.data()[j]+=noise(generator);
        const auto fit=img::fitEllipticalGaussian(dn,variance);
        success(fit);
        const Eigen::Vector2d error{fit.model.center.sample-truth.center.sample,fit.model.center.line-truth.center.line};
        normalized+=error.dot(fit.measurement->covariance.ldlt().solve(error));
        require(fit.parameterCovariance.has_value(),"Missing seven-parameter covariance.");
        const img::EllipticalParameters parameterError{error[0],error[1],fit.model.heightDn-truth.heightDn,
            fit.model.sigmaMajorPixels-truth.sigmaMajorPixels,fit.model.sigmaMinorPixels-truth.sigmaMinorPixels,
            std::remainder(fit.model.angleRadians-truth.angleRadians,std::numbers::pi),fit.model.backgroundDn-truth.backgroundDn};
        fullNormalized+=parameterError.dot(fit.parameterCovariance->ldlt().solve(parameterError));
        reduced+=fit.chiSquared/fit.degreesOfFreedom;
    }
    normalized/=trials;fullNormalized/=trials;reduced/=trials;
    require(normalized>.8 && normalized<3.5,"Centroid covariance does not match seeded scatter.");
    require(reduced>.9 && reduced<1.1,"Incorrect noise residual statistics.");
    require(fullNormalized>4 && fullNormalized<10,"Physical seven-parameter covariance does not match scatter.");
    std::cout << "PASS: " << trials << " seeded noisy fits; mean normalized squared centroid error=" << normalized
              << ", seven-parameter normalized error=" << fullNormalized << ", reduced chi^2=" << reduced << '\n';
}
void failures() {
    const img::Image flat=img::Image::Constant(11,11,12), variance=img::Image::Constant(11,11,4);
    require(img::fitEllipticalGaussian(flat,variance).status==img::FitStatus::NoSignal,"Flat image was fitted.");
    auto invalid=variance;invalid(0,0)=0;
    rejects([&]{ (void)img::fitEllipticalGaussian(flat,invalid); });
    rejects([&]{ (void)img::renderEllipticalGaussian({11,11},{{5,5},1,1,.1,0,0}); });
    rejects([&]{ (void)img::fitEllipticalGaussian(flat,variance,{},nullptr,{}, {0,1e-7,5}); });
    img::PixelMask mask=img::PixelMask::Ones(11,11);
    rejects([&]{ (void)img::fitEllipticalGaussian(flat,variance,{},&mask); });
    const auto dn=img::renderEllipticalGaussian({11,11},{{5.2,5.3},1500,1.2,.6,.4,12});
    const auto unfinished=img::fitEllipticalGaussian(dn,variance,{},nullptr,{}, {1,1e-12,5});
    require(!unfinished.measurement,"Nonconverged fit produced a measurement.");
    std::cout << "PASS: no signal, invalid data/options/widths, unusable masks, nonconvergence\n";
}
}
int main() {
    try { integrationAndPartials();recoveryAndMasks();noisyCovariance();failures(); }
    catch (const std::exception& error) { std::cerr << "FAIL: " << error.what() << '\n';return 1; }
}
