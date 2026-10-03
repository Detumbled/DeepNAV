#include "EllipticalGaussianInternal.hpp"
#include <Eigen/Cholesky>
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <limits>

namespace fd::opnav::image {
namespace {
using Vector = EllipticalParameters;
using Matrix = EllipticalParameterCovariance;

Eigen::Matrix2d precision(const Vector& p) {
    const double a = std::exp(p[3]), b = p[4], c = std::exp(p[5]);
    Eigen::Matrix2d q;
    q << a*a,a*b,a*b,b*b+c*c;
    return q;
}
Vector pack(const EllipticalGaussian& m) {
    const Eigen::Matrix2d lower = detail::precision(m).llt().matrixL();
    return Vector{m.center.sample,m.center.line,std::log(m.heightDn),
        std::log(lower(0,0)),lower(1,0),std::log(lower(1,1)),m.backgroundDn};
}
EllipticalGaussian unpack(const Vector& p) {
    const Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> eigen(precision(p));
    if (eigen.info()!=Eigen::Success || eigen.eigenvalues()[0]<=0)
        throw std::overflow_error("Invalid elliptical precision matrix.");
    const auto axis = eigen.eigenvectors().col(0);
    double angle = std::atan2(axis[1],axis[0]);
    if (angle >= std::numbers::pi/2) angle -= std::numbers::pi;
    if (angle < -std::numbers::pi/2) angle += std::numbers::pi;
    return {{p[0],p[1]},std::exp(p[2]),1/std::sqrt(eigen.eigenvalues()[0]),
        1/std::sqrt(eigen.eigenvalues()[1]),angle,p[6]};
}
bool inWindow(const EllipticalGaussian& m, const Image& dn, PixelCoordinates origin) {
    return m.center.sample >= origin.sample-0.5 && m.center.sample <= origin.sample+dn.cols()-0.5
        && m.center.line >= origin.line-0.5 && m.center.line <= origin.line+dn.rows()-0.5;
}
std::optional<EllipticalGaussian> initialize(const Image& dn, PixelCoordinates origin, const PixelMask* mask) {
    std::vector<double> border, all;
    for (Eigen::Index y = 0; y < dn.rows(); ++y)
        for (Eigen::Index x = 0; x < dn.cols(); ++x) {
            if (mask && (*mask)(y,x)) continue;
            all.push_back(dn(y,x));
            if (y==0 || x==0 || y==dn.rows()-1 || x==dn.cols()-1) border.push_back(dn(y,x));
        }
    auto& values = border.empty() ? all : border;
    auto middle = values.begin()+values.size()/2;
    std::nth_element(values.begin(),middle,values.end());
    const double background = *middle;
    double sum = 0;
    Eigen::Vector2d center = Eigen::Vector2d::Zero();
    for (Eigen::Index y = 0; y < dn.rows(); ++y)
        for (Eigen::Index x = 0; x < dn.cols(); ++x) {
            if (mask && (*mask)(y,x)) continue;
            const double weight = std::max(0.0,dn(y,x)-background);
            sum += weight; center += weight*Eigen::Vector2d{double(x),double(y)};
        }
    if (!(sum>0) || !std::isfinite(sum)) return std::nullopt;
    center /= sum;
    Eigen::Matrix2d moment = Eigen::Matrix2d::Zero();
    for (Eigen::Index y = 0; y < dn.rows(); ++y)
        for (Eigen::Index x = 0; x < dn.cols(); ++x) {
            if (mask && (*mask)(y,x)) continue;
            const double weight = std::max(0.0,dn(y,x)-background);
            const Eigen::Vector2d offset = Eigen::Vector2d{double(x),double(y)}-center;
            moment.noalias() += weight*offset*offset.transpose();
        }
    moment /= sum;
    moment.diagonal().array() -= 1.0/12; // Remove approximate square-pixel broadening.
    const Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> eigen(moment);
    if (eigen.info()!=Eigen::Success) throw std::overflow_error("Elliptical initializer overflowed.");
    const double maximum = 0.5*std::max(dn.rows(),dn.cols());
    const double major = std::clamp(std::sqrt(std::max(0.16,eigen.eigenvalues()[1])),0.4,maximum);
    const double minor = std::clamp(std::sqrt(std::max(0.16,eigen.eigenvalues()[0])),0.4,maximum);
    const auto axis = eigen.eigenvectors().col(1);
    return EllipticalGaussian{{origin.sample+center[0],origin.line+center[1]},
        sum/(2*std::numbers::pi*major*minor),major,minor,std::atan2(axis[1],axis[0]),background};
}
struct Normal {
    Matrix information{Matrix::Zero()};
    Vector gradient{Vector::Zero()};
    double cost{0};
};
Normal evaluate(const Image& dn, const Image& variance, const PixelMask* mask, PixelCoordinates origin,
                const Vector& p) {
    const auto q = precision(p);
    const double height = std::exp(p[2]), a = std::exp(p[3]), b = p[4], c = std::exp(p[5]);
    Normal result;
    for (Eigen::Index y = 0; y < dn.rows(); ++y)
        for (Eigen::Index x = 0; x < dn.cols(); ++x) {
            if (mask && (*mask)(y,x)) continue;
            auto pixel = detail::integrateEllipse({p[0],p[1]},height,p[6],q,{origin.sample+x,origin.line+y});
            const auto shape = pixel.partials.segment<3>(3).eval();
            pixel.partials[2] *= height;
            pixel.partials[3] = shape[0]*2*a*a+shape[1]*a*b;
            pixel.partials[4] = shape[1]*a+shape[2]*2*b;
            pixel.partials[5] = shape[2]*2*c*c;
            const double residual = dn(y,x)-pixel.dn, weight = 1/variance(y,x);
            result.cost += weight*residual*residual;
            result.gradient.noalias() += weight*residual*pixel.partials;
            result.information.noalias() += weight*pixel.partials*pixel.partials.transpose();
        }
    if (!std::isfinite(result.cost) || !result.gradient.allFinite() || !result.information.allFinite())
        throw std::overflow_error("Elliptical normal equations overflowed.");
    return result;
}
}

EllipticalGaussianFitResult fitEllipticalGaussian(const Image& dn, const Image& varianceDn2,
    PixelCoordinates origin, const PixelMask* mask, std::optional<EllipticalGaussian> initial,
    GaussianFitOptions options) {
    detail::validateOrigin(origin);
    if (dn.rows()<3 || dn.cols()<3 || dn.rows()!=varianceDn2.rows() || dn.cols()!=varianceDn2.cols() ||
        (mask && (mask->rows()!=dn.rows() || mask->cols()!=dn.cols())))
        throw std::invalid_argument("Fit needs matching image/variance/mask grids, at least 3x3.");
    if (options.maxIterations==0 || !std::isfinite(options.stepTolerance) || options.stepTolerance<=0 ||
        !std::isfinite(options.minimumHeightSnr) || options.minimumHeightSnr<0)
        throw std::invalid_argument("Invalid Gaussian fit options.");
    std::size_t usable = 0;
    for (Eigen::Index i = 0; i < dn.size(); ++i) {
        if (mask && mask->data()[i]) continue;
        if (!std::isfinite(dn.data()[i]) || !std::isfinite(varianceDn2.data()[i]) ||
            varianceDn2.data()[i]<=0 || !std::isfinite(1/varianceDn2.data()[i]))
            throw std::invalid_argument("Unmasked pixels need finite DN and positive finite variance.");
        ++usable;
    }
    if (usable<=7) throw std::invalid_argument("Seven-parameter fit needs more than seven usable pixels.");
    EllipticalGaussianFitResult result;
    result.degreesOfFreedom = usable-7;
    if (!initial) initial = initialize(dn,origin,mask);
    if (!initial) { result.status=FitStatus::NoSignal; return result; }
    detail::validate(*initial);
    if (!inWindow(*initial,dn,origin)) throw std::invalid_argument("Initial center must lie inside the fitting window.");
    Vector p = pack(*initial);
    auto normal = evaluate(dn,varianceDn2,mask,origin,p);
    double damping = 1e-3;
    for (std::size_t iteration = 0; iteration < options.maxIterations; ++iteration) {
        result.iterations = iteration+1;
        if ((normal.information.diagonal().array()<=0).any()) { result.status=FitStatus::Singular; break; }
        const Vector scale = normal.information.diagonal().array().sqrt().inverse();
        const Matrix scaled = scale.asDiagonal()*normal.information*scale.asDiagonal();
        const Vector gradient = scale.asDiagonal()*normal.gradient;
        if (gradient.norm()<=options.stepTolerance) { result.status=FitStatus::Converged; break; }
        bool accepted = false;
        for (int trial = 0; trial < 16; ++trial) {
            Matrix damped = scaled;
            damped.diagonal().array() += damping;
            const Vector step = scale.asDiagonal()*damped.ldlt().solve(gradient);
            const Vector candidateP = p+step;
            if (candidateP.allFinite() && std::isfinite(std::exp(candidateP[2])) && std::exp(candidateP[2])>0 &&
                precision(candidateP).allFinite() && precision(candidateP).determinant()>0) {
                const auto candidate = unpack(candidateP);
                if (candidate.sigmaMinorPixels>=0.2 && candidate.sigmaMajorPixels<=2*std::max(dn.rows(),dn.cols()) &&
                    inWindow(candidate,dn,origin)) {
                    const auto next = evaluate(dn,varianceDn2,mask,origin,candidateP);
                    if (next.cost<=normal.cost) {
                        const double backgroundStep = std::abs(step[6])/(1+std::abs(p[6]));
                        const double shearStep = std::abs(step[4])/(1+std::abs(p[4]));
                        p=candidateP; normal=next; accepted=true;
                        damping=std::max(1e-12,damping/3);
                        if (step.head<4>().cwiseAbs().maxCoeff()<=options.stepTolerance &&
                            std::abs(step[5])<=options.stepTolerance && shearStep<=options.stepTolerance &&
                            backgroundStep<=options.stepTolerance) result.status=FitStatus::Converged;
                        break;
                    }
                }
            }
            damping *= 10;
        }
        if (!accepted || result.status==FitStatus::Converged) break;
    }
    result.model=unpack(p); result.chiSquared=normal.cost;
    if (result.status!=FitStatus::Converged) return result;
    const Vector scale = normal.information.diagonal().array().sqrt().inverse();
    const Matrix scaled = scale.asDiagonal()*normal.information*scale.asDiagonal();
    const Eigen::SelfAdjointEigenSolver<Matrix> eigen(scaled,Eigen::EigenvaluesOnly);
    if (!scale.allFinite() || eigen.info()!=Eigen::Success || eigen.eigenvalues()[0]<=1e-10*eigen.eigenvalues()[6]) {
        result.status=FitStatus::Singular; return result;
    }
    const Matrix internalCov = scale.asDiagonal()*scaled.ldlt().solve(Matrix::Identity())*scale.asDiagonal();
    // Convert to [sample,line,height,Qxx,Qxy,Qyy,background], valid at circularity.
    Matrix transform = Matrix::Identity();
    transform(2,2)=result.model.heightDn;
    transform.block<3,3>(3,3).setZero();
    const double a=std::exp(p[3]), b=p[4], c=std::exp(p[5]);
    transform(3,3)=2*a*a; transform(4,3)=a*b; transform(4,4)=a;
    transform(5,4)=2*b; transform(5,5)=2*c*c;
    Matrix covariance = transform*internalCov*transform.transpose();
    covariance=(0.5*covariance+0.5*covariance.transpose()).eval();
    if (!covariance.allFinite() || (covariance.diagonal().array()<=0).any()) {
        result.status=FitStatus::Singular; return result;
    }
    if (result.model.heightDn/std::sqrt(covariance(2,2))<options.minimumHeightSnr) {
        result.status=FitStatus::NoSignal; return result;
    }
    result.measurement = PixelMeasurement{result.model.center,covariance.topLeftCorner<2,2>()};
    const double major=result.model.sigmaMajorPixels, minor=result.model.sigmaMinorPixels;
    const double gap=1/(major*major)-1/(minor*minor);
    if (std::abs(gap)>1e-6/(minor*minor)) {
        const double cs=std::cos(result.model.angleRadians), sn=std::sin(result.model.angleRadians);
        Matrix physical = Matrix::Identity();
        physical.block<3,3>(3,3).setZero();
        physical.block<1,3>(3,3) = -0.5*major*major*major*Eigen::RowVector3d{cs*cs,2*cs*sn,sn*sn};
        physical.block<1,3>(4,3) = -0.5*minor*minor*minor*Eigen::RowVector3d{sn*sn,-2*cs*sn,cs*cs};
        physical.block<1,3>(5,3) = Eigen::RowVector3d{-cs*sn,cs*cs-sn*sn,cs*sn}/gap;
        Matrix full=physical*covariance*physical.transpose();
        full=(0.5*full+0.5*full.transpose()).eval();
        if (full.allFinite()) result.parameterCovariance=full;
    }
    return result;
}

EllipticalGaussianFitResult fitBrightestEllipticalGaussian(const Image& dn, const Image& varianceDn2,
    int windowRadius, PixelCoordinates origin, const PixelMask* mask, GaussianFitOptions options) {
    detail::validateOrigin(origin);
    if (windowRadius<1 || dn.rows()<3 || dn.cols()<3 || dn.rows()!=varianceDn2.rows() || dn.cols()!=varianceDn2.cols() ||
        (mask && (mask->rows()!=dn.rows() || mask->cols()!=dn.cols())))
        throw std::invalid_argument("Source search needs matching grids, at least 3x3, and positive radius.");
    double peak=-std::numeric_limits<double>::infinity();
    Eigen::Index sample=0,line=0;
    for (Eigen::Index y=0;y<dn.rows();++y)
        for (Eigen::Index x=0;x<dn.cols();++x) {
            if (mask && (*mask)(y,x)) continue;
            if (!std::isfinite(dn(y,x)) || !std::isfinite(varianceDn2(y,x)) || varianceDn2(y,x)<=0 ||
                !std::isfinite(1/varianceDn2(y,x)))
                throw std::invalid_argument("Unmasked pixels need finite DN and positive finite variance.");
            if (dn(y,x)>peak) { peak=dn(y,x); sample=x; line=y; }
        }
    if (!std::isfinite(peak)) throw std::invalid_argument("Source search has no usable pixels.");
    const Eigen::Index width=std::min(2*Eigen::Index(windowRadius)+1,dn.cols());
    const Eigen::Index height=std::min(2*Eigen::Index(windowRadius)+1,dn.rows());
    const Eigen::Index x0=std::clamp(sample-windowRadius,Eigen::Index{0},dn.cols()-width);
    const Eigen::Index y0=std::clamp(line-windowRadius,Eigen::Index{0},dn.rows()-height);
    const Image window=dn.block(y0,x0,height,width), variance=varianceDn2.block(y0,x0,height,width);
    PixelMask localMask;
    if (mask) localMask=mask->block(y0,x0,height,width);
    return fitEllipticalGaussian(window,variance,{origin.sample+x0,origin.line+y0},mask ? &localMask : nullptr,
        std::nullopt,options);
}
} // namespace fd::opnav::image
