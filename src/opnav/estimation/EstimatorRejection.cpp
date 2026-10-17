#include "opnav/estimation/EstimatorRejection.hpp"
#include <Eigen/Cholesky>
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace fd::opnav {
namespace {
bool covariance(const Eigen::Matrix2d& p,bool positive) {
    if(!p.allFinite() || (p-p.transpose()).norm()>1e-10*std::max(1.,p.norm())) return false;
    const double minimum=Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d>(p).eigenvalues().minCoeff();
    return positive ? minimum>0 : minimum>=-1e-12*std::max(1.,p.norm());
}
bool sameObservation(const OpticalImage& image,const BodyFitResult& fit) {
    return fit.body==image.body && std::isfinite(image.time.secondsPastJ2000) &&
        std::isfinite(fit.time.secondsPastJ2000) && std::abs(fit.time.secondsPastJ2000-image.time.secondsPastJ2000)<=1e-6;
}
double sigma(const Eigen::Matrix2d& p) {
    return std::sqrt(Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d>(p).eigenvalues().maxCoeff());
}
}

const char* rejectionReasonName(RejectionReason reason) {
    switch(reason) {
        case RejectionReason::FitFailed:return "fit_failed";
        case RejectionReason::InvalidInput:return "invalid_input";
        case RejectionReason::InvalidCovariance:return "invalid_covariance";
        case RejectionReason::PoorResiduals:return "poor_residuals";
        case RejectionReason::InsufficientSignal:return "insufficient_signal";
        case RejectionReason::ExcessiveMasking:return "excessive_masking";
        case RejectionReason::InsufficientBoundaryMargin:return "insufficient_boundary_margin";
        case RejectionReason::ExcessiveCenterUncertainty:return "excessive_center_uncertainty";
        case RejectionReason::AmbiguousSolution:return "ambiguous_solution";
        case RejectionReason::InconsistentWithNavigation:return "inconsistent_with_navigation";
        case RejectionReason::MissingRequiredCheck:return "missing_required_check";
    }
    return "unknown";
}

EstimatorRejection::EstimatorRejection(RejectionOptions options):options_(options) {
    for(double positive:{options.maximumReducedChiSquared,options.minimumPeakSignalToNoise,options.minimumSourceSignalToNoise,
        options.psfSupportSigmas,options.maximumFormalCenterSigmaPixels,options.maximumMeasurementCenterSigmaPixels,
        options.ambiguityCenterSeparationPixels,options.maximumInnovationSquared})
        if(!std::isfinite(positive) || positive<=0) throw std::invalid_argument("Invalid rejection threshold.");
    if(!std::isfinite(options.maximumMaskedFraction) || options.maximumMaskedFraction<0 || options.maximumMaskedFraction>=1 ||
       options.minimumValidPixels<8 || !std::isfinite(options.minimumBoundaryMarginPixels) || options.minimumBoundaryMarginPixels<0 ||
       !std::isfinite(options.ambiguityDeltaChiSquared) || options.ambiguityDeltaChiSquared<0)
        throw std::invalid_argument("Invalid rejection options.");
}

RejectionDecision EstimatorRejection::evaluate(const OpticalImage& image,const BodyFitResult& fit,const RejectionContext& context) const {
    RejectionDecision result;result.evaluated=true;
    const auto state=[&](RejectionCheck check,CheckState value) {result.checks[static_cast<std::size_t>(check)]=value;};
    const auto reject=[&](RejectionReason reason,RejectionCheck check) {
        if(std::find(result.reasons.begin(),result.reasons.end(),reason)==result.reasons.end())result.reasons.push_back(reason);
        state(check,CheckState::Failed);
    };
    if(fit.status!=image::FitStatus::Converged || !fit.measurement) {
        reject(RejectionReason::FitFailed,RejectionCheck::FitValidity);return result;
    }
    const auto& p=fit.parameters;
    const int rows=fit.modelDn.rows(),cols=fit.modelDn.cols();
    const auto matched=[&](const image::Image& a) {return a.rows()==rows && a.cols()==cols;};
    if(!sameObservation(image,fit) || !p.vector().allFinite() || p.radiusPixels<=0 || p.amplitudeDn<=0 ||
       p.phaseRadians<0 || p.phaseRadians>std::acos(-1.) || image.dn.rows()<1 || image.dn.cols()<1 ||
       context.exposureSamples<1 || context.exposureSamples>101 || !image.smearPixels.allFinite() ||
       !std::isfinite(fit.origin.sample) || !std::isfinite(fit.origin.line) ||
       fit.origin.sample<0 || fit.origin.line<0 || fit.origin.sample+cols>image.dn.cols() || fit.origin.line+rows>image.dn.rows() ||
       fit.origin.sample!=std::floor(fit.origin.sample) || fit.origin.line!=std::floor(fit.origin.line) ||
       rows<1 || cols<1 || !matched(fit.residualDn) || !matched(fit.varianceDn2) ||
       (image.mask.size()!=0 && (image.mask.rows()!=image.dn.rows() || image.mask.cols()!=image.dn.cols())) ||
       !std::isfinite(fit.measurement->center.sample) || !std::isfinite(fit.measurement->center.line) ||
       std::hypot(fit.measurement->center.sample-p.center.sample,fit.measurement->center.line-p.center.line)>1e-8) {
        reject(RejectionReason::InvalidInput,RejectionCheck::FitValidity);return result;
    }
    state(RejectionCheck::FitValidity,CheckState::Passed);
    if(!covariance(fit.formalCenterCovariance,true) || !covariance(fit.measurement->covariance,true)) {
        reject(RejectionReason::InvalidCovariance,RejectionCheck::Uncertainty);return result;
    }
    if(!std::isfinite(fit.reducedChiSquared) || fit.reducedChiSquared<0 ||
       fit.reducedChiSquared>options_.maximumReducedChiSquared || fit.modelMismatch)
        reject(RejectionReason::PoorResiduals,RejectionCheck::Residuals);
    else state(RejectionCheck::Residuals,CheckState::Passed);

    if(options_.checkUncertainty) {
        result.formalCenterSigmaPixels=sigma(fit.formalCenterCovariance);
        result.measurementCenterSigmaPixels=sigma(fit.measurement->covariance);
        if(*result.formalCenterSigmaPixels>options_.maximumFormalCenterSigmaPixels ||
           *result.measurementCenterSigmaPixels>options_.maximumMeasurementCenterSigmaPixels)
            reject(RejectionReason::ExcessiveCenterUncertainty,RejectionCheck::Uncertainty);
        else state(RejectionCheck::Uncertainty,CheckState::Passed);
    }
    const CameraCalibration* calibration=context.calibration;
    if(calibration && (!std::isfinite(calibration->psfSigmaPixels) || calibration->psfSigmaPixels<=0 ||
        !std::isfinite(calibration->saturationDn))) {
        reject(RejectionReason::InvalidInput,RejectionCheck::FitValidity);return result;
    }
    if((options_.checkSignal || options_.checkMasking || options_.checkBoundary) && !calibration) {
        reject(RejectionReason::MissingRequiredCheck,RejectionCheck::FitValidity);return result;
    }
    if(context.renderer && context.renderer->model().body!=image.body) {
        reject(RejectionReason::InvalidInput,RejectionCheck::FitValidity);return result;
    }
    if(options_.checkBoundary && !context.renderer) reject(RejectionReason::MissingRequiredCheck,RejectionCheck::Boundary);
    double support=p.radiusPixels*(context.renderer?context.renderer->model().axes.maxCoeff():1.);
    if(calibration) support+=options_.psfSupportSigmas*calibration->psfSigmaPixels;
    if(context.exposureSamples>1) support+=image.smearPixels.norm()/2;
    if(!std::isfinite(support) || support<=0) {
        reject(RejectionReason::InvalidInput,RejectionCheck::FitValidity);return result;
    }
    if(options_.checkBoundary && context.renderer) {
        const double margin=std::min({p.center.sample+.5,image.dn.cols()-.5-p.center.sample,
            p.center.line+.5,image.dn.rows()-.5-p.center.line,p.center.sample-fit.origin.sample+.5,
            fit.origin.sample+cols-.5-p.center.sample,p.center.line-fit.origin.line+.5,
            fit.origin.line+rows-.5-p.center.line})-support;
        result.boundaryMarginPixels=margin;
        if(margin<options_.minimumBoundaryMarginPixels) reject(RejectionReason::InsufficientBoundaryMargin,RejectionCheck::Boundary);
        else state(RejectionCheck::Boundary,CheckState::Passed);
    }

    const int x0=int(fit.origin.sample),y0=int(fit.origin.line);
    const auto valid=[&](int y,int x) {
        const double observed=image.dn(y+y0,x+x0);
        return std::isfinite(observed) && observed>=0 && (!calibration || observed<calibration->saturationDn) &&
            (image.mask.size()==0 || image.mask(y+y0,x+x0)==0);
    };
    int usable=0,aperture=0,bad=0;double peak=0,numerator=0,information=0;
    for(int y=0;y<rows;++y)for(int x=0;x<cols;++x) {
        const bool inside=std::hypot(x+x0-p.center.sample,y+y0-p.center.line)<=support;
        if(inside) {++aperture;if(!valid(y,x))++bad;}
        if(!valid(y,x))continue;
        if(!std::isfinite(fit.modelDn(y,x)) || !std::isfinite(fit.residualDn(y,x)) ||
           !std::isfinite(fit.varianceDn2(y,x)) || fit.varianceDn2(y,x)<=0) {
            reject(RejectionReason::InvalidInput,RejectionCheck::FitValidity);return result;
        }
        const double expectedResidual=image.dn(y+y0,x+x0)-fit.modelDn(y,x);
        if(std::abs(expectedResidual-fit.residualDn(y,x))>1e-8*std::max(1.,std::abs(image.dn(y+y0,x+x0)))) {
            reject(RejectionReason::InvalidInput,RejectionCheck::FitValidity);return result;
        }
        ++usable;
        const double source=std::max(0.,fit.modelDn(y,x)-p.backgroundDn),observed=image.dn(y+y0,x+x0)-p.backgroundDn;
        const double variance=fit.varianceDn2(y,x);
        if(inside)peak=std::max(peak,observed/std::sqrt(variance));
        numerator+=source*observed/variance;information+=source*source/variance;
    }
    if(usable<options_.minimumValidPixels) reject(RejectionReason::InsufficientSignal,RejectionCheck::Signal);
    if(options_.checkMasking) {
        result.maskedFraction=aperture ? double(bad)/aperture : 1.;
        if(*result.maskedFraction>options_.maximumMaskedFraction) reject(RejectionReason::ExcessiveMasking,RejectionCheck::Masking);
        else state(RejectionCheck::Masking,CheckState::Passed);
    }
    if(options_.checkSignal) {
        result.peakSignalToNoise=peak;
        result.sourceSignalToNoise=information>0 ? numerator/std::sqrt(information) : 0.;
        if(!std::isfinite(*result.sourceSignalToNoise) || !std::isfinite(peak) ||
           usable<options_.minimumValidPixels || peak<options_.minimumPeakSignalToNoise ||
           *result.sourceSignalToNoise<options_.minimumSourceSignalToNoise)
            reject(RejectionReason::InsufficientSignal,RejectionCheck::Signal);
        else state(RejectionCheck::Signal,CheckState::Passed);
    }

    if(!context.alternatives.empty() && (!context.renderer || !calibration || context.exposureSamples<1 || context.exposureSamples>101))
        reject(RejectionReason::MissingRequiredCheck,RejectionCheck::Ambiguity);
    else if(!context.alternatives.empty()) {
        const auto renderer=image.bodyToView?context.renderer->reoriented(*image.bodyToView):*context.renderer;
        const auto cost=[&](const BodyParameters& candidate) {
            const auto rendered=renderer.render({cols,rows},fit.origin,candidate,calibration->psfSigmaPixels,
                context.exposureSamples>1?image.smearPixels:Eigen::Vector2d::Zero(),context.exposureSamples).dn;
            double value=0;
            for(int y=0;y<rows;++y)for(int x=0;x<cols;++x)if(valid(y,x))
                value+=std::pow(image.dn(y+y0,x+x0)-rendered(y,x),2)/fit.varianceDn2(y,x);
            return value;
        };
        const double reference=cost(p);
        for(const auto& alternative:context.alternatives) {
            if(!sameObservation(image,alternative) || !alternative.parameters.vector().allFinite() ||
                alternative.parameters.radiusPixels<=0 || alternative.parameters.amplitudeDn<0 ||
                alternative.parameters.phaseRadians<0 || alternative.parameters.phaseRadians>std::acos(-1.)) {
                reject(RejectionReason::InvalidInput,RejectionCheck::Ambiguity);continue;
            }
            if(alternative.status!=image::FitStatus::Converged || !alternative.measurement)continue;
            ++result.alternativesCompared;
            const double separation=std::hypot(alternative.parameters.center.sample-p.center.sample,alternative.parameters.center.line-p.center.line);
            if(separation>options_.ambiguityCenterSeparationPixels && cost(alternative.parameters)<=reference+options_.ambiguityDeltaChiSquared) {
                result.ambiguousCenterSeparationPixels=separation;
                reject(RejectionReason::AmbiguousSolution,RejectionCheck::Ambiguity);
            }
        }
        if(result.alternativesCompared>0 && result.checks[static_cast<std::size_t>(RejectionCheck::Ambiguity)]!=CheckState::Failed)
            state(RejectionCheck::Ambiguity,CheckState::Passed);
    }
    if(options_.requireAmbiguityCheck && result.checks[static_cast<std::size_t>(RejectionCheck::Ambiguity)]==CheckState::NotRun)
        reject(RejectionReason::MissingRequiredCheck,RejectionCheck::Ambiguity);

    if(context.prediction) {
        const auto& prediction=*context.prediction;
        if(prediction.body!=image.body || !std::isfinite(prediction.time.secondsPastJ2000) ||
           std::abs(prediction.time.secondsPastJ2000-image.time.secondsPastJ2000)>1e-6 ||
           !std::isfinite(prediction.center.sample) || !std::isfinite(prediction.center.line) ||
           !covariance(prediction.covariancePixels2,false) || !prediction.measurementCrossCovariance.allFinite())
            reject(RejectionReason::InvalidInput,RejectionCheck::Navigation);
        else {
            const Eigen::Matrix2d cross=prediction.measurementCrossCovariance;
            const Eigen::Matrix2d innovation=fit.measurement->covariance+prediction.covariancePixels2-cross-cross.transpose();
            Eigen::Matrix4d joint;joint<<fit.measurement->covariance,cross,cross.transpose(),prediction.covariancePixels2;
            const bool validJoint=Eigen::SelfAdjointEigenSolver<Eigen::Matrix4d>(joint).eigenvalues().minCoeff()>=-1e-12*std::max(1.,joint.norm());
            if(!validJoint || !covariance(innovation,true))reject(RejectionReason::InvalidCovariance,RejectionCheck::Navigation);
            else {
                const Eigen::Vector2d error{p.center.sample-prediction.center.sample,p.center.line-prediction.center.line};
                result.innovationSquared=error.dot(innovation.ldlt().solve(error));
                if(*result.innovationSquared>options_.maximumInnovationSquared)reject(RejectionReason::InconsistentWithNavigation,RejectionCheck::Navigation);
                else state(RejectionCheck::Navigation,CheckState::Passed);
            }
        }
    } else if(options_.requireNavigationCheck)reject(RejectionReason::MissingRequiredCheck,RejectionCheck::Navigation);
    return result;
}
} // namespace fd::opnav
