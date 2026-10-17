#include "opnav/estimation/EstimatorRejection.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

using namespace fd::opnav;
void require(bool value,const char* message) {if(!value)throw std::runtime_error(message);}
bool has(const RejectionDecision& d,RejectionReason reason) {
    return std::find(d.reasons.begin(),d.reasons.end(),reason)!=d.reasons.end();
}

int main() {try {
    const BodyRenderer renderer;
    CameraCalibration calibration;calibration.intrinsics={2003.44,{31,31},-Eigen::Matrix2d::Identity()/.012};
    OpticalImage image;image.body="ENCELADUS";image.time=TdbEpoch(100);
    const BodyParameters truth{{31.2,30.8},6,1.2,.5,800,105};
    image.dn=renderer.render({64,64},{},truth,.8).dn;
    image.mask=image::PixelMask::Zero(64,64);
    const auto fit=BodyEstimator(renderer,calibration).fit(image);
    require(fit.passesFitChecks(),"Test fixture did not converge.");
    RejectionContext context;context.calibration=&calibration;context.renderer=&renderer;
    const EstimatorRejection rejection;
    require(rejection.evaluate(image,fit,context).accepted(),"Good source rejected.");
    require(!RejectionDecision{}.accepted(),"Unevaluated decision accepted.");
    auto bad=fit;bad.status=image::FitStatus::NonConverged;
    require(has(rejection.evaluate(image,bad,context),RejectionReason::FitFailed),"Failed optimizer accepted.");
    bad=fit;bad.reducedChiSquared=4;
    require(has(rejection.evaluate(image,bad,context),RejectionReason::PoorResiduals),"Poor residuals accepted.");
    bad=fit;bad.measurement->covariance(0,0)=-1;
    require(has(rejection.evaluate(image,bad,context),RejectionReason::InvalidCovariance),"Negative covariance accepted.");
    bad=fit;bad.formalCenterCovariance=Eigen::Matrix2d::Identity()*.26*.26;
    require(has(rejection.evaluate(image,bad,context),RejectionReason::ExcessiveCenterUncertainty),"Uncertain extraction accepted.");
    bad=fit;bad.varianceDn2(0,0)=-1;
    require(has(rejection.evaluate(image,bad,context),RejectionReason::InvalidInput),"Invalid variance accepted.");
    bad=fit;bad.time=TdbEpoch(101);
    require(has(rejection.evaluate(image,bad,context),RejectionReason::InvalidInput),"Stale fit accepted.");
    require(has(rejection.evaluate(image,fit),RejectionReason::MissingRequiredCheck),"Missing calibration/shape ignored.");

    auto masked=image;masked.mask.block(26,26,11,11).setOnes();
    require(has(rejection.evaluate(masked,fit,context),RejectionReason::ExcessiveMasking),"Masked source accepted.");
    // Build consistent diagnostic arrays for weak and edge-clipped sources.
    const auto fixture=[&](BodyParameters parameters,PixelCoordinates origin) {
        OpticalImage observed=image;observed.dn=renderer.render({64,64},{},parameters,.8).dn;
        BodyFitResult result=fit;result.parameters=parameters;result.measurement->center=parameters.center;
        result.origin=origin;result.modelDn=renderer.render({37,37},origin,parameters,.8).dn;
        result.residualDn=image::Image::Zero(37,37);
        result.varianceDn2=(result.modelDn.array()-100).max(0)/10+.09+1./12;
        return std::pair(observed,result);
    };
    auto weakParameters=truth;weakParameters.amplitudeDn=1;
    auto [weak,weakFit]=fixture(weakParameters,{13,13});
    require(has(rejection.evaluate(weak,weakFit,context),RejectionReason::InsufficientSignal),"Weak source accepted.");
    auto edgeParameters=truth;edgeParameters.center={4,31};
    auto [edge,edgeFit]=fixture(edgeParameters,{0,13});
    require(has(rejection.evaluate(edge,edgeFit,context),RejectionReason::InsufficientBoundaryMargin),"Clipped body accepted.");

    RejectionOptions required;required.requireNavigationCheck=true;
    require(has(EstimatorRejection(required).evaluate(image,fit,context),RejectionReason::MissingRequiredCheck),"Required navigation silently skipped.");
    NavigationPrediction prediction;prediction.body=image.body;prediction.time=image.time;
    prediction.center=fit.parameters.center;prediction.covariancePixels2=Eigen::Matrix2d::Identity()*.09;
    context.prediction=prediction;
    require(rejection.evaluate(image,fit,context).accepted(),"Consistent prediction rejected.");
    context.prediction->center.sample+=10;
    require(has(rejection.evaluate(image,fit,context),RejectionReason::InconsistentWithNavigation),"Inconsistent prediction accepted.");
    context.prediction=prediction;context.prediction->center.sample+=.5;
    context.prediction->measurementCrossCovariance=Eigen::Matrix2d::Identity()*.03;
    auto correlated=fit;correlated.measurement->covariance=Eigen::Matrix2d::Identity()*.04;
    auto decision=rejection.evaluate(image,correlated,context);
    require(decision.innovationSquared && std::abs(*decision.innovationSquared-.25/.07)<1e-10,"Innovation cross covariance incorrect.");
    context.prediction->measurementCrossCovariance=-Eigen::Matrix2d::Identity()*10;
    require(has(rejection.evaluate(image,correlated,context),RejectionReason::InvalidCovariance),"Impossible cross covariance accepted.");
    context.prediction=prediction;context.prediction->time=TdbEpoch(101);
    require(has(rejection.evaluate(image,fit,context),RejectionReason::InvalidInput),"Stale prediction accepted.");
    context.prediction.reset();

    RejectionOptions ambiguityOptions;ambiguityOptions.checkSignal=false;ambiguityOptions.requireAmbiguityCheck=true;
    const EstimatorRejection ambiguity(ambiguityOptions);
    require(has(ambiguity.evaluate(weak,weakFit,context),RejectionReason::MissingRequiredCheck),"Missing ambiguity test ignored.");
    auto alternate=weakFit;alternate.parameters.center.sample+=.5;alternate.measurement->center=alternate.parameters.center;
    // Different reported windows/chi2 do not affect the common-pixel comparison.
    alternate.origin={0,0};alternate.reducedChiSquared=100;
    const std::array<BodyFitResult,1> alternatives{alternate};context.alternatives=alternatives;
    decision=ambiguity.evaluate(weak,weakFit,context);
    require(has(decision,RejectionReason::AmbiguousSolution) && decision.alternativesCompared==1,"Equal-quality distinct center accepted.");
    auto poor=fit;poor.parameters.center.sample+=3;poor.measurement->center=poor.parameters.center;poor.reducedChiSquared=0;
    const std::array<BodyFitResult,1> poorAlternatives{poor};context.alternatives=poorAlternatives;
    require(ambiguity.evaluate(image,fit,context).accepted(),"Poor alternative trusted its own chi2 instead of common pixels.");
    bool invalid=false;try{RejectionOptions o;o.maximumMaskedFraction=2;(void)EstimatorRejection(o);}catch(const std::invalid_argument&){invalid=true;}
    require(invalid,"Invalid options accepted.");
    std::cout<<"PASS: quality, masking, boundary, covariance, common-pixel ambiguity, navigation and missing-input gates.\n";
} catch(const std::exception& error){std::cerr<<"FAIL: "<<error.what()<<'\n';return 1;} }
