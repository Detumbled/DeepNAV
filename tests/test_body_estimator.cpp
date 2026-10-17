#include "opnav/estimation/BodyEstimator.hpp"
#include <Eigen/Geometry>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <random>
#include <sstream>

using namespace fd::opnav;
void require(bool condition,const char* message){if(!condition)throw std::runtime_error(message);}
template<class F>void rejects(F f){bool rejected=false;try{f();}catch(const std::exception&){rejected=true;}require(rejected,"Invalid fit input accepted.");}
int main(){try{
    const BodyRenderer renderer({},16);const BodyParameters truth{{20.23,19.74},7.15,1.19,.35,850,105};
    const auto evaluated=renderer.render({41,41},{},truth,.8,{},1,true);
    const auto p=truth.vector();
    for(int j=0;j<7;++j){auto plus=p,minus=p;const double step=j==5 ? .01 : 1e-5;plus[j]+=step;minus[j]-=step;
        const auto a=renderer.render({41,41},{},BodyParameters::fromVector(plus),.8);
        const auto b=renderer.render({41,41},{},BodyParameters::fromVector(minus),.8);
        const image::Image numeric=(a.dn-b.dn)/(2*step);
        require((Eigen::Map<const Eigen::VectorXd>(numeric.data(),numeric.size())-evaluated.jacobian.col(j)).norm()/std::max(1.,numeric.norm())<2e-5,"Body renderer analytic derivative failed.");}
    CameraCalibration calibration;calibration.intrinsics={2003.44,{20.23,19.74},-Eigen::Matrix2d::Identity()/.012};
    OpticalImage image;image.body="ENCELADUS";image.dn=evaluated.dn;image.mask=image::PixelMask::Zero(41,41);
    image.varianceDn2=image::Image::Ones(41,41);BodyFitOptions options;options.inferVarianceFromDn=false;
    const BodyEstimator estimator(renderer,calibration,options);
    const auto fit=estimator.fit(image);
    require(fit.passesFitChecks(),"Matched sphere fit did not converge.");
    require(std::hypot(fit.parameters.center.sample-truth.center.sample,fit.parameters.center.line-truth.center.line)<1e-4,"Matched-model center recovery failed.");
    BodyModel ellipsoid;ellipsoid.axes={1.08,.97,.94};ellipsoid.bodyToView=Eigen::AngleAxisd(.6,Eigen::Vector3d::UnitZ()).toRotationMatrix();
    const BodyRenderer ellipseRenderer(ellipsoid,16);image.dn=ellipseRenderer.render({41,41},{},truth,.8).dn;
    const auto ellipseFit=BodyEstimator(ellipseRenderer,calibration,options).fit(image);
    require(ellipseFit.passesFitChecks() && std::hypot(ellipseFit.parameters.center.sample-truth.center.sample,ellipseFit.parameters.center.line-truth.center.line)<1e-3,"Ellipsoid recovery failed.");
    // Exposure integration and a rotated non-spherical surface must retain correct derivatives.
    const Eigen::Vector2d smear{1.8,.8};
    const auto smeared=ellipseRenderer.render({41,41},{},truth,.8,smear,9,true);
    for(int j=0;j<7;++j){auto plus=p,minus=p;const double step=j==5 ? .01 : 1e-5;plus[j]+=step;minus[j]-=step;
        const image::Image numeric=(ellipseRenderer.render({41,41},{},BodyParameters::fromVector(plus),.8,smear,9).dn-
            ellipseRenderer.render({41,41},{},BodyParameters::fromVector(minus),.8,smear,9).dn)/(2*step);
        require((Eigen::Map<const Eigen::VectorXd>(numeric.data(),numeric.size())-smeared.jacobian.col(j)).norm()/std::max(1.,numeric.norm())<2e-5,"Integrated ellipsoid derivative failed.");}
    image.dn=smeared.dn;image.smearPixels=smear;options.exposureSamples=9;
    const auto integrated=BodyEstimator(ellipseRenderer,calibration,options).fit(image);
    require(integrated.passesFitChecks() && std::hypot(integrated.parameters.center.sample-truth.center.sample,integrated.parameters.center.line-truth.center.line)<1e-3,"Exposure-integrated recovery failed.");
    options.exposureSamples=1;image.smearPixels.setZero();
    // Nominal camera/attitude budgets must enlarge measurement covariance, not fit covariance.
    image.dn=evaluated.dn;image.attitudeCovariance=Eigen::Matrix3d::Identity()*1e-12;
    calibration.principalPointCovariance=Eigen::Matrix2d::Identity()*.01;
    const auto uncertain=BodyEstimator(renderer,calibration,options).fit(image);
    require(uncertain.measurement && (uncertain.measurement->covariance-uncertain.formalCenterCovariance).trace()>.02,"Nuisance covariance was not propagated.");
    image.body="TITAN";rejects([&]{(void)estimator.fit(image);});image.body="ENCELADUS";
    image.mask=image::PixelMask::Ones(41,41);require(estimator.fit(image).status==image::FitStatus::NoSignal,"All-masked frame must fail.");
    image.mask=image::PixelMask::Zero(41,41);image.dn.setConstant(105);require(estimator.fit(image).status==image::FitStatus::NoSignal,"Blank frame must fail.");
    OpenCVDistortionCoefficients distortion;distortion.K1=100;distortion.P1=.01;
    const CameraModel camera(calibration.intrinsics,std::make_shared<OpenCVDistortionModel>(distortion));
    for(const Eigen::Vector3d v:{Eigen::Vector3d(.001,-.002,1),Eigen::Vector3d(-.003,.001,1)})
        require((camera.unproject(camera.project(v))-v.normalized()).norm()<1e-12,"Camera inverse distortion roundtrip failed.");
    // Independent prior C++ ray-traced frame / Python quadrature regression, when local data exist.
    const auto path=std::filesystem::path(DEEPNAV_SOURCE_DIR)/"results/cassini_synthetic_optical/simulated_sequence/frame_0.csv";
    if(std::filesystem::exists(path)){
        std::ifstream input(path);std::string line;std::getline(input,line);image.dn=image::Image::Zero(128,128);image.mask=image::PixelMask::Zero(128,128);
        while(std::getline(input,line)){std::stringstream row(line);std::vector<double> cells;std::string cell;while(std::getline(row,cell,','))cells.push_back(std::stod(cell));
            image.dn(int(cells[1]),int(cells[0]))=cells[2];image.mask(int(cells[1]),int(cells[0]))=int(cells[4]);}
        BodyFitOptions inferred;calibration.intrinsics.principalPoint={63.23,62.74};
        const auto actual=BodyEstimator(BodyRenderer({},24),calibration,inferred).fit(image);
        std::cout<<"Ray-traced fixture center: "<<actual.parameters.center.sample<<", "<<actual.parameters.center.line<<" chi2="<<actual.reducedChiSquared<<"\n";
        require(actual.passesFitChecks() && std::hypot(actual.parameters.center.sample-63.23,actual.parameters.center.line-62.74)<.08,"Independent ray-traced sphere fit failed.");
    }
    std::cout<<"PASS: analytic derivatives, sphere/ellipsoid recovery, uncertainty propagation, invalid images, inverse distortion.\n";
}catch(const std::exception& e){std::cerr<<"FAIL: "<<e.what()<<'\n';return 1;}}
