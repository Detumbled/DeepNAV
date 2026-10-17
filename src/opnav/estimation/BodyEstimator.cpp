#include "opnav/estimation/BodyEstimator.hpp"
#include <Eigen/Cholesky>
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <cmath>
#include <numbers>
#include <stdexcept>
#include <vector>

namespace fd::opnav {
namespace {
using Parameters=Eigen::Matrix<double,7,1>;
double median(std::vector<double> v) {
    if(v.empty()) throw std::invalid_argument("No valid pixels.");
    const auto m=v.begin()+v.size()/2;std::nth_element(v.begin(),m,v.end());return *m;
}
void require(bool condition,const char* message) {if(!condition)throw std::invalid_argument(message);}
template<class Matrix> bool covariance(const Matrix& p) {
    return p.allFinite() && (p-p.transpose()).norm()<1e-10 &&
        Eigen::SelfAdjointEigenSolver<Matrix>(p).eigenvalues().minCoeff()>=-1e-12;
}
Parameters bounded(Parameters p,int x0,int y0,int width,int height) {
    p[0]=std::clamp(p[0],double(x0),double(x0+width-1));p[1]=std::clamp(p[1],double(y0),double(y0+height-1));
    p[2]=std::clamp(p[2],.5,30.);p[3]=std::clamp(p[3],.04,std::numbers::pi-.04);
    p[5]=std::clamp(p[5],1e-6,1e6);p[6]=std::clamp(p[6],0.,1e6);return p;
}
}
BodyEstimator::BodyEstimator(BodyRenderer renderer,CameraCalibration calibration,BodyFitOptions options)
    :renderer_(std::move(renderer)),calibration_(std::move(calibration)),options_(options) {
    require(calibration_.intrinsics.focalLength>0 && std::isfinite(calibration_.intrinsics.focalLength) &&
            calibration_.intrinsics.focalPlaneToPixelMatrix_PixelPermm.allFinite() &&
            std::abs(calibration_.intrinsics.focalPlaneToPixelMatrix_PixelPermm.determinant())>1e-15 &&
            std::isfinite(calibration_.intrinsics.principalPoint.sample) && std::isfinite(calibration_.intrinsics.principalPoint.line),"Invalid camera intrinsics.");
    const Eigen::Matrix2d k=calibration_.intrinsics.focalPlaneToPixelMatrix_PixelPermm;
    require((k*k.transpose()/std::abs(k.determinant())-Eigen::Matrix2d::Identity()).norm()<1e-10,"Body estimator currently requires square pixels without shear.");
    require(calibration_.psfSigmaPixels>=.05 && calibration_.psfSigmaPixels<=10 && std::isfinite(calibration_.psfSigmaPixels) &&
            calibration_.electronsPerDn>0 && std::isfinite(calibration_.electronsPerDn) &&
            calibration_.readNoiseElectrons>=0 && std::isfinite(calibration_.readNoiseElectrons) &&
            std::isfinite(calibration_.biasDn) && calibration_.saturationDn>calibration_.biasDn &&
            std::isfinite(calibration_.saturationDn) && covariance(calibration_.principalPointCovariance) &&
            calibration_.relativeFocalSigma>=0 && std::isfinite(calibration_.relativeFocalSigma) &&
            calibration_.centerModelSigmaPixels>=0 && std::isfinite(calibration_.centerModelSigmaPixels),"Invalid detector/calibration uncertainty.");
    require(options_.windowRadius>=4 && options_.windowRadius<=40 && options_.maxIterations>0 &&
            options_.exposureSamples>=1 && options_.exposureSamples<=101 && options_.maximumReducedChiSquared>0 &&
            std::isfinite(options_.maximumReducedChiSquared),"Invalid body fit options.");
}
BodyFitResult BodyEstimator::fit(const OpticalImage& image,std::optional<BodyParameters> initialGuess) const {
    require(image.body==renderer_.model().body,"Observation body does not match registered renderer.");
    require(std::isfinite(image.time.secondsPastJ2000) && std::isfinite(image.exposureSeconds) && image.exposureSeconds>=0 &&
            image.smearPixels.allFinite() && covariance(image.attitudeCovariance),"Invalid observation metadata.");
    const int rows=image.dn.rows(),cols=image.dn.cols();
    require(rows>0 && cols>0,"Empty body image.");
    require(image.mask.size()==0 || (image.mask.rows()==rows && image.mask.cols()==cols),"Mask dimensions do not match image.");
    if(!options_.inferVarianceFromDn) require(image.varianceDn2.rows()==rows && image.varianceDn2.cols()==cols,"Supplied variance dimensions do not match image.");
    if(initialGuess) require(initialGuess->vector().allFinite() && initialGuess->radiusPixels>0 && initialGuess->amplitudeDn>0,"Invalid initial guess.");
    BodyFitResult answer;answer.body=image.body;answer.time=image.time;
    const BodyRenderer renderer=image.bodyToView ? renderer_.reoriented(*image.bodyToView) : renderer_;
    int peakX=0,peakY=0;double peak=-1e100;
    const auto validPixel=[&](int y,int x) {return std::isfinite(image.dn(y,x)) && image.dn(y,x)>=0 && image.dn(y,x)<calibration_.saturationDn &&
        (image.mask.size()==0 || image.mask(y,x)==0);};
    for(int y=0;y<rows;++y)for(int x=0;x<cols;++x)if(validPixel(y,x) && image.dn(y,x)>peak) {peak=image.dn(y,x);peakX=x;peakY=y;}
    if(peak<0) {answer.status=image::FitStatus::NoSignal;return answer;}
    if(initialGuess) {peakX=int(std::round(initialGuess->center.sample));peakY=int(std::round(initialGuess->center.line));}
    const int x0=std::max(0,peakX-options_.windowRadius),y0=std::max(0,peakY-options_.windowRadius);
    const int width=std::min(cols-1,peakX+options_.windowRadius)-x0+1,height=std::min(rows-1,peakY+options_.windowRadius)-y0+1;
    require(width>=5 && height>=5,"Initial guess does not intersect a usable image window.");
    answer.origin={double(x0),double(y0)};
    const image::Image observed=image.dn.block(y0,x0,height,width);
    std::vector<double> border,noise,validDn;std::vector<int> usable;
    image::Image variance(height,width);
    for(int y=0;y<height;++y)for(int x=0;x<width;++x) {
        const double dn=observed(y,x);
        double v=options_.inferVarianceFromDn ? std::max(0.,dn-calibration_.biasDn)/calibration_.electronsPerDn+
            std::pow(calibration_.readNoiseElectrons/calibration_.electronsPerDn,2)+(calibration_.digitized ? 1./12 : 0) : image.varianceDn2(y+y0,x+x0);
        variance(y,x)=std::max(v,1e-8);
        if(!validPixel(y+y0,x+x0) || !std::isfinite(v) || v<0) continue;
        usable.push_back(y*width+x);noise.push_back(std::sqrt(variance(y,x)));validDn.push_back(dn);
        if(y==0 || x==0 || y==height-1 || x==width-1) border.push_back(dn);
    }
    if(usable.size()<30) {answer.status=image::FitStatus::NoSignal;return answer;}
    const double bg=border.empty() ? median(validDn) : median(border),sigma=median(noise);
    if(peak-bg<5*sigma) {answer.status=image::FitStatus::NoSignal;return answer;}
    Eigen::Vector2d moment=Eigen::Vector2d::Zero();double signalSum=0;
    int minX=cols,minY=rows,maxX=0,maxY=0;
    for(int index:usable) {const int y=index/width,x=index%width;const double s=std::max(0.,observed(y,x)-bg-3*sigma);
        moment+=s*Eigen::Vector2d(x+x0,y+y0);signalSum+=s;if(s>0){minX=std::min(minX,x+x0);maxX=std::max(maxX,x+x0);minY=std::min(minY,y+y0);maxY=std::max(maxY,y+y0);}}
    if(signalSum<=0) {answer.status=image::FitStatus::NoSignal;return answer;}
    moment/=signalSum;Eigen::Matrix2d moments=Eigen::Matrix2d::Zero();
    for(int index:usable) {const int y=index/width,x=index%width;const double s=std::max(0.,observed(y,x)-bg-3*sigma);
        const Eigen::Vector2d d=Eigen::Vector2d(x+x0,y+y0)-moment;moments+=s*d*d.transpose();}
    const Eigen::Vector2d minor=Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d>(moments).eigenvectors().col(0);
    const double angle=std::atan2(minor.y(),minor.x()),radius=std::clamp(.5*std::max(maxX-minX,maxY-minY),3.,12.);
    Parameters scales;scales<<1,1,1,.5,.5,std::max(1.,peak-bg),10;
    const Eigen::Vector2d smear=options_.exposureSamples>1 ? image.smearPixels : Eigen::Vector2d::Zero();
    const auto evaluate=[&](const Parameters& p,bool derivatives,Eigen::VectorXd& residual,Eigen::MatrixXd& jac) {
        const auto model=renderer.render({width,height},answer.origin,BodyParameters::fromVector(p),calibration_.psfSigmaPixels,smear,options_.exposureSamples,derivatives);
        residual.resize(usable.size());if(derivatives)jac.resize(usable.size(),7);
        for(std::size_t i=0;i<usable.size();++i) {const int k=usable[i],y=k/width,x=k%width;const double weight=1/std::sqrt(variance(y,x));
            residual[i]=(model.dn(y,x)-observed(y,x))*weight;
            if(derivatives)jac.row(i)=model.jacobian.row(k)*weight*scales.asDiagonal();}
    };
    double bestCost=1e300;Parameters best;bool bestConverged=false;int bestIterations=0;
    for(int start=0;start<(initialGuess ? 1 : 2);++start) {
        Parameters p;
        if(initialGuess)p=initialGuess->vector();else {
            const double azimuth=angle+start*std::numbers::pi;
            p<<moment.x()-.35*radius*std::cos(azimuth),moment.y()-.35*radius*std::sin(azimuth),radius,65*std::numbers::pi/180,azimuth,1.2*(peak-bg),bg;
        }
        p=bounded(p,x0,y0,width,height);double lambda=1e-3;bool converged=false;int iteration=0;
        Eigen::VectorXd residual;Eigen::MatrixXd jac;evaluate(p,true,residual,jac);double cost=residual.squaredNorm();
        for(;iteration<options_.maxIterations;++iteration) {
            const Eigen::Matrix<double,7,7> information=jac.transpose()*jac;
            Eigen::Matrix<double,7,7> damped=information;
            damped.diagonal().array()+=lambda*information.diagonal().cwiseMax(1e-8).array();
            const Parameters step=damped.ldlt().solve(-jac.transpose()*residual);
            if(!step.allFinite())break;
            const Parameters trial=bounded(p+scales.cwiseProduct(step),x0,y0,width,height);
            Eigen::VectorXd next;Eigen::MatrixXd unused;evaluate(trial,false,next,unused);const double nextCost=next.squaredNorm();
            if(nextCost<cost) {
                const double improvement=(cost-nextCost)/std::max(1.,cost);p=trial;cost=nextCost;lambda=std::max(1e-12,lambda/3);
                evaluate(p,true,residual,jac);
                if(step.norm()<1e-5 || improvement<1e-9){converged=true;break;}
            } else {lambda*=10;if(step.norm()<1e-6){converged=true;break;}if(lambda>1e15)break;}
        }
        if(cost<bestCost) {bestCost=cost;best=p;bestConverged=converged;bestIterations=iteration+1;}
    }
    answer.parameters=BodyParameters::fromVector(best);answer.iterations=bestIterations;
    answer.reducedChiSquared=bestCost/(usable.size()-7);answer.modelMismatch=answer.reducedChiSquared>options_.maximumReducedChiSquared;
    answer.status=bestConverged ? image::FitStatus::Converged : image::FitStatus::NonConverged;
    const auto final=renderer.render({width,height},answer.origin,answer.parameters,calibration_.psfSigmaPixels,smear,options_.exposureSamples);
    answer.modelDn=final.dn;answer.residualDn=observed-final.dn;answer.varianceDn2=variance;
    Eigen::VectorXd residual;Eigen::MatrixXd jac;evaluate(best,true,residual,jac);
    const Eigen::Matrix<double,7,7> information=jac.transpose()*jac;
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double,7,7>> decomposition(information);
    if(decomposition.eigenvalues().minCoeff()<decomposition.eigenvalues().maxCoeff()*1e-12) {answer.status=image::FitStatus::Singular;return answer;}
    const Eigen::Matrix<double,7,7> parameterCov=scales.asDiagonal()*information.ldlt().solve(Eigen::Matrix<double,7,7>::Identity())*scales.asDiagonal();
    answer.formalCenterCovariance=parameterCov.topLeftCorner<2,2>();
    Eigen::Matrix2d effective=answer.formalCenterCovariance+calibration_.principalPointCovariance+
        Eigen::Matrix2d::Identity()*std::pow(calibration_.centerModelSigmaPixels,2);
    const Eigen::Vector2d offset{best[0]-calibration_.intrinsics.principalPoint.sample,best[1]-calibration_.intrinsics.principalPoint.line};
    effective+=std::pow(calibration_.relativeFocalSigma,2)*offset*offset.transpose();
    const Eigen::Matrix2d focal=calibration_.intrinsics.focalLength*calibration_.intrinsics.focalPlaneToPixelMatrix_PixelPermm;
    const Eigen::Vector2d plane=focal.inverse()*offset;const Eigen::Vector3d ray{plane.x(),plane.y(),1};
    Eigen::Matrix3d cross;cross<<0,-1,ray.y(),1,0,-ray.x(),-ray.y(),ray.x(),0;
    Eigen::Matrix<double,2,3> projection;projection<<1,0,-ray.x(),0,1,-ray.y();
    const Eigen::Matrix<double,2,3> attitudeJac=-focal*projection*cross;
    effective+=attitudeJac*image.attitudeCovariance*attitudeJac.transpose();
    if(bestConverged) answer.measurement=image::PixelMeasurement{answer.parameters.center,effective};
    return answer;
}
} // namespace fd::opnav
