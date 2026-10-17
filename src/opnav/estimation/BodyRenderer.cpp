#include "opnav/estimation/BodyEstimator.hpp"
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <cmath>
#include <numbers>
#include <stdexcept>

namespace fd::opnav {
Eigen::Matrix<double,7,1> BodyParameters::vector() const {
    Eigen::Matrix<double,7,1> v;v<<center.sample,center.line,radiusPixels,phaseRadians,sunAngleRadians,amplitudeDn,backgroundDn;return v;
}
BodyParameters BodyParameters::fromVector(const Eigen::Matrix<double,7,1>& v) {
    return {{v[0],v[1]},v[2],v[3],v[4],v[5],v[6]};
}
namespace {
std::pair<int,int> footprint(double position,double origin,double sigma,int pixels) {
    if(!std::isfinite(position)) throw std::invalid_argument("Nonfinite projected body quadrature.");
    return {int(std::clamp(std::ceil(position-origin-8*sigma-.5),0.,double(pixels))),
            int(std::clamp(std::floor(position-origin+8*sigma+.5),-1.,double(pixels-1)))};
}
std::pair<Eigen::VectorXd,Eigen::VectorXd> gaussLegendre(int n) {
    Eigen::MatrixXd jac=Eigen::MatrixXd::Zero(n,n);
    for(int i=1;i<n;++i) jac(i,i-1)=jac(i-1,i)=i/std::sqrt(4.*i*i-1);
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(jac);
    return {solver.eigenvalues(),2*solver.eigenvectors().row(0).array().square().matrix().transpose()};
}
void psf(Eigen::MatrixXd& probability,Eigen::MatrixXd& derivative,int pixels,double origin,
         const Eigen::VectorXd& positions,double sigma,bool derivatives) {
    probability.resize(pixels,positions.size()); if(derivatives) derivative.resize(pixels,positions.size());
    const double root2=std::sqrt(2.), norm=1/(sigma*std::sqrt(2*std::numbers::pi));
    probability.setZero(); if(derivatives) derivative.setZero();
    for(int k=0;k<positions.size();++k) {
      // An eight-sigma footprint omits < 2e-15 of the Gaussian probability.
      const auto [first,last]=footprint(positions[k],origin,sigma,pixels);
      for(int i=first;i<=last;++i) {
        const double lo=(origin+i-.5-positions[k])/sigma,hi=lo+1/sigma;
        probability(i,k)=.5*(std::erf(hi/root2)-std::erf(lo/root2));
        if(derivatives) derivative(i,k)=norm*(std::exp(-.5*lo*lo)-std::exp(-.5*hi*hi));
      }
    }
}
}
BodyRenderer::BodyRenderer(BodyModel model,int order):model_(std::move(model)),order_(order) {
    if(model_.body.empty() || order<8 || order>64 || !model_.axes.allFinite() || (model_.axes.array()<=0).any() ||
       !model_.bodyToView.allFinite() || (model_.bodyToView*model_.bodyToView.transpose()-Eigen::Matrix3d::Identity()).norm()>1e-10)
        throw std::invalid_argument("Invalid registered body model or quadrature order.");
    const Eigen::Matrix3d q=model_.bodyToView*model_.axes.cwiseInverse().cwiseAbs2().asDiagonal()*model_.bodyToView.transpose();
    const Eigen::Matrix2d projected=q.topLeftCorner<2,2>()-q.topRightCorner<2,1>()*q.bottomLeftCorner<1,2>()/q(2,2);
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> solver(projected);
    const Eigen::Matrix2d ellipse=solver.eigenvectors()*solver.eigenvalues().array().sqrt().inverse().matrix().asDiagonal();
    const auto [r,w]=gaussLegendre(order); const int angular=4*order,n=order*angular;
    x_.resize(n);y_.resize(n);nx_.resize(n);ny_.resize(n);nz_.resize(n);weights_.resize(n);
    for(int i=0;i<order;++i)for(int j=0;j<angular;++j) {
        const int k=i*angular+j;const double radius=(r[i]+1)/2,angle=(j+.5)*2*std::numbers::pi/angular;
        const Eigen::Vector2d xy=ellipse*Eigen::Vector2d(radius*std::cos(angle),radius*std::sin(angle));
        const double z=(-q.bottomLeftCorner<1,2>().dot(xy)+std::sqrt(q(2,2)*(1-radius*radius)))/q(2,2);
        const Eigen::Vector3d normal=(q*Eigen::Vector3d(xy.x(),xy.y(),z)).normalized();
        x_[k]=xy.x();y_[k]=xy.y();nx_[k]=normal.x();ny_[k]=normal.y();nz_[k]=normal.z();
        weights_[k]=w[i]/2*radius*2*std::numbers::pi/angular*std::abs(ellipse.determinant());
    }
}
RenderedBody BodyRenderer::render(ImageSize size,PixelCoordinates origin,const BodyParameters& p,double sigma,
                                 Eigen::Vector2d smear,int exposureSamples,bool derivatives) const {
    if(size.width<1 || size.height<1 || std::int64_t(size.width)*size.height>16777216 ||
       !std::isfinite(origin.sample) || !std::isfinite(origin.line) || !p.vector().allFinite() ||
       p.radiusPixels<=0 || p.amplitudeDn<0 || p.phaseRadians<0 || p.phaseRadians>std::numbers::pi ||
       !std::isfinite(sigma) || sigma<.05 || sigma>10 || exposureSamples<1 || exposureSamples>101 || !smear.allFinite())
        throw std::invalid_argument("Invalid body rendering request.");
    RenderedBody result;result.dn=image::Image::Constant(size.height,size.width,p.backgroundDn);
    if(derivatives) result.jacobian=Eigen::Matrix<double,Eigen::Dynamic,7>::Zero(size.width*size.height,7);
    const Eigen::VectorXd horizontal=nx_*std::cos(p.sunAngleRadians)+ny_*std::sin(p.sunAngleRadians);
    const Eigen::VectorXd lit=std::sin(p.phaseRadians)*horizontal+std::cos(p.phaseRadians)*nz_;
    const Eigen::VectorXd active=(lit.array()>0).cast<double>();
    const Eigen::VectorXd base=p.radiusPixels*p.radiusPixels*weights_.array()*lit.cwiseMax(0).array();
    const Eigen::VectorXd flux=p.amplitudeDn*base;
    const Eigen::VectorXd phase=p.amplitudeDn*p.radiusPixels*p.radiusPixels*weights_.array()*active.array()*
        (std::cos(p.phaseRadians)*horizontal-std::sin(p.phaseRadians)*nz_).array();
    const Eigen::VectorXd azimuth=p.amplitudeDn*p.radiusPixels*p.radiusPixels*weights_.array()*active.array()*std::sin(p.phaseRadians)*
        (-nx_*std::sin(p.sunAngleRadians)+ny_*std::cos(p.sunAngleRadians)).array();
    for(int t=0;t<exposureSamples;++t) {
        const double fraction=(t+.5)/exposureSamples-.5;
        const Eigen::VectorXd sx=(p.center.sample+fraction*smear.x()+p.radiusPixels*x_.array()).matrix();
        const Eigen::VectorXd sy=(p.center.line+fraction*smear.y()+p.radiusPixels*y_.array()).matrix();
        Eigen::MatrixXd ex,ey,dx,dy;psf(ex,dx,size.width,origin.sample,sx,sigma,derivatives);psf(ey,dy,size.height,origin.line,sy,sigma,derivatives);
        // Accumulate compact PSF footprints directly. Dense matrix products do
        // mostly zero work here and make multi-image validation unnecessarily slow.
        const int count=size.width*size.height;
        for(int k=0;k<flux.size();++k) {
            const auto [x0,x1]=footprint(sx[k],origin.sample,sigma,size.width);
            const auto [y0,y1]=footprint(sy[k],origin.line,sigma,size.height);
            const double f=flux[k]/exposureSamples;
            for(int y=y0;y<=y1;++y)for(int x=x0;x<=x1;++x) {
                const int index=y*size.width+x;
                const double probability=ey(y,k)*ex(x,k),value=f*probability;
                result.dn.data()[index]+=value;
                if(derivatives) {
                    const double ds=f*ey(y,k)*dx(x,k),dl=f*dy(y,k)*ex(x,k);
                    double* jac=result.jacobian.data()+index;
                    jac[0*count]+=ds;jac[1*count]+=dl;
                    jac[2*count]+=2*value/p.radiusPixels+ds*x_[k]+dl*y_[k];
                    jac[3*count]+=phase[k]*probability/exposureSamples;
                    jac[4*count]+=azimuth[k]*probability/exposureSamples;
                    jac[5*count]+=base[k]*probability/exposureSamples;
                }
            }
        }
    }
    if(derivatives) result.jacobian.col(6).setOnes();
    return result;
}
} // namespace fd::opnav
