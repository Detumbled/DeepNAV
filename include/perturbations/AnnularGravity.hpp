#pragma once
#include <Eigen/Core>
#include <array>
#include <cmath>
#include <stdexcept>

namespace fd::perturbations {
// Uniform, infinitesimally thin annulus. Exterior Legendre expansion of the
// Newton integral; degree 80. Valid only r >= 1.2 * outer radius.
// GM in km^3/s^2, radii/positions in km, acceleration in km/s^2.
// This models separate ring mass, not extra Saturn gravity coefficients.
class AnnularGravity {
    double gm_, outer_;
    std::array<double,81> coefficients_{};
public:
    AnnularGravity(double gm, double inner, double outer):gm_(gm),outer_(outer) {
        if (!std::isfinite(gm) || gm<0 || !std::isfinite(inner) ||
            !std::isfinite(outer) || inner<0 || !(outer>inner))
            throw std::invalid_argument("Invalid annulus");
        const double q=inner/outer;
        double p0=1;
        coefficients_[0]=1;
        for (int l=2;l<=80;l+=2) {
            p0*=-(l-1.0)/l; // P_l(0)
            // Area-weighted radial moment / outer^l, stable for narrow annuli.
            const double moment=q==0 ? 1/(l/2.0+1) :
                std::expm1((l+2)*std::log(q))/std::expm1(2*std::log(q))/(l/2.0+1);
            coefficients_[l]=p0*moment;
        }
    }
    // Optionally return only the non-monopole field for origin reactions.
    Eigen::Vector3d acceleration(const Eigen::Vector3d& position,
                                const Eigen::Vector3d& pole, bool shapeOnly=false) const {
        const double r=position.norm();
        if (!position.allFinite() || !pole.allFinite() ||
            std::abs(pole.squaredNorm()-1)>1e-12 || !std::isfinite(r) || r<1.2*outer_)
            throw std::invalid_argument("Annulus requires exterior point and unit pole");
        const Eigen::Vector3d n=position/r;
        const double u=n.dot(pole), q=outer_/r;
        Eigen::Vector3d sum=shapeOnly ? Eigen::Vector3d::Zero().eval() : (-n).eval();
        double pm2=1, pm1=u, dm2=0, dm1=1, power=q;
        for (int l=2;l<=80;++l) {
            const double p=((2*l-1)*u*pm1-(l-1)*pm2)/l;
            const double dp=((2*l-1)*(pm1+u*dm1)-(l-1)*dm2)/l;
            power*=q;
            if (l%2==0) sum+=coefficients_[l]*power*
                (-(l+1)*p*n+dp*(pole-u*n));
            pm2=pm1;pm1=p;dm2=dm1;dm1=dp;
        }
        return gm_/(r*r)*sum;
    }
};
} // namespace fd::perturbations
