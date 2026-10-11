#pragma once
#include <Eigen/Core>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace fd::perturbations {
struct RelativisticBody {
    double gm;
    Eigen::Vector3d position, velocity;
};

// GR (beta=gamma=1), Moyer Eq. 4-26. Barycentric km, km/s, TDB seconds.
// Returns only the first post-Newtonian acceleration, not Newtonian gravity.
// Prescribed point masses; source accelerations may be evaluated Newtonianly.
// No variational partials are supplied by this acceleration-only interface.
class RelativisticSystem {
    std::vector<RelativisticBody> bodies_;
    std::vector<double> potentials_;
    std::vector<Eigen::Vector3d> accelerations_;
public:
    explicit RelativisticSystem(std::vector<RelativisticBody> bodies)
        : bodies_(std::move(bodies)), potentials_(bodies_.size(), 0),
          accelerations_(bodies_.size(), Eigen::Vector3d::Zero()) {
        for (const auto& b : bodies_)
            if (!(b.gm > 0) || !std::isfinite(b.gm) ||
                !b.position.allFinite() || !b.velocity.allFinite())
                throw std::invalid_argument("Invalid relativistic source");
        for (std::size_t j=0;j<bodies_.size();++j)
            for (std::size_t k=0;k<j;++k) {
                const Eigen::Vector3d d=bodies_[k].position-bodies_[j].position;
                const double r=d.norm();
                if (!(r>0)) throw std::invalid_argument("Coincident relativistic sources");
                potentials_[j]+=bodies_[k].gm/r;
                potentials_[k]+=bodies_[j].gm/r;
                accelerations_[j]+=bodies_[k].gm/(r*r*r)*d;
                accelerations_[k]-=bodies_[j].gm/(r*r*r)*d;
            }
    }
    // For a massive source target, omit its index from the outer sum only;
    // its mass remains in the other sources' potentials and accelerations.
    Eigen::Vector3d correction(const Eigen::Vector3d& position,
                              const Eigen::Vector3d& velocity,
                              int omit=-1) const {
        if (!position.allFinite() || !velocity.allFinite() || omit < -1 ||
            omit >= static_cast<int>(bodies_.size()))
            throw std::invalid_argument("Invalid relativistic target");
        constexpr double c=299792.458;
        double potential=0;
        for (std::size_t j=0;j<bodies_.size();++j) if (static_cast<int>(j)!=omit) {
            const double distance=(position-bodies_[j].position).norm();
            if (!(distance>0)) throw std::invalid_argument("Target on relativistic source");
            potential+=bodies_[j].gm/distance;
        }
        Eigen::Vector3d result=Eigen::Vector3d::Zero();
        for (std::size_t j=0;j<bodies_.size();++j) if (static_cast<int>(j)!=omit) {
            const auto& b=bodies_[j];
            const Eigen::Vector3d d=b.position-position;
            const double r=d.norm(), projection=d.dot(b.velocity)/r;
            const double bracket=-4*potential-potentials_[j]+velocity.squaredNorm()
                +2*b.velocity.squaredNorm()-4*velocity.dot(b.velocity)
                -1.5*projection*projection+0.5*d.dot(accelerations_[j]);
            result+=b.gm/(c*c*r*r*r)*(d*bracket
                +(-d).dot(4*velocity-3*b.velocity)*(velocity-b.velocity));
            result+=3.5*b.gm/(c*c*r)*accelerations_[j];
        }
        return result;
    }
};
} // namespace fd::perturbations
