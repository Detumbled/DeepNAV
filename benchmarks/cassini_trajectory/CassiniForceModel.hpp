#pragma once
#include "CassiniAttitude.hpp"
#include "Sat427Gravity.hpp"
#include "perturbations/AnnularGravity.hpp"
#include "perturbations/FlatPlateSRP.hpp"
#include "perturbations/Gravitational.hpp"
#include "perturbations/J2.hpp"
#include "perturbations/Relativistic.hpp"
#include "perturbations/ThermalRecoil.hpp"
#include "utils/CSPICE/SpiceError.hpp"
#include <vector>
#include <string>

namespace fd::simulation {
using CassiniState = Eigen::Matrix<double,6,1>;
inline CassiniState cassiniRelativeState(const std::string& body,double t) {
    double x[6],lt;
    spkezr_c(body.c_str(),t,"J2000","NONE","SATURN",x,&lt);
    od::throwIfSpiceFailed("Saturn-relative state of "+body);
    return Eigen::Map<CassiniState>(x);
}
inline double cassiniKernelGM(const std::string& body) {
    double mu;SpiceInt n;bodvrd_c(body.c_str(),"GM",1,&n,&mu);
    od::throwIfSpiceFailed("GM of "+body);return mu;
}
inline std::vector<std::pair<std::string,double>> cassiniPerturbers(const Sat427Gravity& gravity) {
    std::vector<std::pair<std::string,double>> bodies;
    for (const char* body:{"SUN","1","2","3","4","5","7","8","9","MIMAS","ENCELADUS","TETHYS","DIONE","RHEA","TITAN","HYPERION","IAPETUS","PHOEBE","HELENE"}) {
        double mu=gravity.resolvedGM(body,std::string(body)=="HELENE" ? 0.0 : cassiniKernelGM(body));
        if (std::string(body)=="SUN") mu=gravity.value("GM10")-cassiniKernelGM("1")-cassiniKernelGM("2")-cassiniKernelGM("3")-cassiniKernelGM("4");
        bodies.emplace_back(body,mu);
    }
    return bodies;
}
// Shared forces for the retained open-loop comparison and the two-way WLS.
// km, km/s, TDB. Optional SRP/recoil scales default to the retained nominal
// forces. evaluate() never resets a spacecraft state from SPK.
class CassiniForceModel {
    const Sat427Gravity& gravity_;
    const CassiniAttitude& attitude_;
    std::vector<std::pair<std::string,double>> bodies_;
    double mu_,ringGM_,originGM_,mass_,area_,nu_,shadowRadius_;
    int mask_;
    perturbations::ThermalRecoil thermal_;
    std::vector<perturbations::AnnularGravity> rings_;
public:
    struct Evaluation : perturbations::AccelerationEvaluation {
        Eigen::Vector3d srpScaleJacobian=Eigen::Vector3d::Zero();
        Eigen::Vector3d thermalScaleJacobian=Eigen::Vector3d::Zero();
    };
    struct GravityTerms {
        Eigen::Vector3d newtonian=Eigen::Vector3d::Zero(),indirectJ2=Eigen::Vector3d::Zero(),ring=Eigen::Vector3d::Zero(),relativistic=Eigen::Vector3d::Zero(),sun;
        Eigen::Matrix3d newtonianJacobian=Eigen::Matrix3d::Zero();
    };
    CassiniForceModel(const Sat427Gravity& gravity,const CassiniAttitude& attitude,
                     std::vector<std::pair<std::string,double>> bodies,int mask,double mass,
                     double area,double cr,const Eigen::Vector3d& thermalReference,double thermalEpoch,double shadowRadius)
        :gravity_(gravity),attitude_(attitude),bodies_(std::move(bodies)),mu_(gravity.massGM("SATURN")),
        ringGM_(gravity.value("RINGM1")+gravity.value("RINGM2")+gravity.value("RINGM4")),originGM_(mu_+((mask&2)?ringGM_:0)),
        mass_(mass),area_(area),nu_((cr-1)/2),shadowRadius_(shadowRadius),mask_(mask),thermal_(thermalReference,thermalEpoch,2781),
        rings_{{gravity.value("RINGM1"),gravity.value("C_ring_RI"),gravity.value("C_ring_RO")},
               {gravity.value("RINGM2"),gravity.value("B_ring_RI"),gravity.value("B_ring_RO")},
               {gravity.value("RINGM4"),gravity.value("A_ring_RI"),gravity.value("A_ring_RO")}} {}
    GravityTerms gravityTerms(double et,const Eigen::Vector3d& r,const Eigen::Vector3d& v,bool auditAll=false) const {
        GravityTerms terms;const auto pole=gravity_.pole(et);
        Eigen::Vector3d originVelocity=Eigen::Vector3d::Zero();
        std::vector<perturbations::RelativisticBody> pn;
        if (auditAll || (mask_&4)) {
            double x[6],lt;spkezr_c("SATURN",et,"J2000","NONE","SOLAR SYSTEM BARYCENTER",x,&lt);
            od::throwIfSpiceFailed("Saturn barycentric velocity");originVelocity=Eigen::Map<CassiniState>(x).tail<3>();
            pn.push_back({mu_,Eigen::Vector3d::Zero(),originVelocity});
        }
        for (std::size_t i=0;i<bodies_.size();++i) {
            const auto& body=bodies_[i];const auto x=cassiniRelativeState(body.first,et);const Eigen::Vector3d rb=x.head<3>();
            if (i==0) terms.sun=rb;
            const auto tide=perturbations::thirdBodyGravity(body.second,r,rb);
            terms.newtonian+=tide.acceleration;terms.newtonianJacobian+=tide.positionJacobian;
            if (auditAll || (mask_&4)) pn.push_back({body.second,rb,x.tail<3>()+originVelocity});
            if (i>=9) {
                if (auditAll || (mask_&1)) terms.indirectJ2+=(mu_/originGM_)*perturbations::j2OriginCorrection(mu_,body.second,gravity_.value("RADIUS"),gravity_.value("J602"),rb,pole).acceleration;
                if (auditAll || (mask_&2)) for (const auto& ring:rings_) terms.ring+=body.second/(mu_+ringGM_)*ring.acceleration(rb,pole,true);
            }
        }
        if (auditAll || (mask_&2)) for (const auto& ring:rings_) terms.ring+=ring.acceleration(r,pole);
        if (auditAll || (mask_&4)) {
            const perturbations::RelativisticSystem system(std::move(pn));
            terms.relativistic=system.correction(r,v+originVelocity)-system.correction(Eigen::Vector3d::Zero(),originVelocity,0);
        }
        return terms;
    }
    Evaluation evaluate(double et,const Eigen::Vector3d& r,const Eigen::Vector3d& v,bool withPartials=false,double srpScale=1,double thermalScale=1) const {
        if (withPartials && (mask_&4)) throw std::invalid_argument("Variational Cassini model requires relativity disabled");
        Evaluation total;static_cast<perturbations::AccelerationEvaluation&>(total)=perturbations::pointMassGravity(mu_)(et,r,v);const auto terms=gravityTerms(et,r,v);
        total.acceleration+=terms.newtonian;total.positionJacobian+=terms.newtonianJacobian;
        if (mask_&1) total.acceleration+=terms.indirectJ2;
        if (mask_&2) total.acceleration+=terms.ring;
        if (mask_&4) total.acceleration+=terms.relativistic;
        const auto pole=gravity_.pole(et);const auto j2=perturbations::j2Gravity(mu_,gravity_.value("RADIUS"),gravity_.value("J602"),r,pole);
        total.acceleration+=j2.acceleration;total.positionJacobian+=j2.positionJacobian;
        const Eigen::Matrix3d rotation=attitude_.rotation(et);
        total.thermalScaleJacobian=thermal_.evaluate(et,mass_,rotation).acceleration;
        total.acceleration+=thermalScale*total.thermalScaleJacobian;
        const auto srp=perturbations::flatPlateSRPWithShadow(r,terms.sun,-rotation.col(2),area_,mass_,nu_,shadowRadius_,695700);
        total.srpScaleJacobian=srp.acceleration;
        total.acceleration+=srpScale*srp.acceleration;total.positionJacobian+=srpScale*srp.positionJacobian;
        if (withPartials && (mask_&2)) for (int j=0;j<3;++j) {
            // Exterior ring field is smooth; 1 km is tiny versus arc radii.
            Eigen::Vector3d plus=r,minus=r;plus[j]+=1;minus[j]-=1;
            for (const auto& ring:rings_) total.positionJacobian.col(j)+=(ring.acceleration(plus,pole)-ring.acceleration(minus,pole))/2;
        }
        return total;
    }
};
} // namespace fd::simulation
