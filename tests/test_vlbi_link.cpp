#include "dynamics/LinearStateProvider.hpp"
#include "observations/radiometric/VlbiMeasurement.hpp"
#include <iostream>
#include <cmath>
#include <stdexcept>
namespace r=fd::observations::radiometric;
using Linear=fd::dynamics::LinearStateProvider;
void near(double x,double y,double e,const char* what) {if (!std::isfinite(x)||std::abs(x-y)>e) {std::cerr<<what<<": "<<x<<" vs "<<y<<'\n';throw std::runtime_error(what);}}
Linear state(double x,double v=0,double t=0) {return {fd::opnav::TdbEpoch{t},{{x,0,0},{v,0,0}},"J2000","SSB"};}
int main() {try {
    const double c=r::speedOfLightKmPerSecond;const r::VlbiLink link(r::libraryReceptionSolver());
    for (double epoch:{0.,1.5e8,-1.5e8}) {
        const auto probe=state(1e9,0,epoch),one=state(1e8,0,epoch),two=state(1e8-6000,0,epoch);
        const auto a=link.solve(epoch,probe,one,two);near(a.differentialRangeKm,6000,1e-8,"Static interplanetary baseline");
        near(link.solve(epoch,probe,two,one).differentialRangeKm,-6000,1e-8,"Baseline sign");
        near(a.first.emissionTdb,a.second.emissionTdb,0,"Common emission event");
        const auto movingProbe=state(1e9,50,epoch);
        near(link.solve(epoch,movingProbe,one,two).differentialRangeKm,6000,1e-8,"Moving spacecraft must share one emission event");
        const auto moving=state(1e8-6000,30,epoch);
        near(link.solve(epoch,probe,one,moving).differentialRangeKm,6000*c/(c+30),5e-8,"Moving receiving station at large epoch");
        const r::DelayModel first=[](const auto&,const auto&,const auto&){return 1e-5;},second=[](const auto&,const auto&,const auto&){return 2e-5;};
        const r::VlbiLink delayed(r::libraryReceptionSolver({},first),second);
        near(delayed.solve(epoch,probe,one,two).differentialRangeKm,6000+c*1e-5,1e-8,"Differential propagation correction");
    }
    const Linear one{fd::opnav::TdbEpoch{0},{{0,0,0},{0,0,0}},"J2000","SSB"},two{fd::opnav::TdbEpoch{0},{{6000,3000,0},{0,0,0}},"J2000","SSB"};
    Eigen::VectorXd x(8);x<<1e9,2e8,1e8,0,0,0,1,1;
    const r::TrajectoryFactory factory=[](double epoch,const Eigen::VectorXd& p) {return std::make_unique<Linear>(fd::opnav::TdbEpoch{epoch},fd::dynamics::CartesianState{p.head<3>(),p.segment<3>(3)},"J2000","SSB");};
    const r::VlbiMeasurement m(link,one,two,0,factory,Eigen::VectorXd::Constant(8,1));
    const auto prediction=m(0,x);const Eigen::RowVector3d expected=((x.head<3>()-two.stateAt(fd::opnav::TdbEpoch{0}).positionKm).normalized()-x.head<3>().normalized()).transpose();
    if ((prediction.jacobian.leftCols(3)-expected).norm()>2e-9) throw std::runtime_error("VLBI transverse sensitivity");
    near(prediction.jacobian(0,6),0,0,"Unused coefficient sensitivity");
    bool rejected=false;try {(void)link.solve(0,*factory(0,x),one,one);} catch(const std::invalid_argument&) {rejected=true;}if (!rejected) throw std::runtime_error("Identical stations accepted");
    const Linear wrong{fd::opnav::TdbEpoch{0},{{1,0,0},{0,0,0}},"J2000","SUN"};rejected=false;try {(void)link.solve(0,*factory(0,x),one,wrong);} catch(const std::invalid_argument&) {rejected=true;}if (!rejected) throw std::runtime_error("Mixed origins accepted");
    std::cout<<"VLBI common-wavefront, corrections, precision and Jacobian checks passed.\n";return 0;
} catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}}
