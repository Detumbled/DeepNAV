// Compare native GMAT observables with DeepNAV on GMAT-exported trajectories.
#include "observations/radiometric/TwoWayLink.hpp"
#include <SpiceUsr.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>
using namespace fd;
namespace radio = observations::radiometric;
struct Sample { double t; dynamics::CartesianState state; };
class Table final : public dynamics::StateProvider {
public:
    std::vector<Sample> rows;
    std::string_view referenceFrame() const noexcept override { return "J2000"; }
    std::string_view referenceOrigin() const noexcept override { return "SSB"; }
    dynamics::CartesianState stateAt(opnav::TdbEpoch epoch) const override {
        const double t = epoch.secondsPastJ2000;
        auto it = std::upper_bound(rows.begin(), rows.end(), t,
            [](double x, const Sample& y) { return x < y.t; });
        if (it == rows.begin() || it == rows.end())
            throw std::runtime_error("Epoch outside exported trajectory coverage");
        const auto& a = *std::prev(it); const auto& b = *it;
        const double h = b.t-a.t, u = (t-a.t)/h;
        dynamics::CartesianState result;
        result.positionKm = (2*u*u*u-3*u*u+1)*a.state.positionKm
          + (u*u*u-2*u*u+u)*h*a.state.velocityKmPerSec
          + (-2*u*u*u+3*u*u)*b.state.positionKm
          + (u*u*u-u*u)*h*b.state.velocityKmPerSec;
        result.velocityKmPerSec = ((6*u*u-6*u)/h)*a.state.positionKm
          + (3*u*u-4*u+1)*a.state.velocityKmPerSec
          + ((-6*u*u+6*u)/h)*b.state.positionKm
          + (3*u*u-2*u)*b.state.velocityKmPerSec;
        return result;
    }
};
double seconds(long double mjd) { return static_cast<double>((mjd-21545.L)*86400.L); }
int main(int argc, char** argv) {
 try {
    if (argc != 5) throw std::runtime_error("Usage: compare_gmat_two_way states.csv measurements.gmd delay_seconds naif0012.tls");
    erract_c("SET",0,const_cast<char*>("RETURN"));
    furnsh_c(argv[4]);
    if (failed_c()) throw std::runtime_error("Cannot load leap-second kernel");
    double arcStart;
    str2et_c("2004-10-10T18:13:06.046", &arcStart);
    if (failed_c()) throw std::runtime_error("Cannot convert the benchmark arc start");
    Table sc, gs;
    std::ifstream states(argv[1]); if (!states) throw std::runtime_error("Cannot open states");
    std::string line; std::getline(states,line);
    while (std::getline(states,line)) {
        std::replace(line.begin(),line.end(),',',' '); std::istringstream in(line);
        long double mjd; Sample a{},b{};
        if (!(in>>mjd)) throw std::runtime_error("Malformed state epoch");
        a.t=b.t=seconds(mjd);
        for(auto* row : {&a,&b}) {
            for(int j=0;j<3;++j) in>>row->state.positionKm[j];
            for(int j=0;j<3;++j) in>>row->state.velocityKmPerSec[j];
        }
        if (!in || !a.state.allFinite() || !b.state.allFinite()) throw std::runtime_error("Malformed state row");
        if (!sc.rows.empty() && a.t<=sc.rows.back().t) throw std::runtime_error("State epochs must increase");
        sc.rows.push_back(a); gs.rows.push_back(b);
    }
    if (sc.rows.size()<2) throw std::runtime_error("Missing trajectory samples");
    radio::TwoWayLink link(radio::libraryReceptionSolver(),{880./749.,std::stod(argv[3])});
    constexpr double f=7.175e9; const double reference=(880./749.)*f;
    std::ifstream data(argv[2]); if (!data) throw std::runtime_error("Cannot open measurements");
    std::cout<<std::setprecision(16)<<"tai_mjd,elapsed_tdb_s,type,gmat_converted,deepnav,error_m_or_hz\n";
    double maxRange=0,maxDoppler=0; int nr=0,nd=0;
    while(std::getline(data,line)) {
        if(line.empty() || line[0]=='%') continue;
        std::istringstream in(line); long double mjd; std::string type;
        if(!(in>>mjd>>type)) throw std::runtime_error("Malformed measurement");
        // GMAT's GMD tags are TAI modified Julian days, offset JD 2430000.0.
        double t=unitim_c(seconds(mjd),"TAI","TDB");
        if (failed_c()) throw std::runtime_error("SPICE time conversion failed");
        const auto last=line.find_last_of(" \t");
        double value=std::stod(line.substr(last+1)), predicted, error;
        if (!std::isfinite(value)) throw std::runtime_error("Non-finite measurement");
        if(type=="Range") {
            value*=.5; predicted=link.solve(t,gs,sc).rangeKm();
            error=(predicted-value)*1000.; maxRange=std::max(maxRange,std::abs(error)); ++nr;
        } else if(type=="DSN_TCP") {
            value=-value-reference;
            predicted=link.count(t,60.,gs,sc,{t,f,0},reference).residualHz;
            error=predicted-value; maxDoppler=std::max(maxDoppler,std::abs(error)); ++nd;
        } else throw std::runtime_error("Unexpected measurement type");
        std::cout<<mjd<<','<<t-arcStart<<','<<type<<','<<value<<','<<predicted<<','<<error<<'\n';
    }
    if(nr==0 || nr!=nd) throw std::runtime_error("Expected non-empty paired ranges and Doppler counts; check elevation/arc");
    std::cerr<<std::setprecision(9)<<"Max range error: "<<maxRange<<" m; max Doppler error: "<<maxDoppler<<" Hz\n";
    if(failed_c()) throw std::runtime_error("SPICE time conversion failed");
    return 0;
 } catch(const std::exception& e) { std::cerr<<e.what()<<'\n'; return 1; }
}
