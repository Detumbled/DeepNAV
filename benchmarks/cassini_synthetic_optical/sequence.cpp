#include "../cassini_trajectory/SavedSimulation.hpp"
#include "opnav/synth/SyntheticImageGenerator.hpp"
#include "opnav/core/ApparentDirection.hpp"
#include "dynamics/SpiceKernelSet.hpp"
#include "dynamics/SpiceStateProvider.hpp"
#include "simulation/CassiniRadioKernels.hpp"
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numbers>

namespace fs=std::filesystem;
using namespace fd::opnav;
std::ofstream output(const fs::path& path) {
    std::ofstream out(path); out.exceptions(std::ios::badbit|std::ios::failbit);
    out << std::setprecision(17); return out;
}
int main(int argc,char** argv) {
    try {
        const fs::path root=DEEPNAV_SOURCE_DIR;
        fs::path directory=root/"results/cassini_synthetic_optical/simulated_sequence";
        for(int i=1;i<argc;++i) {
            const std::string key=argv[i];
            if(key=="--help") {std::cout << "cassini_synthetic_sequence [--output DIR]\n";return 0;}
            if(key!="--output" || ++i==argc) throw std::invalid_argument("Expected --output DIR.");
            directory=fs::absolute(argv[i]);
        }
        od::SpiceErrorModeGuard errors("RETURN","NONE");
        auto paths=fd::simulation::cassiniRadioKernels(root/"Kernels");
        for(const char* file:{"Cassini/fk/cas_v40.tf","Cassini/sclk/cas00172.tsc",
            "Cassini/ck/04256_04292py_as_flown.bc","Cassini/ck/04292_04320py_as_flown.bc",
            "Cassini/ck/04281_04286ra.bc","Cassini/ck/04286_04291ra.bc","Cassini/ck/04291_04296ra.bc",
            "Cassini/ck/04296_04301ra.bc","Cassini/ck/04293_04298ra.bc"}) paths.push_back(root/"Kernels"/file);
        const fd::dynamics::SpiceKernelSet kernels(paths);
        // Same retained-truth provider and force parameters as cassini_one_way_ekf.
        // CK attitude is needed only for force-model interpolation derivatives;
        // optical camera pointing below is independently declared and synthetic.
        const CassiniAttitude forceAttitude(paths);
        const Sat427Gravity gravity(root/"Kernels/Cassini/spk/200128RU_SCPSE_04251_04303.bsp.lbl");
        double recoilEpoch; str2et_c("2005-08-13T00:00:00",&recoilEpoch);
        double radii[3]; SpiceInt count; bodvrd_c("SATURN","RADII",3,&count,radii);
        od::throwIfSpiceFailed("Sequence force constants");
        const fd::simulation::CassiniForceModel force(gravity,forceAttitude,fd::simulation::cassiniPerturbers(gravity),
            2,3220,12,1.42,{-2.725e-13,-2.725e-13,-5.232237415805135e-12},recoilEpoch,radii[0]);
        const fd::dynamics::SpiceStateProvider saturn("SATURN"),target("ENCELADUS"),spkCassini("CASSINI");
        const fs::path truthPath=root/"results/cassini_two_way_wls/pass_schedule/truth/trajectory.csv";
        const fd::simulation::SavedSimulation truth(truthPath,force,saturn);
        fs::create_directories(directory);
        auto observations=output(directory/"observations.csv");
        observations << "tdb_seconds,target,observer,sample,line,cov_ss,cov_sl,cov_ll,seed";
        for(int i=0;i<3;++i)for(int j=0;j<3;++j) observations << ",rotation_" << i << j;
        observations << '\n';
        auto evaluated=output(directory/"truth.csv");
        evaluated << "tdb_seconds,ssb_x_km,ssb_y_km,ssb_z_km,ssb_vx_km_s,ssb_vy_km_s,ssb_vz_km_s,saturn_x_km,saturn_y_km,saturn_z_km,saturn_vx_km_s,saturn_vy_km_s,saturn_vz_km_s,spk_position_difference_km\n";
        auto summary=output(directory/"summary.json");
        summary << "{\n\"observer_source\": \"SavedSimulation: multibody_j2_srp_thermal; same provider as Cassini EKF\","
                << "\n\"truth_csv\": " << std::quoted(truthPath.string())
                << ",\n\"cadence_utc_days\": 1, \"exposure_seconds\": 0.06,"
                << "\n\"camera\": {\"focal_length_mm\": 2003.44, \"pixel_pitch_mm\": 0.012, \"principal_point\": [63.23,62.74]},"
                << "\n\"pointing\": \"predeclared at each exposure from simulated truth; fixed during exposure; not a function of EKF estimate\","
                << "\n\"frames\": [\n";
        for(int index=0;index<4;++index) {
            const std::string utcEpoch="2004-10-"+std::to_string(10+index)+"T00:00:00";
            double t; str2et_c(utcEpoch.c_str(),&t); od::throwIfSpiceFailed("Sequence epoch");
            const TdbEpoch epoch(t);
            if(t+0.03>truth.endEpoch()) throw std::runtime_error("Sequence outside retained truth coverage.");
            const auto state=truth.stateAt(epoch);
            const auto direction=core::computeApparentDirection(epoch,state,target);
            const Eigen::Vector3d z=direction.unitDirection,x=z.unitOrthogonal(),y=z.cross(x);
            Eigen::Matrix3d rotation;rotation.row(0)=x.transpose();rotation.row(1)=y.transpose();rotation.row(2)=z.transpose();
            const CameraAttitude attitude(rotation);
            const CameraModel camera({2003.44,{63.23,62.74},-Eigen::Matrix2d::Identity()/0.012});
            synth::SyntheticImageOptions options; options.psfSigmaPixels=0.8;
            const synth::SyntheticImageGenerator generator(camera,[attitude](TdbEpoch){return attitude;},options);
            synth::OpticalDetectorNoise noise;noise.electronsPerDn=10;noise.backgroundElectrons=50;noise.biasDn=100;noise.seed=4200+index;
            const auto frame=generator.createImage("ENCELADUS",epoch,noise,truth,"SIMULATED_CASSINI");
            if((frame.truth.observerPositionKm-state.positionKm).norm()>1e-10 ||
               (frame.truth.observerVelocityKmPerSecond-state.velocityKmPerSec).norm()>1e-13)
                throw std::runtime_error("Image observer differs from simulated truth.");
            const auto fit=image::fitBrightestCircularGaussian(frame.dn,frame.varianceDn2,12,{},&frame.mask);
            if(!fit.measurement || !frame.truth.brightnessCenter) throw std::runtime_error("Missing source or fit in sequence.");
            const auto pixel=generator.createMeasurement("ENCELADUS",epoch,{0.15,std::uint64_t(5000+index)},truth,"SIMULATED_CASSINI");
            observations << t << ",ENCELADUS,SIMULATED_CASSINI," << pixel.measurement.center.sample << ',' << pixel.measurement.center.line
                         << ",0.0225,0,0.0225," << pixel.truth.seed;
            for(int i=0;i<3;++i)for(int j=0;j<3;++j) observations << ',' << rotation(i,j);
            observations << '\n';
            const auto relative=truth.relative(t);
            const double difference=(spkCassini.stateAt(epoch).positionKm-state.positionKm).norm();
            evaluated << t;
            for(int j=0;j<6;++j) evaluated << ',' << (j<3 ? state.positionKm[j] : state.velocityKmPerSec[j-3]);
            for(int j=0;j<6;++j) evaluated << ',' << relative[j];
            evaluated << ',' << difference << '\n';
            const std::string filename="frame_"+std::to_string(index)+".csv";
            auto pixels=output(directory/filename);
            pixels << "sample,line,dn,variance_dn2,mask,source_electrons\n";
            for(int l=0;l<frame.dn.rows();++l)for(int s=0;s<frame.dn.cols();++s)
                pixels << s << ',' << l << ',' << frame.dn(l,s) << ',' << frame.varianceDn2(l,s) << ',' << int(frame.mask(l,s)) << ',' << frame.expectedSourceElectrons(l,s) << '\n';
            SpiceChar utc[64];et2utc_c(t,"ISOC",6,sizeof(utc),utc);od::throwIfSpiceFailed("Sequence UTC");
            summary << (index ? ",\n" : "") << "{\"file\": " << std::quoted(filename) << ", \"utc\": " << std::quoted(utc)
                    << ", \"tdb_seconds\": " << t << ", \"seed\": " << noise.seed
                    << ", \"phase_degrees\": " << frame.truth.phaseAngleRadians*180/std::numbers::pi
                    << ", \"geometric_center\": [" << frame.truth.geometricCenter.sample << ',' << frame.truth.geometricCenter.line << ']'
                    << ", \"brightness_center\": [" << frame.truth.brightnessCenter->sample << ',' << frame.truth.brightnessCenter->line << ']'
                    << ", \"fitted_center\": [" << fit.model.center.sample << ',' << fit.model.center.line << ']'
                    << ", \"fit_covariance\": [[" << fit.measurement->covariance(0,0) << ',' << fit.measurement->covariance(0,1)
                    << "],[" << fit.measurement->covariance(1,0) << ',' << fit.measurement->covariance(1,1) << "]]"
                    << ", \"fit_reduced_chi_squared\": " << fit.chiSquared/fit.degreesOfFreedom
                    << ", \"spk_position_difference_km\": " << difference << '}';
            std::cout << utc << ": simulated/SPK spacecraft separation " << difference << " km\n";
        }
        summary << "\n],\n\"loaded_kernels\": [";
        for(std::size_t i=0;i<paths.size();++i) summary << (i ? "," : "") << std::quoted(paths[i].string());
        summary << "]\n}\n";
        std::cout << "Saved four simulated-truth frames and direct pixel observations in " << directory << '\n';
    } catch(const std::exception& e) {std::cerr << "ERROR: " << e.what() << '\n';return 1;}
}
