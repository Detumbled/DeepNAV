#include "opnav/synth/SyntheticImageGenerator.hpp"
#include "dynamics/SpiceKernelSet.hpp"
#include "dynamics/SpiceStateProvider.hpp"
#include "opnav/core/ApparentDirection.hpp"
#include "utils/CSPICE/SpiceError.hpp"
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numbers>

namespace fs = std::filesystem;
using namespace fd::opnav;
namespace {
std::ofstream output(const fs::path& path) {
    std::ofstream out(path); out.exceptions(std::ios::badbit|std::ios::failbit);
    out << std::setprecision(17); return out;
}
}
int main(int argc, char** argv) {
    try {
        const fs::path root = DEEPNAV_SOURCE_DIR;
        fs::path directory = root/"results/cassini_synthetic_optical";
        for (int i=1; i<argc; ++i) {
            const std::string arg=argv[i];
            if (arg=="--help") { std::cout << "cassini_synthetic_optical [--output DIR]\n"; return 0; }
            if (arg!="--output" || ++i==argc) throw std::invalid_argument("Expected --output DIR.");
            directory=fs::absolute(argv[i]);
        }
        od::SpiceErrorModeGuard guard("RETURN","NONE");
        fs::current_path(root);
        const fs::path meta=root/"benchmarks/cassini_synthetic_optical/kernels.tm";
        const fd::dynamics::SpiceKernelSet kernels({meta});
        double et; str2et_c("2004-10-10T18:13:06.046",&et); od::throwIfSpiceFailed("Synthetic epoch");
        const TdbEpoch epoch(et);
        const fd::dynamics::SpiceStateProvider cassini("CASSINI"), enceladus("ENCELADUS");
        const auto direction=core::computeApparentDirection(epoch,cassini.stateAt(epoch),enceladus);
        // One fixed synthetic attitude, set once from the reference epoch.
        // It is not recomputed from each observation or an estimated spacecraft state.
        const Eigen::Vector3d z=direction.unitDirection, x=z.unitOrthogonal(), y=z.cross(x);
        Eigen::Matrix3d rotation; rotation.row(0)=x.transpose(); rotation.row(1)=y.transpose(); rotation.row(2)=z.transpose();
        const CameraAttitude attitude(rotation);
        const CameraModel camera({2003.44,{63.23,62.74},-Eigen::Matrix2d::Identity()/0.012});
        synth::SyntheticImageOptions options;
        options.psfSigmaPixels=0.8;
        synth::SyntheticImageGenerator generator(camera,[attitude](TdbEpoch){ return attitude; },options);
        synth::OpticalDetectorNoise noise;
        noise.electronsPerDn=10; noise.backgroundElectrons=50; noise.biasDn=100; noise.seed=42;
        const auto frame=generator.create("ENCELADUS",epoch,noise,"CASSINI");
        const auto fit=image::fitBrightestCircularGaussian(frame.dn,frame.varianceDn2,12,{},&frame.mask);
        if (!fit.measurement) throw std::runtime_error("Synthetic image centroid fit failed.");
        const auto& truth=frame.truth;
        if (!truth.brightnessCenter) throw std::runtime_error("No source flux in synthetic detector.");
        // Independent CSPICE reference checks center geometry; same synthetic attitude.
        double vector[3],lt; spkpos_c("ENCELADUS",et,"J2000","CN+S","CASSINI",vector,&lt);
        od::throwIfSpiceFailed("Independent synthetic reference");
        const auto reference=camera.project(Eigen::Map<Eigen::Vector3d>(vector),attitude);
        const double referenceError=std::hypot(reference.sample-truth.geometricCenter.sample,reference.line-truth.geometricCenter.line);
        if (referenceError>0.01) throw std::runtime_error("Synthetic center disagrees with CSPICE.");
        fs::create_directories(directory);
        auto pixels=output(directory/"frame.csv");
        pixels << "sample,line,dn,variance_dn2,masked,source_electrons,expected_dn\n";
        for (int l=0;l<frame.dn.rows();++l) for(int s=0;s<frame.dn.cols();++s)
            pixels << s << ',' << l << ',' << frame.dn(l,s) << ',' << frame.varianceDn2(l,s) << ',' << int(frame.mask(l,s)) << ','
                   << frame.expectedSourceElectrons(l,s) << ',' << noise.biasDn+(frame.expectedSourceElectrons(l,s)+noise.backgroundElectrons)/noise.electronsPerDn << '\n';
        auto measurements=output(directory/"measurements.csv");
        measurements << "offset_seconds,tdb_seconds,truth_sample,truth_line,sample,line,variance_sample,variance_line,seed\n";
        for(int i=0;i<61;++i) {
            const auto observation=generator.createMeasurement("ENCELADUS",TdbEpoch(et+i-30),{0.15,std::uint64_t(1000+i)},"CASSINI");
            measurements << i-30 << ',' << et+i-30 << ',' << observation.truth.geometricCenter.sample << ',' << observation.truth.geometricCenter.line << ','
                         << observation.measurement.center.sample << ',' << observation.measurement.center.line << ','
                         << observation.measurement.covariance(0,0) << ',' << observation.measurement.covariance(1,1) << ',' << observation.truth.seed << '\n';
        }
        auto summary=output(directory/"summary.json");
        summary << "{\n\"body\": \"ENCELADUS\", \"observer\": \"CASSINI\",\n\"mid_time_utc\": \"2004-10-10T18:13:06.046\","
                << "\n\"tdb_seconds\": " << et << ",\n\"seed\": " << noise.seed
                << ",\n\"attitude_source\": \"fixed synthetic attitude; no CK pointing used\",\n\"meta_kernel\": " << std::quoted(meta.string())
                << ",\n\"geometric_center\": [" << truth.geometricCenter.sample << ',' << truth.geometricCenter.line << ']'
                << ",\n\"brightness_center\": [" << truth.brightnessCenter->sample << ',' << truth.brightnessCenter->line << ']'
                << ",\n\"fitted_center\": [" << fit.model.center.sample << ',' << fit.model.center.line << ']'
                << ",\n\"fit_reduced_chi_squared\": " << fit.chiSquared/fit.degreesOfFreedom
                << ",\n\"range_km\": " << truth.rangeKm << ",\n\"radius_km\": " << truth.radiusKm
                << ",\n\"light_time_seconds\": " << truth.lightTimeSeconds
                << ",\n\"phase_angle_degrees\": " << truth.phaseAngleRadians*180/std::numbers::pi
                << ",\n\"cspice_center_difference_pixels\": " << referenceError
                << ",\n\"source_electrons\": " << frame.expectedSourceElectrons.sum()
                << ",\n\"full_phase_electrons\": " << options.fullPhaseElectrons
                << ",\n\"psf_sigma_pixels\": " << options.psfSigmaPixels
                << ",\n\"exposure_seconds\": " << options.exposureSeconds
                << ",\n\"exposure_samples\": " << options.exposureSamples << ",\n\"subpixels\": " << options.subpixels
                << ",\n\"gain_electrons_per_dn\": " << noise.electronsPerDn
                << ",\n\"focal_length_mm\": 2003.44, \"pixel_pitch_mm\": 0.012"
                << ",\n\"principal_point\": [63.23,62.74]"
                << ",\n\"background_electrons\": " << noise.backgroundElectrons
                << ",\n\"read_noise_electrons\": " << noise.readNoiseElectrons
                << ",\n\"bias_dn\": " << noise.biasDn << ",\n\"saturation_dn\": " << noise.saturationDn
                << ",\n\"shot_noise\": true, \"quantize\": true,\n\"observer_state_j2000_ssb\": [";
        for(int i=0;i<6;++i) summary << (i ? "," : "") << (i<3 ? truth.observerPositionKm[i] : truth.observerVelocityKmPerSecond[i-3]);
        summary << "],\n\"inertial_to_camera\": [";
        for(int i=0;i<3;++i) { summary << (i ? "," : "") << '['; for(int j=0;j<3;++j) summary << (j ? "," : "") << rotation(i,j); summary << ']'; }
        summary << "],\n\"loaded_kernels\": [";
        SpiceInt count; ktotal_c("ALL",&count);
        for (SpiceInt i=0;i<count;++i) {
            SpiceChar file[2048],type[32],source[2048]; SpiceInt handle; SpiceBoolean found;
            kdata_c(i,"ALL",sizeof(file),sizeof(type),sizeof(source),file,type,source,&handle,&found);
            od::throwIfSpiceFailed("Synthetic kernel provenance");
            if (!found) throw std::runtime_error("Missing kernel inventory entry.");
            summary << (i ? "," : "") << std::quoted(fs::absolute(file).string());
        }
        summary << "]\n}\n";
        std::cout << "Synthetic geometry vs CSPICE: " << referenceError << " px\n"
                  << "Geometric center: " << truth.geometricCenter.sample << ", " << truth.geometricCenter.line << '\n'
                  << "Brightness center: " << truth.brightnessCenter->sample << ", " << truth.brightnessCenter->line << '\n'
                  << "Circular fit: " << fit.model.center.sample << ", " << fit.model.center.line << '\n'
                  << "Saved " << directory << '\n';
    } catch(const std::exception& e) { std::cerr << "ERROR: " << e.what() << '\n'; return 1; }
}
