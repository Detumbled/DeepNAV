#include "DP853Integrator.hpp"
#include "CassiniAttitude.hpp"
#include "CassiniForceModel.hpp"
#include "Sat427Gravity.hpp"
#include "perturbations/ThermalRecoil.hpp"
#include "perturbations/FlatPlateSRP.hpp"
#include "dynamics/SpiceKernelSet.hpp"
#include "perturbations/Gravitational.hpp"
#include "perturbations/J2.hpp"
#include "perturbations/AnnularGravity.hpp"
#include "perturbations/Relativistic.hpp"
#include "utils/CSPICE/SpiceError.hpp"
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

using State = od::DP853Integrator::State6;
namespace {
State state(const std::string& body, double t) {
    double x[6], lt;
    spkezr_c(body.c_str(), t, "J2000", "NONE", "SATURN", x, &lt);
    od::throwIfSpiceFailed("State of " + body);
    return Eigen::Map<State>(x);
}
double gm(const std::string& body) {
    double mu; SpiceInt n;
    bodvrd_c(body.c_str(), "GM", 1, &n, &mu);
    od::throwIfSpiceFailed("GM of " + body);
    return mu;
}
double epoch(const std::string& utc) {
    double t; str2et_c(utc.c_str(), &t);
    od::throwIfSpiceFailed("UTC epoch"); return t;
}
struct Model {
    std::string name; bool refined;
    State x; double step{30};
};
}
int main(int argc, char** argv) {
  try {
    std::filesystem::path root(DEEPNAV_SOURCE_DIR);
    auto output = root / "results/cassini_trajectory";
    std::string utc = "2004-10-09T00:00:00";
    double days = 14, cadence = 600, area = 12, mass = 3220, cr = 1.42;
    // Bounded first-seven-day SPK calibration; see RECOIL_UPDATE.md.
    // Reference acceleration at 2781 kg, not an October thermal measurement.
    double thermalX = -2.725e-13, thermalY = -2.725e-13, thermalZ = -5.232237415805135e-12;
    bool gravityAudit=false, studyOnly=false;
    // Smallest position RMS among all eight combinations on the retained arc.
    // Bits: 1 indirect J2; 2 rings; 4 point-mass 1PN. Direct J2 is always on.
    constexpr int retainedMask=2;
    int gravityMask=retainedMask;
    for (int i=1; i<argc; ++i) {
      const std::string key(argv[i]);
      if (key == "--gravity-audit") { gravityAudit=true; continue; }
      if (key == "--study-only") { studyOnly=true; continue; }
      if (key == "--help") {
        std::cout << "--output DIR --start-utc UTC --days D --cadence-s S --area-m2 A --mass-kg M --cr CR --thermal-x-km-s2 AX --thermal-y-km-s2 AY --thermal-z-km-s2 AZ --gravity-audit\n"
                  << "--gravity-mask 0..7 (bits: 1 indirect J2, 2 rings, 4 relativity); --study-only skips numerical refinement for temporary ablation runs\n";
        return 0;
      }
      if (++i == argc) throw std::invalid_argument("Missing option value");
      if (key == "--gravity-mask") {
        const std::string value=argv[i];
        std::size_t parsed=0;
        gravityMask=std::stoi(value,&parsed);
        if (parsed!=value.size() || gravityMask<0 || gravityMask>7)
          throw std::invalid_argument("Gravity mask must be an integer from 0 to 7");
      }
      else if (key == "--output") output = argv[i];
      else if (key == "--start-utc") utc = argv[i];
      else {
        const double value = std::stod(argv[i]);
        if (!std::isfinite(value) || (key.rfind("--thermal-",0) != 0 && value <= 0)) throw std::invalid_argument("Options must be positive and finite");
        if (key == "--days") days = value;
        else if (key == "--cadence-s") cadence = value;
        else if (key == "--area-m2") area = value;
        else if (key == "--mass-kg") mass = value;
        else if (key == "--cr") cr = value;
        else if (key == "--thermal-x-km-s2") thermalX = value;
        else if (key == "--thermal-y-km-s2") thermalY = value;
        else if (key == "--thermal-z-km-s2") thermalZ = value;
        else throw std::invalid_argument("Unknown option: " + key);
      }
    }
    od::SpiceErrorModeGuard errors("RETURN", "NONE");
    std::vector<std::filesystem::path> paths;
    for (const char* name : {"naif0012.tls", "pck00011.tpc", "gm_de440.tpc", "de442.bsp",
          "Cassini/spk/200128RU_SCPSE_04251_04303.bsp"}) paths.push_back(root / "Kernels" / name);
    for (const char* name : {"Cassini/fk/cas_v40.tf", "Cassini/sclk/cas00172.tsc",
          "Cassini/ck/04256_04292py_as_flown.bc", "Cassini/ck/04292_04320py_as_flown.bc",
          "Cassini/ck/04281_04286ra.bc", "Cassini/ck/04286_04291ra.bc",
          "Cassini/ck/04291_04296ra.bc", "Cassini/ck/04296_04301ra.bc",
          "Cassini/ck/04293_04298ra.bc"}) paths.push_back(root / "Kernels" / name);
    fd::dynamics::SpiceKernelSet kernels(paths);
    const CassiniAttitude attitude(paths);
    const auto gravityLabel=root / "Kernels/Cassini/spk/200128RU_SCPSE_04251_04303.bsp.lbl";
    const Sat427Gravity gravity(gravityLabel);
    const double t0 = epoch(utc), duration = days*86400, mu = gravity.massGM("SATURN");
    if (gravityAudit) {
        std::filesystem::create_directories(output);
        std::ofstream audit(output / "gravity_audit.csv");
        audit.exceptions(std::ios::failbit | std::ios::badbit);
        audit << "epoch_tdb,px,py,pz\n" << std::setprecision(17);
        for (double et : {0.0, -36525*86400.0, 36525*86400.0, t0, t0+duration/2, t0+duration}) {
            const auto pole=gravity.pole(et);
            audit << et << "," << pole.x() << "," << pole.y() << "," << pole.z() << "\n";
        }
        return 0;
    }
    const double maneuver = epoch("2004-10-23T06:16:00"); // Earlier OTM-004 time in NASA PDS mass table; conservatively bracket event-report 07:29:28.
    // Validate coverage, including every gap midpoint, before writing results.
    (void)attitude.rotation(t0); (void)attitude.rotation(t0+duration);
    for (std::size_t i=1;i<attitude.coverage.size();++i) {
        const double left=attitude.coverage[i-1].second, right=attitude.coverage[i].first;
        if (left < t0+duration && right > t0) (void)attitude.rotation((left+right)/2);
    }
    const double thermalEpoch = epoch("2005-08-13T00:00:00");
    const fd::perturbations::ThermalRecoil thermal({thermalX,thermalY,thermalZ},thermalEpoch,2781);
    if (cr < 1 || cr > 5.0/3) throw std::invalid_argument("Diffuse HGA plate requires CR between 1 and 5/3");
    const double nu = (cr-1)/2;
    const State initial = state("CASSINI", t0);
    (void)state("CASSINI", t0+duration);
    const auto bodies=fd::simulation::cassiniPerturbers(gravity);
    for (const auto& body:bodies) { (void)state(body.first,t0); (void)state(body.first,t0+duration); }
    double radii[3]; SpiceInt n;
    bodvrd_c("SATURN", "RADII", 3, &n, radii);
    od::throwIfSpiceFailed("Saturn radii");
    const fd::simulation::CassiniForceModel force(gravity,attitude,bodies,gravityMask,mass,area,cr,
        {thermalX,thermalY,thermalZ},thermalEpoch,radii[0]);
    const double ringGM=gravity.value("RINGM1")+gravity.value("RINGM2")+gravity.value("RINGM4");
    const auto gravityTerms = [&](double et,const Eigen::Vector3d& r,const Eigen::Vector3d& v) {
        return force.gravityTerms(et,r,v,true);
    };
    std::vector<Model> models{
      {"multibody_j2_srp_thermal", false, initial},
      {"refined_full_model", true, initial}
    };
    if (studyOnly) models.resize(1);
    std::filesystem::create_directories(output);
    std::ofstream csv(output / "trajectory.csv"), validation(output / "validation.csv");
    csv.exceptions(std::ios::failbit | std::ios::badbit);
    validation.exceptions(std::ios::failbit | std::ios::badbit);
    const std::string header = "elapsed_s,epoch_tdb,utc,model,x_km,y_km,z_km,vx_km_s,vy_km_s,vz_km_s,dx_km,dy_km,dz_km,position_error_km,velocity_error_m_s,srp_acceleration_m_s2,thermal_acceleration_m_s2,attitude_source\n";
    csv << header << std::setprecision(17);
    validation << header << std::setprecision(17);
    std::ofstream audit(output / "force_audit.csv");
    audit.exceptions(std::ios::failbit | std::ios::badbit);
    audit << "elapsed_s,indirect_j2_x_km_s2,indirect_j2_y_km_s2,indirect_j2_z_km_s2,ring_x_km_s2,ring_y_km_s2,ring_z_km_s2,relativistic_x_km_s2,relativistic_y_km_s2,relativistic_z_km_s2\n" << std::setprecision(17);
    double previous = 0;
    for (std::size_t sample=0;; ++sample) {
      const double elapsed = std::min(sample*cadence, duration), t = t0+elapsed;
      const State reference = state("CASSINI", t);
      const auto assessed=gravityTerms(t,reference.head<3>(),reference.tail<3>());
      audit << elapsed;
      for (const auto& a:{assessed.indirectJ2,assessed.ring,assessed.relativistic})
          for (int j=0;j<3;++j) audit << ',' << a[j];
      audit << '\n';
      char stamp[80]; et2utc_c(t, "ISOC", 3, sizeof(stamp), stamp);
      od::throwIfSpiceFailed("Output UTC");
      auto write = [&](const std::string& name, const State& x, double radiation, double recoil, bool refined = false) {
        auto& stream = refined ? validation : csv;
        const State delta = x-reference;
        stream << elapsed << ',' << t << ',' << stamp << ',' << name;
        for (int j=0; j<6; ++j) stream << ',' << x[j];
        for (int j=0; j<3; ++j) stream << ',' << delta[j];
        stream << ',' << delta.head<3>().norm() << ',' << delta.tail<3>().norm()*1000 << ',' << radiation << ',' << recoil << ',' << attitude.source(t) << '\n';
      };
      write("spk", reference, 0, 0);
      for (auto& model : models) {
        auto dynamics = [&](double dt, State x) -> State {
          const double et = t0+dt;
          const Eigen::Vector3d r=x.head<3>(), v=x.tail<3>();
          const Eigen::Vector3d a=force.evaluate(et,r,v).acceleration;
          State derivative; derivative << v, a; return derivative;
        };
        od::DP853Integrator::Options options;
        options.absoluteTolerance = model.refined ? 1e-13 : 1e-12;
        options.relativeTolerance = model.refined ? 1e-14 : 1e-12;
        options.maximumStep = model.refined ? 30 : 120;
        if (sample) model.x = od::DP853Integrator(options).integrate(previous, model.x, elapsed, dynamics, model.step);
        const double radiation = 1000*fd::perturbations::flatPlateSRPWithShadow(
            model.x.head<3>(),state("SUN",t).head<3>(),-attitude.rotation(t).col(2),area,mass,nu,radii[0],695700).acceleration.norm();
        write(model.name,model.x,radiation,1000*thermal.evaluate(t,mass,attitude.rotation(t)).acceleration.norm(),model.refined);
      }
      previous = elapsed;
      if (sample%500 == 0) std::cout << "Sample " << sample << ": day " << elapsed/86400 << '\n';
      if (elapsed == duration) break;
    }
    std::ofstream meta(output / "metadata.json");
    meta.exceptions(std::ios::failbit | std::ios::badbit);
    meta << std::setprecision(17) << "{\n\"start_utc\":\"" << utc << "\",\n\"duration_days\":" << days
      << ",\n\"cadence_s\":" << cadence << ",\n\"area_m2\":" << area << ",\n\"mass_kg\":" << mass
      << ",\n\"cr\":" << cr << ",\n\"saturn_gm_km3_s2\":" << mu << ",\n\"saturn_radius_km\":" << gravity.value("RADIUS")
      << ",\n\"saturn_shadow_radius_km\":" << radii[0]
      << ",\n\"gravity_release\":\"SAT427.23\",\n\"gravity_label_path\":" << std::quoted(gravityLabel.string())
      << ",\n\"pole_model\":\"SAT427 label analytic series; TDB Julian centuries; RA sine, DEC cosine\""
      << ",\n\"gravity_mask\":" << gravityMask
      << ",\n\"retained_default_gravity_mask\":" << retainedMask
      << ",\n\"numerical_validation\":" << (studyOnly ? "false" : "true")
      << ",\n\"indirect_j2\":" << ((gravityMask&1) ? "true" : "false")
      << ",\n\"rings_enabled\":" << ((gravityMask&2) ? "true" : "false")
      << ",\n\"relativity_enabled\":" << ((gravityMask&4) ? "true" : "false")
      << ",\n\"ring_total_gm_km3_s2\":" << ringGM
      << ",\n\"ring_model\":\"Uniform concentric C/B/A annuli; exterior degree 80; composite-origin shape reaction\""
      << ",\n\"relativity\":\"Barycentric point-mass EIH, beta=gamma=1; spacecraft minus Saturn; 20 sources\""
      << ",\n\"thermal_reference_x_km_s2\":" << thermalX
      << ",\n\"thermal_reference_y_km_s2\":" << thermalY
      << ",\n\"attitude_reconstructed_seconds\":" << CassiniAttitude::secondsBetween(attitude.reconstructed,t0,t0+duration)
      << ",\n\"thermal_reference_z_km_s2\":" << thermalZ
      << ",\n\"thermal_reference_mass_kg\":2781,\n\"thermal_reference_utc\":\"2005-08-13T00:00:00\",\n\"thermal_half_life_years\":87.74"
      << ",\n\"srp_diffuse_nu\":" << nu
      << ",\n\"j2\":" << gravity.value("J602") << ",\n\"maneuver_elapsed_s\":" << maneuver-t0 << ",\n\"bodies\":[";
    for (std::size_t i=0;i<bodies.size();++i) meta << (i ? "," : "") << "{\"name\":\"" << bodies[i].first << "\",\"gm_km3_s2\":" << bodies[i].second << '}';
    meta << "],\n\"kernel_paths\":[";
    for (std::size_t i=0;i<paths.size();++i) meta << (i ? "," : "") << std::quoted(paths[i].string());
    meta << "]\n}\n";
    std::cout << "Wrote " << output << '\n';
  } catch (const std::exception& error) { std::cerr << error.what() << '\n'; return 1; }
}
