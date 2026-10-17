// Compare production optical geometry with direct CSPICE, using independent CK pointing.
#include "dynamics/SpiceKernelSet.hpp"
#include "dynamics/SpiceStateProvider.hpp"
#include "opnav/CameraModel.hpp"
#include "opnav/core/ApparentDirection.hpp"
#include "opnav/image/CassiniImage.hpp"
#include "utils/CSPICE/SpiceError.hpp"
#include <SpiceUsr.h>
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <numbers>
#include <stdexcept>
#include <string>

namespace fs = std::filesystem;
using namespace fd::opnav;
namespace {
constexpr double arcsec = 180.0 * 3600 / std::numbers::pi;
void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}
std::ofstream output(const fs::path& path) {
    std::ofstream stream(path);
    if (!stream) throw std::runtime_error("Cannot write " + path.string());
    stream.exceptions(std::ios::badbit | std::ios::failbit);
    stream << std::setprecision(17);
    return stream;
}
double poolScalar(const char* name) {
    double value;
    SpiceInt count;
    SpiceBoolean found;
    gdpool_c(name, 0, 1, &count, &value, &found);
    od::throwIfSpiceFailed(name);
    require(found && count == 1 && std::isfinite(value), "Missing instrument parameter.");
    return value;
}
Eigen::Vector3d reference(double epoch, const char* correction, double& lightTime) {
    SpiceDouble vector[3];
    spkpos_c("ENCELADUS", epoch, "J2000", correction, "CASSINI", vector, &lightTime);
    od::throwIfSpiceFailed("Enceladus reference direction / retarded SPK coverage");
    return Eigen::Map<Eigen::Vector3d>(vector);
}
double angle(const Eigen::Vector3d& a, const Eigen::Vector3d& b) {
    const Eigen::Vector3d u = a.normalized(), v = b.normalized();
    return std::atan2(u.cross(v).norm(), u.dot(v));
}
}
int main(int argc, char** argv) {
    try {
        const fs::path root = DEEPNAV_SOURCE_DIR;
        fs::path input = root / "Input Optical/N1476124698_2.IMG";
        fs::path kernel = root / "Kernels/Cassini/cassini_2004_opnav.tm";
        fs::path directory = root / "results/cassini_optical_validation";
        for (int i = 1; i < argc; ++i) {
            const std::string option = argv[i];
            if (option == "--help") {
                std::cout << "cassini_optical_validation [--image IMG] [--kernels TM] [--output DIR]\n";
                return 0;
            }
            require(option == "--image" || option == "--kernels" || option == "--output", "Unknown option.");
            require(++i < argc, "Missing option value.");
            const fs::path path = fs::absolute(argv[i]);
            if (option == "--image") input = path;
            else if (option == "--kernels") kernel = path;
            else directory = path;
        }
        const auto image = image::loadCassiniImage(input);
        // This benchmark deliberately supports only the original full-frame NAC exposure.
        // Reject other modes instead of silently applying its calibration to them.
        auto label = input; label.replace_extension(".LBL");
        std::ifstream labels(label);
        const std::string text((std::istreambuf_iterator<char>(labels)), {});
        require(text.find("INSTRUMENT_ID = \"ISSNA\"") != std::string::npos &&
                text.find("INSTRUMENT_MODE_ID = \"FULL\"") != std::string::npos &&
                image.targetDescription == "ENCELADUS" && image.dn.rows() == 1024 && image.dn.cols() == 1024,
                "Requires a full-frame ISS NAC Enceladus exposure.");
        od::SpiceErrorModeGuard guard("RETURN", "NONE");
        // Standalone executable: start from an empty pool, with explicit load precedence.
        kclear_c();
        fs::current_path(root); // The local meta-kernel contains repository-relative paths.
        const fd::dynamics::SpiceKernelSet kernels({kernel});
        SpiceDouble epoch;
        str2et_c(image.imageMidTimeUtc.c_str(), &epoch);
        od::throwIfSpiceFailed("Image midpoint UTC to TDB");
        const double focal = poolScalar("INS-82360_FOCAL_LENGTH");
        const double pitch = poolScalar("INS-82360_PIXEL_SIZE") / 1000; // micrometers -> mm
        SpiceDouble center[2]; SpiceInt count; SpiceBoolean found;
        gdpool_c("INS-82360_FOV_CENTER_PIXEL", 0, 2, &count, center, &found);
        od::throwIfSpiceFailed("NAC zero-based principal point");
        require(found && count == 2 && focal > 0 && pitch > 0, "Invalid NAC calibration.");
        // IK: sample and line increase toward -X and -Y of CASSINI_ISS_NAC.
        const CameraModel camera({focal, {center[0], center[1]}, -Eigen::Matrix2d::Identity() / pitch});
        const fd::dynamics::SpiceStateProvider observer("CASSINI"), target("ENCELADUS");
        fs::create_directories(directory);
        auto rows = output(directory / "predictions.csv");
        rows << "offset_seconds,correction,sample,line,light_time_seconds\n";
        double maxCnAngle = 0, maxApparentPixels = 0, maxProjectionPixels = 0, maxLtError = 0;
        const double halfExposure = image.exposureMilliseconds / 2000;
        require(halfExposure > 0 && std::isfinite(halfExposure), "Invalid exposure duration.");
        // CK continuous interpolation coverage must contain the entire exposure,
        // not just its midpoint. Check all loaded CKs, without a pointing tolerance.
        SpiceInt kernelCount;
        ktotal_c("CK", &kernelCount);
        od::throwIfSpiceFailed("CK inventory");
        bool covered = false;
        SPICEDOUBLE_CELL(coverage, 20000);
        for (SpiceInt i = 0; i < kernelCount; ++i) {
            SpiceChar file[2048], type[32], source[2048]; SpiceInt handle;
            kdata_c(i, "CK", sizeof(file), sizeof(type), sizeof(source), file, type, source, &handle, &found);
            od::throwIfSpiceFailed("CK inventory entry");
            require(found, "Missing CK inventory entry.");
            scard_c(0, &coverage);
            ckcov_c(file, -82000, SPICEFALSE, "INTERVAL", 0, "TDB", &coverage);
            od::throwIfSpiceFailed("Reconstructed spacecraft CK interpolation coverage");
            covered = covered || wnincd_c(epoch-halfExposure, epoch+halfExposure, &coverage);
        }
        require(covered, "No continuous CK interval covers the full exposure.");
        for (double offset : {-halfExposure, 0.0, halfExposure}) {
            const double t = epoch + offset;
            SpiceDouble rotation[3][3];
            pxform_c("J2000", "CASSINI_ISS_NAC", t, rotation);
            od::throwIfSpiceFailed("NAC CK/FK pointing at exposure sample");
            Eigen::Matrix3d matrix;
            for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) matrix(i,j) = rotation[i][j];
            const CameraAttitude attitude(matrix);
            const auto ours = core::computeApparentDirection(TdbEpoch(t), observer.stateAt(TdbEpoch(t)), target);
            for (const char* correction : {"NONE", "CN", "CN+S"}) {
                double lt;
                const auto vector = reference(t, correction, lt);
                const auto pixel = camera.project(vector, attitude);
                // Independent projection with CSPICE rotation and explicit IK equation.
                SpiceDouble cameraVector[3];
                mxv_c(rotation, vector.data(), cameraVector);
                od::throwIfSpiceFailed("Reference camera rotation");
                const double ds = pixel.sample - (center[0] - focal/pitch*cameraVector[0]/cameraVector[2]);
                const double dl = pixel.line - (center[1] - focal/pitch*cameraVector[1]/cameraVector[2]);
                maxProjectionPixels = std::max(maxProjectionPixels, std::hypot(ds, dl));
                rows << offset << ',' << correction << ',' << pixel.sample << ',' << pixel.line << ',' << lt << '\n';
                if (std::string(correction) == "CN") {
                    maxCnAngle = std::max(maxCnAngle, angle(vector, ours.lightTime.geometricLineOfSightKm));
                    maxLtError = std::max(maxLtError, std::abs(lt-ours.lightTime.lightTimeSeconds));
                }
                if (std::string(correction) == "CN+S") {
                    const auto prediction = camera.project(ours.vectorKm, attitude);
                    maxApparentPixels = std::max(maxApparentPixels,
                        std::hypot(pixel.sample-prediction.sample, pixel.line-prediction.line));
                }
            }
            const auto prediction = camera.project(ours.vectorKm, attitude);
            rows << offset << ",DeepNAV," << prediction.sample << ',' << prediction.line << ','
                 << ours.lightTime.lightTimeSeconds << '\n';
        }
        // Numerical checks are separate from real-image residuals. Approximate
        // stellar aberration gets a declared 0.01 pixel budget, not exact equality.
        require(maxCnAngle < 1e-10 && maxLtError < 1e-7, "Light-time geometry disagrees with CSPICE CN.");
        require(maxProjectionPixels < 1e-9, "Camera projection disagrees with explicit IK equation.");
        require(maxApparentPixels < 0.01, "Approximate apparent direction exceeds 0.01 pixel budget.");
        auto metadata = output(directory / "geometry.json");
        metadata << "{\n\"image\": " << std::quoted(input.string())
                 << ",\n\"mid_time_utc\": " << std::quoted(image.imageMidTimeUtc)
                 << ",\n\"epoch_tdb_seconds\": " << epoch
                 << ",\n\"exposure_ms\": " << image.exposureMilliseconds
                 << ",\n\"focal_length_mm\": " << focal << ",\n\"pixel_pitch_mm\": " << pitch
                 << ",\n\"principal_point\": [" << center[0] << ',' << center[1] << ']'
                 << ",\n\"max_cn_angle_arcsec\": " << maxCnAngle*arcsec
                 << ",\n\"max_light_time_error_seconds\": " << maxLtError
                 << ",\n\"max_projection_error_pixels\": " << maxProjectionPixels
                 << ",\n\"max_apparent_error_pixels\": " << maxApparentPixels
                 << ",\n\"apparent_error_budget_pixels\": 0.01"
                 << ",\n\"continuous_ck_exposure_coverage\": true,\n\"distortion\": \"none; nominal IK pinhole calibration only\",\n\"loaded_kernels\": [";
        ktotal_c("ALL", &kernelCount);
        for (SpiceInt i = 0; i < kernelCount; ++i) {
            SpiceChar file[2048], type[32], source[2048]; SpiceInt handle;
            kdata_c(i, "ALL", sizeof(file), sizeof(type), sizeof(source), file, type, source, &handle, &found);
            od::throwIfSpiceFailed("Kernel provenance");
            require(found, "Missing kernel provenance entry.");
            metadata << (i ? "," : "") << std::quoted(fs::absolute(file).string());
        }
        metadata << "]\n}\n";
        std::cout << std::setprecision(10) << "PASS: continuous CK exposure coverage, retarded SPK queries, CN light time, camera projection.\n"
                  << "Maximum DeepNAV vs CN+S difference: " << maxApparentPixels << " pixels\n"
                  << "Saved " << directory << '\n';
    } catch (const std::exception& error) {
        std::cerr << "ERROR: " << error.what() << '\n';
        return 1;
    }
}
