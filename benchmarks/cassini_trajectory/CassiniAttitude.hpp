#pragma once

#include "utils/CSPICE/SpiceError.hpp"
#include <Eigen/Geometry>
#include <algorithm>
#include <filesystem>
#include <stdexcept>
#include <vector>

// Kernels must already be loaded: as-flown first, reconstructed last.
// INTERVAL coverage is checked, rather than trusting each CK's overall span.
class CassiniAttitude {
public:
    using Interval = std::pair<double,double>;
    std::vector<Interval> coverage, reconstructed;
    explicit CassiniAttitude(const std::vector<std::filesystem::path>& kernels) {
        for (const auto& path : kernels) {
            if (path.extension() != ".bc") continue;
            SPICEDOUBLE_CELL(window, 20000);
            // CSPICE's cell macro declares static storage; CKCOV unions
            // into its output, so clear it before reading each individual CK.
            scard_c(0, &window);
            ckcov_c(path.c_str(), -82000, SPICEFALSE, "INTERVAL", 0, "TDB", &window);
            od::throwIfSpiceFailed("Cassini CK interval coverage");
            for (SpiceInt i=0; i<wncard_c(&window); ++i) {
                double a,b; wnfetd_c(&window,i,&a,&b);
                coverage.emplace_back(a,b);
                if (path.filename().string().ends_with("ra.bc")) reconstructed.emplace_back(a,b);
            }
        }
        merge(coverage); merge(reconstructed);
    }
    static void merge(std::vector<Interval>& intervals) {
        std::sort(intervals.begin(),intervals.end());
        std::vector<Interval> result;
        for (const auto& interval : intervals) {
            if (!result.empty() && interval.first <= result.back().second)
                result.back().second = std::max(result.back().second,interval.second);
            else result.push_back(interval);
        }
        intervals = std::move(result);
    }
    static bool contains(const std::vector<Interval>& intervals, double t) {
        auto it = std::upper_bound(intervals.begin(),intervals.end(),t,
            [](double value,const Interval& interval){return value < interval.first;});
        return it != intervals.begin() && t <= std::prev(it)->second;
    }
    static double secondsBetween(const std::vector<Interval>& intervals, double start, double end) {
        double seconds=0;
        for (const auto& interval : intervals)
            seconds += std::max(0.0,std::min(end,interval.second)-std::max(start,interval.first));
        return seconds;
    }
    static Eigen::Matrix3d direct(double t) {
        double matrix[3][3]; pxform_c("CASSINI_SC_COORD","J2000",t,matrix);
        od::throwIfSpiceFailed("Cassini spacecraft attitude");
        Eigen::Matrix3d result;
        for (int i=0;i<3;++i) for (int j=0;j<3;++j) result(i,j)=matrix[i][j];
        return result;
    }
    [[nodiscard]] Eigen::Matrix3d rotation(double t) const {
        if (contains(coverage,t)) return direct(t);
        auto it = std::upper_bound(coverage.begin(),coverage.end(),t,
            [](double value,const Interval& interval){return value < interval.first;});
        if (it == coverage.begin() || it == coverage.end())
            throw std::runtime_error("Cassini attitude epoch outside CK coverage");
        const double left=std::prev(it)->second, right=it->first;
        if (right-left > 10)
            throw std::runtime_error("Cassini attitude gap exceeds the explicit 10 s bridge limit");
        // Only bridge the documented short boundary gap. No extrapolation
        // across telemetry outages: those use the loaded as-flown history.
        return Eigen::Quaterniond(direct(left)).slerp((t-left)/(right-left),
            Eigen::Quaterniond(direct(right))).toRotationMatrix();
    }
    [[nodiscard]] const char* source(double t) const {
        return contains(reconstructed,t) ? "reconstructed" :
            contains(coverage,t) ? "as_flown" : "short_gap_slerp";
    }
};
