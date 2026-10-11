#pragma once
#include <Eigen/Core>
#include <array>
#include <cmath>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <map>
#include <regex>
#include <stdexcept>
#include <string>

// Read the release constants paired with the actual reconstructed SPK.
// This analytic pole series uses degrees and Julian centuries of TDB from
// POLTIM. RA uses sine; declination uses cosine (Jacobson 2022, Eqs. 1–2).
// It is the SAT427 series in the label, not the later SAT441 coefficients.
class Sat427Gravity {
    std::map<std::string,double> constants_, masses_;
    double epochJD_, ra0_, raRate_, dec0_, decRate_;
    std::array<double,3> raAmplitude_, raPhase_, raFrequency_, decAmplitude_, decPhase_, decFrequency_;
public:
    explicit Sat427Gravity(const std::filesystem::path& label) {
        std::ifstream file(label);
        if (!file) throw std::runtime_error("Missing SAT427 SPK label: " + label.string());
        const std::string text((std::istreambuf_iterator<char>(file)), {});
        if (text.find("SAT427") == std::string::npos)
            throw std::runtime_error("Expected SAT427 release label");
        const std::regex constant(R"(\b([A-Z][A-Za-z0-9_]*)\s+([+-]?[0-9]+\.[0-9]+E[+-][0-9]+))");
        for (auto it=std::sregex_iterator(text.begin(),text.end(),constant); it!=std::sregex_iterator(); ++it)
            constants_[(*it)[1].str()] = std::stod((*it)[2].str());
        const std::regex mass(R"(\b(Mimas|Enceladus|Tethys|Dione|Rhea|Titan|Hyperion|Iapetus|Phoebe|Helene|Saturn)\s+[0-9]+\s+([0-9]+\.[0-9]+E[+-][0-9]+))");
        for (auto it=std::sregex_iterator(text.begin(),text.end(),mass); it!=std::sregex_iterator(); ++it) {
            std::string name=(*it)[1].str();
            for (auto& c:name) c=static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
            masses_[name]=std::stod((*it)[2].str());
        }
        epochJD_=value("POLTIM"); ra0_=value("ZACPL6"); raRate_=value("DACPL6");
        dec0_=value("ZDEPL6"); decRate_=value("DDEPL6"); (void)massGM("SATURN");
        for (int i=0;i<3;++i) {
            const auto suffix=std::to_string(i+1);
            raAmplitude_[i]=value("RA_AMPL"+suffix); raPhase_[i]=value("RA_FAZE"+suffix); raFrequency_[i]=value("RA_FREQ"+suffix);
            decAmplitude_[i]=value("DE_AMPL"+suffix); decPhase_[i]=value("DE_FAZE"+suffix); decFrequency_[i]=value("DE_FREQ"+suffix);
        }
    }
    double value(const std::string& name) const { return constants_.at(name); }
    double massGM(const std::string& name) const { return masses_.at(name); }
    double resolvedGM(const std::string& name, double fallback) const {
        if (masses_.count(name)) return massGM(name);
        if (name=="5" || name=="7" || name=="8") return value("GM"+name);
        return fallback;
    }
    Eigen::Vector3d pole(double et) const {
        constexpr double rad=3.14159265358979323846/180;
        const double T=(et/86400 + (2451545.0-epochJD_))/36525;
        double ra=ra0_+T*raRate_, dec=dec0_+T*decRate_;
        for (int i=0;i<3;++i) {
            ra += raAmplitude_[i]*std::sin(rad*(raPhase_[i]+T*raFrequency_[i]));
            dec += decAmplitude_[i]*std::cos(rad*(decPhase_[i]+T*decFrequency_[i]));
        }
        return {std::cos(dec*rad)*std::cos(ra*rad),std::cos(dec*rad)*std::sin(ra*rad),std::sin(dec*rad)};
    }
};
