#pragma once

#include "opnav/Types.hpp"

namespace fd::opnav {

        class DistortionModel {
        public:
            virtual ~DistortionModel() = default;

            [[nodiscard]]
            virtual NormalizedCoordinates distort(
                const NormalizedCoordinates& undistorted) const = 0;
        };

        class NoDistortion final : public DistortionModel {
        public:
            [[nodiscard]]
            NormalizedCoordinates distort(
                const NormalizedCoordinates& ideal) const override {
                return ideal;
            }
        };


    /**
     * 
     * Implementation of the simple OpenCV distorion model
     * 
     */

    
    /**
     * 
     * OpenCV distortion model coefficients, this should be selected after calibraiton
     * 
     * 
     */
    struct OpenCVDistortionCoefficients{
        //Radial distortion coefficients
        double K1 {0.0};
        double K2 {0.0};
        double K3 {0.0};
        double K4 {0.0};
        double K5 {0.0};
        double K6 {0.0};

        //Tangential distortion coefficients
        double P1 {0.0};
        double P2 {0.0};




    };

    class OpenCVDistortionModel final : public DistortionModel {
        public:
            explicit OpenCVDistortionModel(
                OpenCVDistortionCoefficients coefficients)
                : coefficients_(coefficients) {}

            [[nodiscard]]
            NormalizedCoordinates distort(
                const NormalizedCoordinates& undistorted) const override;

        private:
            OpenCVDistortionCoefficients coefficients_;
        };

}