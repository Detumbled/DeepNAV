#pragma once

#include<Eigen/Dense>

namespace fd::opnav {

struct TdbEpoch {
    double secondsPastJ2000;

    explicit constexpr TdbEpoch(double seconds)
        : secondsPastJ2000(seconds) {}


};


//Next geometric types for the camera itsfel:
struct PixelCoordinates {
    double sample {0.0};
    double line {0.0};
    

};

using FocalPlaneCoordinates = Eigen::Vector2d;


using NormalizedCoordinates = Eigen::Vector2d;

struct ImageSize {
    int width;
    int height;
};

struct CameraIntrinsics {
    double focalLength; // in mm

    PixelCoordinates principalPoint; // in pixels

    /**
     * Matrix to convert form focal plane coordinates to pixel coordinates offsets 
     * 
     *    s     s_0   K_x k_xy      x
     *      =       *           *
     *    l     l_0   k_yx K_y      y
     * 
     * 
     */
    Eigen::Matrix2d focalPlaneToPixelMatrix_PixelPermm;//Pixels/mm


   using pixelOffset_t = Eigen::Vector2d;


};
} // namespace fd::opnav