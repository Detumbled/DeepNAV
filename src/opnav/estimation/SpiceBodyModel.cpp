#include "opnav/estimation/SpiceBodyModel.hpp"
#include "utils/CSPICE/SpiceError.hpp"
#include <cmath>
namespace fd::opnav {
BodyModel makeSpiceBodyModel(const std::string& body,TdbEpoch epoch,const CameraAttitude& camera,const CameraIntrinsics& intrinsics) {
    if(!std::isfinite(epoch.secondsPastJ2000))throw std::invalid_argument("Invalid body orientation epoch.");
    const Eigen::Matrix2d k=intrinsics.focalPlaneToPixelMatrix_PixelPermm;
    const double scale=std::sqrt(std::abs(k.determinant()));
    if(!(scale>0) || (k*k.transpose()/std::pow(scale,2)-Eigen::Matrix2d::Identity()).norm()>1e-10)
        throw std::invalid_argument("Body renderer requires square pixels without focal-plane shear.");
    od::SpiceErrorModeGuard guard;SpiceDouble radii[3],rotation[3][3];SpiceInt n,frame;SpiceChar name[128];SpiceBoolean found;
    bodvrd_c(body.c_str(),"RADII",3,&n,radii);cnmfrm_c(body.c_str(),sizeof(name),&frame,name,&found);
    od::throwIfSpiceFailed("Body model registration");if(n!=3 || !found)throw std::invalid_argument("Missing body radii/frame.");
    pxform_c(name,"J2000",epoch.secondsPastJ2000,rotation);od::throwIfSpiceFailed("Body orientation");
    Eigen::Matrix3d inertial;for(int i=0;i<3;++i)for(int j=0;j<3;++j)inertial(i,j)=rotation[i][j];
    Eigen::Matrix3d view=Eigen::Matrix3d::Zero();view.topLeftCorner<2,2>()=k/scale;view(2,2)=-1;
    BodyModel model;model.body=body;model.axes=Eigen::Map<Eigen::Vector3d>(radii);model.axes/=model.axes.mean();
    model.bodyToView=view*camera.matrix()*inertial;return model;
}
}
