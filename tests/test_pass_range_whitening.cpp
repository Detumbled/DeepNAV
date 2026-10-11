#include "observations/radiometric/PassRangeWhitening.hpp"
#include <iostream>
#include <stdexcept>
using fd::observations::radiometric::PassRangeWhitening;
int main() {try {
    Eigen::VectorXd sigma(7);sigma<<2,3,2,4,2,5,2;
    const std::vector<std::vector<Eigen::Index>> groups{{0,2,4},{6}};
    const double bias=10;
    PassRangeWhitening whitening(sigma,groups,bias);
    Eigen::MatrixXd R=sigma.array().square().matrix().asDiagonal();
    for(const auto& g:groups) for(auto i:g) for(auto j:g) R(i,j)+=bias*bias;
    const Eigen::MatrixXd W=whitening.matrix(Eigen::MatrixXd::Identity(7,7));
    if((W.transpose()*W-R.inverse()).norm()>1e-12) throw std::runtime_error("Whitening disagrees with dense covariance inverse");
    Eigen::MatrixXd H=Eigen::MatrixXd::Random(7,3);Eigen::VectorXd residual=Eigen::VectorXd::LinSpaced(7,-2,4);
    const auto A=whitening.matrix(H);const auto r=whitening.vector(residual);
    if((A.transpose()*A-H.transpose()*R.inverse()*H).norm()>1e-12 ||
       (A.transpose()*r-H.transpose()*R.inverse()*residual).norm()>1e-12)
        throw std::runtime_error("Correlated normal equations disagree with dense reference");
    for(int n:{1,18,100}) {
        std::vector<Eigen::Index> rows;for(int i=0;i<n;++i) rows.push_back(i);
        PassRangeWhitening p(Eigen::VectorXd::Ones(n),{rows},10);
        const double information=p.vector(Eigen::VectorXd::Ones(n)).squaredNorm();
        if(std::abs(information-n/(1.+n*100.))>1e-12) throw std::runtime_error("Pass shared error averaged away");
    }
    bool rejected=false;try {PassRangeWhitening p(sigma,{{0,0}},bias);}catch(const std::invalid_argument&){rejected=true;}
    if(!rejected)throw std::runtime_error("Overlapping covariance rows accepted");
    std::cout<<"Pass covariance whitening and common-error floor verified.\n";return 0;
} catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
