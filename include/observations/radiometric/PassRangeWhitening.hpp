#pragma once
#include <Eigen/Dense>
#include <vector>
#include <stdexcept>
#include <cmath>

namespace fd::observations::radiometric {
// R_pass = sigma_independent^2 I + sigma_bias^2 11^T. Symmetric
// inverse-square-root whitening marginalizes an independent Gaussian bias
// per pass without adding bias parameters to the estimated state.
class PassRangeWhitening {
    Eigen::VectorXd sigma_;
    std::vector<std::vector<Eigen::Index>> groups_;
    std::vector<double> meanFactors_;
public:
    PassRangeWhitening(Eigen::VectorXd sigma,
        std::vector<std::vector<Eigen::Index>> groups,double biasSigma)
        :sigma_(std::move(sigma)),groups_(std::move(groups)) {
        if (sigma_.size()==0 || !sigma_.allFinite() || (sigma_.array()<=0).any() ||
            !std::isfinite(biasSigma) || biasSigma<0)
            throw std::invalid_argument("Invalid pass covariance scales");
        std::vector<bool> used(sigma_.size(),false);
        for (const auto& group:groups_) {
            if (group.empty()) throw std::invalid_argument("Empty range pass");
            double s=0;
            for (auto row:group) {
                if (row<0 || row>=sigma_.size() || used[row]) throw std::invalid_argument("Invalid/overlapping range pass rows");
                used[row]=true;
                if (s && sigma_[row]!=s) throw std::invalid_argument("Range samples within a pass require a common independent sigma");
                s=sigma_[row];
            }
            meanFactors_.push_back(s/std::sqrt(s*s+group.size()*biasSigma*biasSigma)-1);
        }
    }
    Eigen::MatrixXd matrix(const Eigen::MatrixXd& input) const {
        if (input.rows()!=sigma_.size() || !input.allFinite()) throw std::invalid_argument("Invalid whitening input");
        Eigen::MatrixXd out=input.array().colwise()/sigma_.array();
        for (std::size_t j=0;j<groups_.size();++j) {
            Eigen::RowVectorXd mean=Eigen::RowVectorXd::Zero(input.cols());
            for (auto row:groups_[j]) mean+=out.row(row);
            mean*=meanFactors_[j]/groups_[j].size();
            for (auto row:groups_[j]) out.row(row)+=mean;
        }
        return out;
    }
    Eigen::VectorXd vector(const Eigen::VectorXd& input) const {return matrix(Eigen::MatrixXd(input)).col(0);}
};
}
