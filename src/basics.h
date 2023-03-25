#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

double  MultVV(const Eigen::VectorXd & x, const Eigen::VectorXd & y);

Eigen::MatrixXd MultVVoutprod(const Eigen::VectorXd & x);

double MultVVinprod(const Eigen::VectorXd & x);

double CH(const Eigen::MatrixXd & H, double t);

double HAZ(const Eigen::MatrixXd & H, double t);

Eigen::MatrixXd MultVV2outprod(const Eigen::VectorXd & x, const Eigen::VectorXd & y);

Eigen::MatrixXd updatePSLR(const Eigen::MatrixXd & YS, const Eigen::VectorXd & mdata,
                           const Eigen::VectorXd & mdataS, const Eigen::VectorXd & idsum);

Eigen::MatrixXd GetCov(const Eigen::MatrixXd & S);

Rcpp::List GetSE(const int nbeta, const int ntau, const int nSig, 
                 const int ngamma, const int nalpha, const Eigen::MatrixXd & Cov);
  