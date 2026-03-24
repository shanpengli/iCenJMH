// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

using namespace Rcpp;

// [[Rcpp::export]]
double logLik_learn(
    const Eigen::VectorXd& bw,
    const Eigen::VectorXd& Y,
    const Eigen::MatrixXd& X,
    const Eigen::MatrixXd& Z,
    const Eigen::MatrixXd& W,
    const Eigen::MatrixXd& X2,
    const double CH0,
    const double HAZ0,
    const Eigen::VectorXd& beta,
    const Eigen::VectorXd& tau,
    const Eigen::VectorXd& gamma,
    const Eigen::VectorXd& alpha,
    const Eigen::MatrixXd& SigInv,
    const double logdetSig,
    const int status) {
  
  const int p1a = bw.size() - 1;
  
  Eigen::VectorXd b = bw.head(p1a);
  const double w = bw(p1a);
  
  Eigen::VectorXd sigma = (W * tau).array() + w;
  sigma = sigma.array().exp();
  
  Eigen::VectorXd resid = Y - X * beta - Z * b;
  
  const double linpred_surv =
    (X2 * gamma)(0) + alpha.dot(bw);
  
  double total =
    ((resid.array().square()) / (2.0 * sigma.array()) +
    0.5 * sigma.array().log()).sum() +
    CH0 * std::exp(linpred_surv) +
    0.5 * bw.dot(SigInv * bw) +
    0.5 * logdetSig;
  
  if (status == 1) {
    total -= std::log(HAZ0) + linpred_surv;
  }
  
  return total;
}