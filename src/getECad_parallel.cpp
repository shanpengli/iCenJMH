// [[Rcpp::depends(RcppEigen)]]
#include "basics.h"
#include <RcppEigen.h>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

// [[Rcpp::export]]
int getECad_parallel(
    const Eigen::VectorXd & beta, const Eigen::VectorXd & tau,
    const Eigen::VectorXd & gamma, const Eigen::VectorXd & alpha,
    const Eigen::MatrixXd & H0Y, const Eigen::MatrixXd & Sig,
    const Eigen::MatrixXd & X1, const Eigen::MatrixXd & Z,
    const Eigen::MatrixXd & W, const Eigen::VectorXd & Y,
    const Eigen::MatrixXd & X2, const Eigen::VectorXd & survtime,
    const Eigen::VectorXd & status, const Eigen::VectorXd & ni,
    const Eigen::VectorXd & nt, const Eigen::MatrixXd & xsmatrix,
    const Eigen::MatrixXd & wsmatrix, const Eigen::MatrixXd & pSLR,
    const Eigen::MatrixXd & Posbwi, const Eigen::MatrixXd & Poscov,
    Eigen::Map<Eigen::MatrixXd> & Psl, Eigen::Map<Eigen::VectorXd> & FUNENW,
    Eigen::Map<Eigen::MatrixXd> & FUNEBNW, Eigen::Map<Eigen::MatrixXd> & FUNEBSNW,
    Eigen::Map<Eigen::VectorXd> & FUNE, Eigen::Map<Eigen::MatrixXd> & FUNBW,
    Eigen::Map<Eigen::MatrixXd> & FUNBWE,
    Eigen::Map<Eigen::MatrixXd> & FUNBWSE, Eigen::Map<Eigen::MatrixXd> & FUNBWS,
    const double pStol,
    const int nthreads = 1) {
  
  // calculate the square root of inverse random effect covariance matrix
  Eigen::JacobiSVD<Eigen::MatrixXd> svdSig(
      Sig.inverse(), Eigen::ComputeThinU | Eigen::ComputeThinV);
  Eigen::VectorXd eigenSQ = svdSig.singularValues();
  for (int i = 0; i < eigenSQ.size(); i++) {
    eigenSQ(i) = std::sqrt(eigenSQ(i));
  }
  Eigen::MatrixXd SigSQRT =
    svdSig.matrixU() * eigenSQ.asDiagonal() * svdSig.matrixV().transpose();
  
  const int k = ni.rows();
  const int p1a = Z.cols();
  const int point = wsmatrix.rows();
  
  Eigen::VectorXd alphab(p1a);
  for (int i = 0; i < p1a; i++) alphab(i) = alpha(i);
  const double alphaw = alpha(p1a);
  
  // Starting row indices matching the original serial counti/countt logic:
  // - surv rows advance by nt(j)
  // - obs rows advance by ni(j) * nt(j)
  std::vector<int> obs_start(k), surv_start(k);
  obs_start[0] = 0;
  surv_start[0] = 0;
  for (int j = 1; j < k; ++j) {
    obs_start[j] =
      obs_start[j - 1] + static_cast<int>(ni(j - 1) * nt(j - 1));
    surv_start[j] =
      surv_start[j - 1] + static_cast<int>(nt(j - 1));
  }
  
  int error_flag = 0;
  int error_subject = -1;
  
#ifdef _OPENMP
  omp_set_num_threads(nthreads);
#pragma omp parallel for schedule(static)
#endif
  for (int j = 0; j < k; j++) {
    if (error_flag) continue;
    
    const int nij = static_cast<int>(ni(j));
    const int ntj = static_cast<int>(nt(j));
    const int surv0 = surv_start[j];
    const int subj_obs0 = obs_start[j];
    
    double demt = 0.0;
    
    for (int t = 0; t < ntj; t++) {
      const int rowt = surv0 + t;
      const int obs_t0 = subj_obs0 + t * nij;  // critical correction
      double dem = 0.0;
      
      if (pSLR(rowt, 3) <= pStol || std::isnan(pSLR(rowt, 3))) {
        Psl(rowt, 2) = 0.0;
        continue;
      }
      
      const double cuh0 = H0Y(rowt, 2);
      const double haz0 = H0Y(rowt, 3);
      const double xgamma = MultVV(X2.row(rowt), gamma);
      
      Eigen::MatrixXd Hi(p1a + 1, p1a + 1);
      for (int i = 0; i < (p1a + 1); i++) {
        Hi.row(i) = Poscov.row(rowt * (p1a + 1) + i);
      }
      
      Eigen::JacobiSVD<Eigen::MatrixXd> svdHi(
          Hi, Eigen::ComputeThinU | Eigen::ComputeThinV);
      Eigen::VectorXd eigenSQ2 = svdHi.singularValues();
      for (int i = 0; i < eigenSQ2.size(); i++) {
        eigenSQ2(i) = std::sqrt(eigenSQ2(i));
      }
      Eigen::MatrixXd Hi2 =
        svdHi.matrixU() * eigenSQ2.asDiagonal() * svdHi.matrixV().transpose();
      
      for (int db = 0; db < point; db++) {
        Eigen::VectorXd bwi = xsmatrix.row(db);
        Eigen::VectorXd weightbwi = wsmatrix.row(db);
        Eigen::VectorXd bwii = Posbwi.row(rowt);
        Eigen::VectorXd ri = bwii + std::sqrt(2.0) * Hi2 * bwi;
        Eigen::VectorXd rii = SigSQRT * ri;
        
        Eigen::VectorXd bi = ri.head(p1a);
        const double wi = ri(p1a);
        
        double temp = std::exp(10.0);
        
        for (int i = 0; i < nij; i++) {
          const int rowi = obs_t0 + i;
          const double mu = MultVV(X1.row(rowi), beta);
          const double zb = MultVV(Z.row(rowi), bi);
          const double sigma = std::exp(MultVV(W.row(rowi), tau) + wi);
          temp *= 1.0 / std::sqrt(sigma) *
            std::exp(-1.0 / (2.0 * sigma) * std::pow(Y(rowi) - mu - zb, 2.0));
        }
        
        if (status(rowt) == 1) {
          temp *= haz0 * std::exp(xgamma + MultVV(alpha, ri));
        }
        
        temp *= std::exp(-cuh0 * std::exp(xgamma + MultVV(alpha, ri)));
        temp *= pSLR(rowt, 3);
        
        for (int i = 0; i < (p1a + 1); i++) temp *= weightbwi(i);
        
        Eigen::VectorXd bwi2 = xsmatrix.row(db);
        temp *= std::exp(-std::pow(rii.norm(), 2.0) / 2.0) *
          std::exp(std::pow(bwi2.norm(), 2.0));
        
        dem += temp;
        
        // calculate h(theta_i)
        FUNENW(rowt) += temp * std::exp(-wi);
        FUNEBNW.row(rowt) += temp * std::exp(-wi) * bi.transpose();
        
        for (int i = 0; i < p1a; i++) {
          FUNEBSNW(rowt, i) += temp * std::exp(-wi) * std::pow(bi(i), 2.0);
        }
        
        if (p1a == 2) {
          for (int i = 1; i < p1a; i++) {
            for (int q = 0; q < p1a - i; q++) {
              FUNEBSNW(rowt, p1a + q + (i - 1) * (p1a - 1)) +=
                temp * std::exp(-wi) * bi(q) * bi(q + i);
            }
          }
        }
        
        FUNE(rowt) += temp * std::exp(MultVV(ri, alpha));
        FUNBW.row(rowt) += temp * ri.transpose();
        FUNBWE.row(rowt) += temp * std::exp(MultVV(ri, alpha)) * ri.transpose();
        
        for (int i = 0; i < (p1a + 1); i++) {
          FUNBWSE(rowt, i) +=
            temp * std::exp(MultVV(ri, alpha)) * std::pow(ri(i), 2.0);
          FUNBWS(rowt, i) += temp * std::pow(ri(i), 2.0);
        }
        
        if (p1a < 3) {
          for (int i = 1; i < (p1a + 1); i++) {
            for (int q = 0; q < (p1a + 1 - i); q++) {
              FUNBWSE(rowt, p1a + 1 + q + (i - 1) * p1a) +=
                temp * std::exp(MultVV(ri, alpha)) * ri(q) * ri(q + i);
              FUNBWS(rowt, p1a + 1 + q + (i - 1) * p1a) +=
                temp * ri(q) * ri(q + i);
            }
          }
        }
      }
      
      demt += dem;
      Psl(rowt, 2) = dem;
    }
    
    if (demt == 0.0) {
#ifdef _OPENMP
#pragma omp critical
#endif
{
  if (!error_flag) {
    error_flag = 100;
    error_subject = j;
  }
}
      continue;
    }
    
    for (int t = 0; t < ntj; t++) {
      const int rowt = surv0 + t;
      Psl(rowt, 2) = std::exp(std::log(Psl(rowt, 2)) - std::log(demt));
      FUNENW(rowt) /= demt;
      FUNEBSNW.row(rowt) /= demt;
      FUNE(rowt) /= demt;
      FUNEBNW.row(rowt) /= demt;
      FUNBW.row(rowt) /= demt;
      FUNBWE.row(rowt) /= demt;
      FUNBWS.row(rowt) /= demt;
      FUNBWSE.row(rowt) /= demt;
    }
  }
  
  if (error_flag) {
    Rcpp::Rcout << "E step ran into issue for the "
                << error_subject << "th subject. Program stops.\n";
    return error_flag;
  }
  
  return 0;
}