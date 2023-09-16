#include <RcppEigen.h>
#include "basics.h"

// [[Rcpp::depends(RcppEigen)]]

//
// [[Rcpp::export]]
double getES(const Eigen::VectorXd & beta, const Eigen::VectorXd & tau, 
             const Eigen::VectorXd & gamma, const Eigen::VectorXd & alpha, 
             const Eigen::MatrixXd & Sig, 
             Eigen::MatrixXd & Z, Eigen::MatrixXd & X1, 
             Eigen::MatrixXd & W, const Eigen::VectorXd & Y, 
             Eigen::VectorXd & X2, 
             const Eigen::VectorXd & obstime,
             const Eigen::MatrixXd & xsmatrix, 
             const Eigen::MatrixXd & wsmatrix,
             const Eigen::VectorXd & pSLR, 
             const Eigen::VectorXd & Si,
             const Eigen::VectorXd & CH0s,
             const Eigen::VectorXd & CH0u){ 
  
  //calculate the square root of random effect covariance matrix 
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(Sig, Eigen::ComputeThinU | Eigen::ComputeThinV);
  Eigen::VectorXd eigenSQ = svd.singularValues();
  int i,j,q,t,db,uu;
  for (i=0;i<eigenSQ.size();i++) {
    eigenSQ(i) = sqrt(eigenSQ(i));
  }
  Eigen::MatrixXd SigSQRT  = svd.matrixU() * eigenSQ.asDiagonal() * svd.matrixV().transpose();
  
  int p1a=Z.cols();
  int ni = Z.rows();
  double xgamma1,temp,mu,sigma,zb,wi;
  Eigen::VectorXd bwi(p1a+1);
  Eigen::VectorXd bi(p1a);
  Eigen::VectorXd weightbwi(p1a+1);
  Eigen::VectorXd ri(p1a+1);
  Eigen::VectorXd obstimet(ni);
  
  int point=wsmatrix.rows();
  
  double alphaw;
  Eigen::VectorXd alphab(p1a);
  for (i=0;i<p1a;i++) alphab(i) = alpha(i);
  alphaw = alpha(p1a);
  
  int nt=pSLR.size();
  
  double S=0;
  double dem=0;
  for (t=0;t<nt;t++) {
    
    // fill in interval-censored covariate values
    for (i=0;i<ni;i++) {
      obstimet(i) = obstime(i) - Si(t);
      X1(i,1) = Si(t);
      W(i,1) = Si(t);
    }
    X1.col(2) = obstimet;
    W.col(2) = obstimet;
    X2(0) = Si(t);
    xgamma1=MultVV(X2,gamma);
    
    for (db=0;db<point;db++) {
      bwi = xsmatrix.row(db);
      weightbwi = wsmatrix.row(db);
      ri = sqrt(2)*SigSQRT*bwi;
      temp=exp(10);
      
      for (i=0;i<p1a;i++) bi(i)=ri(i);
      wi=ri(p1a);
      
      for (i=0;i<ni;i++) {
        mu=MultVV(X1.row(i),beta);
        zb=MultVV(Z.row(i),bi);
        sigma=exp(MultVV(W.row(i),tau) + wi);
        temp*=1/sqrt(sigma)*exp(-1/(2*sigma)*pow((Y(i) - mu - zb), 2));
      }
      
      temp*=exp(0-CH0s(t)*exp(xgamma1+MultVV(alphab,bi)+alphaw*wi));
      for (i=0;i<(p1a+1);i++) temp*=weightbwi(i);
      temp*=pSLR(t);
      
      dem+=temp;
      
      //calculate survival
      S+=exp(log(temp) + CH0s(t)*exp(xgamma1+MultVV(alphab,bi)+alphaw*wi) - 
        CH0u(t)*exp(xgamma1+MultVV(alphab,bi)+alphaw*wi));
    }
  }
  
  
  if(dem==0) {
    Rprintf("Program stops because of the data issue.\n");
    return ( 100.0 );
  }
  
  S/=dem;
  
  return S;
}