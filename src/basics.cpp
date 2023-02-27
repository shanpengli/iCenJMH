// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>
#include "basics.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppEigen so that the build process will know what to do
//
// [[Rcpp::depends(RcppEigen)]]

// simple example of creating two matrices and
// returning the result of an operatioon on them
//
// via the exports attribute we tell Rcpp to make this function
// available from R
//

// another simple example: outer product of a vector,
// returning a matrix
//
// [[Rcpp::export]]
double  MultVV(const Eigen::VectorXd & x, const Eigen::VectorXd & y) {
    double v = x.transpose() * y;
    return v;
}

// [[Rcpp::export]]
Eigen::MatrixXd MultVVoutprod(const Eigen::VectorXd & x) {
    Eigen::MatrixXd m = x * x.transpose();
    return m;
}

// [[Rcpp::export]]
Eigen::MatrixXd MultVV2outprod(const Eigen::VectorXd & x, const Eigen::VectorXd & y) {
    Eigen::MatrixXd m = x * y.transpose();
    return m;
}

// and the inner product returns a scalar
//
// [[Rcpp::export]]
double  MultVVinprod(const Eigen::VectorXd & x) {
    double v = x.transpose() * x;
    return v;
}

// [[Rcpp::export]]
double  CH(const Eigen::MatrixXd & H, double t) {

    int a = H.rows();
    int i;
    double ch;
    if (t < H(0,0)) ch=0;
    else {
        ch=0;
        i=0;
        do {
            ch+=H(i, 2);
            i+=1;

        } while (i<=a-1 && t>= H(i,0));
        }

    return (ch);

}

// [[Rcpp::export]]
double HAZ(const Eigen::MatrixXd & H, double t) {

    int a = H.rows();
    int i;
    double temp=0;
    for (i=0;i<a;i++) {
        if (t == H(i, 0)) temp = H(i,2);
        }

    return(temp);

    }

// [[Rcpp::export]]
Eigen::MatrixXd MultMM(const Eigen::MatrixXd & x, const Eigen::MatrixXd & y) {
    Eigen::MatrixXd v = x * y;
    return v;
}

// [[Rcpp::export]]
Eigen::MatrixXd updatePSLR(const Eigen::MatrixXd & YS, const Eigen::VectorXd & mdata,
                           const Eigen::VectorXd & mdataS, const Eigen::VectorXd & idsum) {
  
  Eigen::MatrixXd pSLR = YS;
  double phisu;
  int k=idsum.size();
  int i,j,q;
  for (j=0;j<k;j++)
  {
    phisu=0;
    q=mdata(j);
    for (i=0;i<q;i++)
    {
      pSLR(mdataS(j)-1+i, 3)=exp(-phisu)*(1-exp(-YS(mdataS(j)-1+i, 3)))/(1-exp(-idsum(j)));
      phisu += YS(mdataS(j)-1+i, 3);
    }
  }
  
  return pSLR;
}

// [[Rcpp::export]]
Eigen::MatrixXd GetCov(const Eigen::MatrixXd & S) {
  
  Eigen::MatrixXd SS = Eigen::MatrixXd::Zero(S.cols(), S.cols());
  Eigen::VectorXd N = Eigen::VectorXd::Zero(S.cols());
  // Eigen::VectorXd NT = Eigen::VectorXd::Zero(S.cols());
  int i;
  int k = S.rows();
  for (i=0;i<k;i++) {
    N = S.row(i);
    // NT += N;
    SS += MultVVoutprod(N);
  }
  return SS.inverse();
  // return (SS - MultVVoutprod(NT)/k).inverse();
}

// [[Rcpp::export]]
Eigen::MatrixXd GetCov2(const Eigen::MatrixXd & S) {
  
  Eigen::MatrixXd SS = Eigen::MatrixXd::Zero(S.cols(), S.cols());
  Eigen::VectorXd N = Eigen::VectorXd::Zero(S.cols());
  Eigen::VectorXd NT = Eigen::VectorXd::Zero(S.cols());
  int i;
  int k = S.rows();
  for (i=0;i<k;i++) {
    N = S.row(i);
    NT += N;
    //SS += MultVVoutprod(N);
  }
  return MultVVoutprod(NT)/k;
  // return (SS - MultVVoutprod(NT)/k).inverse();
}

// [[Rcpp::export]]
Rcpp::List GetSE(const int nbeta, const int ntau, const int nSig, 
                 const int ngamma, const int nalpha, const Eigen::MatrixXd & Cov) {
  
  Eigen::VectorXd sebeta = Eigen::VectorXd::Zero(nbeta);
  Eigen::VectorXd setau = Eigen::VectorXd::Zero(ntau);
  Eigen::VectorXd segamma = Eigen::VectorXd::Zero(ngamma);
  Eigen::VectorXd sealpha = Eigen::VectorXd::Zero(nalpha);
  Eigen::MatrixXd seSig = Eigen::MatrixXd::Zero(nSig, nSig);
  int t,q;
  for (t=0;t<nbeta;t++) sebeta(t) = sqrt(Cov(t,t));
  for (t=0;t<ntau;t++) setau(t) = sqrt(Cov(nbeta+t,nbeta+t));
  for (t=0;t<nSig;t++) seSig(t,t) = sqrt(Cov(nbeta+ntau+t,nbeta+ntau+t));
  for(q=1;q<nSig;q++)
  {
    for(t=0;t<(nSig-q);t++) {
      seSig(t,q+t) = sqrt(Cov(nbeta+ntau+nSig+t+(q-1)*(nSig-1),nbeta+ntau+nSig+t+(q-1)*(nSig-1)));
      seSig(q+t,t) = seSig(t,q+t);
    }
  }
  for (t=0;t<ngamma;t++) segamma(t) = sqrt(Cov(nbeta+ntau+(nSig+1)*nSig/2+t,
                                 nbeta+ntau+(nSig+1)*nSig/2+t));
  for (t=0;t<nalpha;t++) {
    sealpha(t) = sqrt(Cov(nbeta+ntau+(nSig+1)*nSig/2+ngamma+t,
                      nbeta+ntau+(nSig+1)*nSig/2+ngamma+t));
  }
  
  return Rcpp::List::create(Rcpp::Named("sebeta")=sebeta,
                            Rcpp::Named("setau")=setau,
                            Rcpp::Named("segamma")=segamma,
                            Rcpp::Named("sealpha")=sealpha,
                            Rcpp::Named("seSig")=seSig);
  
}


