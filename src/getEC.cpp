#include <RcppEigen.h>
#include "basics.h"

// [[Rcpp::depends(RcppEigen)]]

//
// [[Rcpp::export]]
int getEC(const Eigen::Map<Eigen::VectorXd> & beta, const Eigen::Map<Eigen::VectorXd> & tau, 
                 const Eigen::Map<Eigen::VectorXd> & gamma, const Eigen::Map<Eigen::VectorXd> & alpha,  
                 const Eigen::Map<Eigen::MatrixXd> & H0Y, const Eigen::Map<Eigen::MatrixXd> & Sig, 
                 const Eigen::Map<Eigen::MatrixXd> & X1, const Eigen::Map<Eigen::MatrixXd> & Z, 
                 const Eigen::Map<Eigen::MatrixXd> & W, const Eigen::Map<Eigen::VectorXd> & Y, 
                 const Eigen::Map<Eigen::MatrixXd> & X2, const Eigen::Map<Eigen::VectorXd> & survtime, 
                 const Eigen::Map<Eigen::VectorXd> & status, const Eigen::Map<Eigen::VectorXd> & ni, 
                 const Eigen::Map<Eigen::VectorXd> & nt, const Eigen::Map<Eigen::MatrixXd> & xsmatrix, 
                 const Eigen::Map<Eigen::MatrixXd> & wsmatrix, const Eigen::Map<Eigen::MatrixXd> & pSLR,
                 Eigen::Map<Eigen::MatrixXd> & Psl, Eigen::Map<Eigen::VectorXd> & FUNENW,
                 Eigen::Map<Eigen::MatrixXd> & FUNEBNW, Eigen::Map<Eigen::MatrixXd> & FUNEBSNW, 
                 Eigen::Map<Eigen::VectorXd> & FUNE, Eigen::Map<Eigen::MatrixXd> & FUNBW,
                 Eigen::Map<Eigen::MatrixXd> & FUNBWE, 
                 Eigen::Map<Eigen::MatrixXd> & FUNBWSE, Eigen::Map<Eigen::MatrixXd> & FUNBWS,
                 const double pStol){ 
  
  //calculate the square root of random effect covariance matrix 
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(Sig, Eigen::ComputeThinU | Eigen::ComputeThinV);
  Eigen::VectorXd eigenSQ = svd.singularValues();
  int i,j,q,t,db;
  for (i=0;i<eigenSQ.size();i++) {
    eigenSQ(i) = sqrt(eigenSQ(i));
  }
  Eigen::MatrixXd SigSQRT  = svd.matrixU() * eigenSQ.asDiagonal() * svd.matrixV().transpose();
  
  int k=ni.rows();
  int p1a=Z.cols();
  
  double dem,demt,cuh0,haz0,xgamma,temp,mu,sigma,zb,wi;
  Eigen::VectorXd bwi(p1a+1);
  Eigen::VectorXd bi(p1a);
  Eigen::VectorXd weightbwi(p1a+1);
  Eigen::VectorXd ri(p1a+1);
  Eigen::VectorXd alphab(p1a);
  double alphaw;
  for (i=0;i<p1a;i++) alphab(i) = alpha(i);
  alphaw = alpha(p1a);
    
  int point=wsmatrix.rows();
  int counti=0;
  int countt=0;
  for(j=0;j<k;j++)
  {
    demt=0;
    for (t=0;t<nt(j);t++)
    {
      dem=0;
      if (pSLR(countt, 3) <= pStol) {
        Psl(countt, 2) = dem;
        counti+=ni(j);
        countt++;
      } else {
        cuh0=H0Y(countt,2);
        haz0=H0Y(countt,3);
        xgamma=MultVV(X2.row(countt),gamma);

        for (db=0;db<point;db++)
        {
          bwi = xsmatrix.row(db);
          weightbwi = wsmatrix.row(db);
          ri = sqrt(2)*SigSQRT*bwi;
          temp=1;
          for (i=0;i<p1a;i++) bi(i)=ri(i);
          wi=ri(p1a);
          for (i=0;i<ni(j);i++)
          {
            mu=MultVV(X1.row(counti),beta);
            zb=MultVV(Z.row(counti),bi);
            sigma=exp(MultVV(W.row(counti),tau) + wi);
            temp*=1/sqrt(sigma)*exp(-1/(2*sigma)*pow((Y(counti) - mu - zb), 2));
            counti++;
          }
          if(status(countt)==1)  temp*=haz0*exp(xgamma+MultVV(alphab, bi)+alphaw*wi);
          temp*=exp(0-cuh0*exp(xgamma+MultVV(alphab,bi)+alphaw*wi));
          temp*=pSLR(countt, 3);
          for (i=0;i<(p1a+1);i++) temp*=weightbwi(i);
          dem+=temp;

          //calculate h(theta_i)
          FUNENW(countt)+=temp*exp(-wi);
          FUNEBNW.row(countt)+=temp*exp(-wi)*bi;

          for (i=0;i<p1a;i++) {
            FUNEBSNW(countt,i)+=temp*exp(-wi)*pow(bi(i),2);
          }
          if (p1a == 2) {
            for(i=1;i<p1a;i++)
            {
              for(q=0;q<p1a-i;q++) {
                FUNEBSNW(countt,p1a+q+(i-1)*(p1a-1))+=temp*exp(-wi)*bi(q)*bi(q+i);
              }
            }
          }

          FUNE(countt)+=temp*exp(MultVV(ri,alpha));
          FUNBW.row(countt)+=temp*ri;
          FUNBWE.row(countt)+=temp*exp(MultVV(ri,alpha))*ri;

          for (i=0;i<(p1a+1);i++) {
            FUNBWSE(countt,i)+=temp*exp(MultVV(ri,alpha))*pow(ri(i),2);
            FUNBWS(countt,i)+=temp*pow(ri(i),2);
          }

          if (p1a < 3) {
            for(i=1;i<(p1a+1);i++)
            {
              for(q=0;q<(p1a+1-i);q++)
              {
                FUNBWSE(countt,p1a+1+q+(i-1)*p1a)+=temp*exp(MultVV(ri,alpha))*ri(q)*ri(q+i);
                FUNBWS(countt,p1a+1+q+(i-1)*p1a)+=temp*ri(q)*ri(q+i);
              }
            }
          }

          if (db != point-1) counti = counti - ni(j);
        }
        demt+=dem;
        Psl(countt, 2) = dem;
        countt++;
      }
      
    }
    countt -= nt(j);
    if(demt==0) {
      Rprintf("E step ran into issue for the %dth subject. Program stops.\n", j);
      return ( 100.0 );
    }
    for (t=0;t<nt(j);t++)
    {
      Psl(countt, 2) /= demt;
      FUNENW(countt) /= demt;
      FUNEBSNW.row(countt) /= demt;
      FUNE(countt) /= demt;
      FUNEBNW.row(countt) /= demt;
      FUNBW.row(countt) /= demt;
      FUNBWE.row(countt) /= demt;
      FUNBWS.row(countt) /= demt;
      FUNBWSE.row(countt) /= demt;
      countt++;
    }

  }
  
  // Eigen::VectorXd TFUNENW = Eigen::VectorXd::Zero(k);
  // Eigen::MatrixXd TFUNEBNW = Eigen::MatrixXd::Zero(k,p1a);
  // Eigen::MatrixXd TFUNEBSNW = Eigen::MatrixXd::Zero(k,p1a*(p1a+1)/2);
  // Eigen::VectorXd TFUNE = Eigen::VectorXd::Zero(k);
  // Eigen::MatrixXd TFUNBW = Eigen::MatrixXd::Zero(k,(p1a+1));
  // Eigen::MatrixXd TFUNBWE = Eigen::MatrixXd::Zero(k,(p1a+1));
  // Eigen::MatrixXd TFUNBWSE = Eigen::MatrixXd::Zero(k,(p1a+2)*(p1a+1)/2);
  // Eigen::MatrixXd TFUNBWS = Eigen::MatrixXd::Zero(k,(p1a+2)*(p1a+1)/2);
  // countt=0;
  // 
  // for (j=0;j<k;j++) {
  //   
  //   for (t=0;t<nt(j);t++) {
  //     
  //     TFUNENW(j)+=FUNENW(countt);
  //     TFUNEBNW.row(j)+=FUNEBNW.row(countt);
  //     TFUNEBSNW.row(j)+=FUNEBSNW.row(countt);
  //     TFUNE(j)+=FUNE(countt);
  //     TFUNBW.row(j)+=FUNBW.row(countt);
  //     TFUNBWE.row(j)+=FUNBWE.row(countt);
  //     TFUNBWSE.row(j)+=FUNBWSE.row(countt);
  //     TFUNBWS.row(j)+=FUNBWS.row(countt);
  //     
  //     countt++;
  //   }
  //   
  // }
  return 0;
  // return Rcpp::List::create(Rcpp::Named("FUNENW")=FUNENW,
  //                           Rcpp::Named("FUNEBNW")=FUNEBNW,
  //                           Rcpp::Named("FUNEBSNW")=FUNEBSNW,
  //                           Rcpp::Named("FUNE")=FUNE,
  //                           Rcpp::Named("FUNBW")=FUNBW,
  //                           Rcpp::Named("FUNBWE")=FUNBWE,
  //                           Rcpp::Named("FUNBWSE")=FUNBWSE,
  //                           Rcpp::Named("FUNBWS")=FUNBWS);

}

