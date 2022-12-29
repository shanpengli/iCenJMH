#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

//
// [[Rcpp::export]]
Eigen::MatrixXd GetrisksetS(const Eigen::MatrixXd & iCendata) {
  
  int k = iCendata.rows();
  int a=0;
  double u=0;
  int i,j;
  Eigen::MatrixXd FH01 = Eigen::MatrixXd::Zero(k, 2);
  
  /* find # events for risk 1*/
  for (j=0;j<k;j++)
  {
    u += iCendata(j, 2);
    if (j == k-1)
    {
      a++;
      FH01(k-a,0) = iCendata(j,1);
      FH01(k-a,1) = u;
      u=0;
    }
    else if (iCendata(j+1,1) != iCendata(j,1))
    {
      a++;
      FH01(k-a,0) = iCendata(j,1);
      FH01(k-a,1) = u;
      u=0;
    }
    else 
    {
      for (j=j+1;j<k;j++)
      {
        u = u + iCendata(j,2);
        if (j == k-1)
        {
          a++;
          FH01(k-a,0) = iCendata(j,1);
          FH01(k-a,1) = u;
          u=0;
          break;
        }
        else if (iCendata(j+1,1) != iCendata(j,1))
        {
          a++;
          FH01(k-a,0) = iCendata(j,1);
          FH01(k-a,1) = u;
          u=0;
          break;
        }
        else continue;
      }
    }
}
  
  Eigen::MatrixXd phiS = Eigen::MatrixXd::Zero(a, 3);
  
  for(i=0;i<3;i++)
  {
    if(i<=1)
    {
      for(j=a;j>0;j--)    phiS(a-j,i) = FH01(k-j,i);
    }
  }
  double totalP = phiS.col(1).sum();
  for(j=0;j<a;j++)
  {
    u = 1 - phiS(j,1)/totalP;
    if (u <= 0) {
      phiS(j,2) = 30;
    } else {
      phiS(j,2) = -log(u);
    }
    totalP -= phiS(j,1);
  }
  phiS(a-1,2) = 30;
  return phiS;
  
}