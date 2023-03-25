#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

//
// [[Rcpp::export]]
Rcpp::List GetrisksetC(const Eigen::MatrixXd & cdata) {
  
  int k = cdata.rows();
  int a=0;
  double u=0;
  int i,j;
  Eigen::MatrixXd FH01 = Eigen::MatrixXd::Zero(k, 2);
  
  /* find # events for risk 1*/
  for (j=0;j<k;j++)
  {
    if (cdata(j,1) == 1)
    {
      u = u + cdata(j,2);
      if (j == k-1)
      {
        a++;
        FH01(k-a,0) = cdata(j,3);
        FH01(k-a,1) = u;
        u=0;
      }
      else if (cdata(j+1,3) != cdata(j,3))
      {
        a++;
        FH01(k-a,0) = cdata(j,3);
        FH01(k-a,1) = u;
        u=0;
      }
      else
      {
        for (j=j+1;j<k;j++)
        {
          if (cdata(j,1) == 1)
          {
            u = u + cdata(j,2);
            if (j == k-1)
            {
              a++;
              FH01(k-a,0) = cdata(j,3);
              FH01(k-a,1) = u;
              u=0;
              break;
            }
            else if (cdata(j+1,3) != cdata(j,3))
            {
              a++;
              FH01(k-a,0) = cdata(j,3);
              FH01(k-a,1) = u;
              u=0;
              break;
            }
            else continue;
          }
          else
          {
            if (j == k-1)
            {
              a++;
              FH01(k-a,0) = cdata(j,3);
              FH01(k-a,1) = u;
              u=0;
              break;
            }
            else if (cdata(j+1,3) != cdata(j,3))
            {
              a++;
              FH01(k-a,0) = cdata(j,3);
              FH01(k-a,1) = u;
              u=0;
              break;
            }
            else continue;
          }
        }
      }
      
    }
    else continue;
  }
  
  if(a==0)
  {
    printf("No failure time information for risk 1; Program exits\n");
    return ( -1.0 );
  } 
  
  Eigen::MatrixXd H01 = Eigen::MatrixXd::Zero(a, 3);
  for(i=0;i<3;i++)
  {
    if(i<=1)
    {
      for(j=a;j>0;j--)    H01(a-j,i) = FH01(k-j,i);
    }
    if(i==2)
    {
      for(j=0;j<a;j++)    H01(j,i) = 0.000001;
    }
  }
  
  return Rcpp::List::create(Rcpp::Named("H01")=H01);
  
}