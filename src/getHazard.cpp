#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

//
// [[Rcpp::export]]
int getHazard(const Eigen::Map<Eigen::VectorXd> & CumuH0,
              const Eigen::Map<Eigen::VectorXd> & survtime,
              const Eigen::Map<Eigen::VectorXd> & status,
              const Eigen::Map<Eigen::MatrixXd> & H0,
              Eigen::Map<Eigen::VectorXd> & CUH0,
              Eigen::Map<Eigen::VectorXd> & HAZ0) {
  
  int a=H0.rows();
  int k=survtime.size();
  
  int risk1_index=a-1;
  
  int j;
  
  for (j=0;j<k;j++)
  {
    if (risk1_index>=0)
    {
      if (survtime(j) >= H0(risk1_index, 0))
      {
        CUH0(j) = CumuH0(risk1_index);
        
      }
      else
      {
        risk1_index--;
        if (risk1_index>=0)
        {
          CUH0(j) = CumuH0(risk1_index);
        }
      }
    }
    else
    {
      risk1_index=0;
    }
  }

  risk1_index = a-1;
  
  for (j=0;j<k;j++)
  {
    if (risk1_index>=0)
    {
      //gsl_vector_set(CUH0,j,gsl_vector_get(CumuH0,risk1_index));
      if (survtime(j) == H0(risk1_index, 0))
      {
        HAZ0(j) = H0(risk1_index, 2);
      }
      if (status(j) == 1)
      {
        if (j == k-1)
        {
          risk1_index--;
        }
        else if (survtime(j+1) != survtime(j))
        {
          risk1_index--;
        }
        else
        {
          for (j=j+1;j<k;j++)
          {
            if (survtime(j) == H0(risk1_index, 0))
            {
              HAZ0(j) = H0(risk1_index, 2);
            }
            if (j == k-1)
            {
              risk1_index--;
              break;
            }
            else if (survtime(j+1) != survtime(j))
            {
              risk1_index--;
              break;
            }
            else continue;
          }
        }
      }
      else continue;
    }
    else continue;
  }
  
  return 0;
}