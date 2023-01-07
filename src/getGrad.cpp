#include <RcppEigen.h>
#include "basics.h"

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
Eigen::MatrixXd getGradT(const Eigen::VectorXd & gamma, const Eigen::VectorXd & alpha, 
                         const Eigen::MatrixXd & H0,
                         const Eigen::MatrixXd & X2, const Eigen::VectorXd & survtime, 
                         const Eigen::VectorXd & status, const Eigen::VectorXd & ni, 
                         const Eigen::VectorXd & nt, const Eigen::VectorXd & PslT,
                         const Eigen::VectorXd & FUNE, 
                         const Eigen::MatrixXd & FUNBW,
                         const Eigen::MatrixXd & FUNBWE) {
  
  int d = gamma.size() + alpha.size();
  int a = H0.rows();
  int ks = X2.rows();
  int p1a = alpha.size();
  int i,p,q,j,t,u;
  
  double temp,temp1;
  
  Eigen::MatrixXd S = Eigen::MatrixXd::Zero(ks,d);
  
  int risk1_index;
  int risk1_index_temp=a-1;
  int risk1_index_ttemp=a-1;
  int risk1_index_tttemp=a-1;
  int risk1_index_vtemp=a-1;
  int risk1_index_vttemp=a-1;
  int risk1_index_vtttemp=a-1;
  
  
  Eigen::VectorXd CumuH0 = Eigen::VectorXd::Zero(a);
  
  temp1=0;
  for (j=0;j<a;j++) {
    temp1 += H0(j, 2);
    CumuH0(j) = temp1;
  }
  
  double epsilon=0;
  double qq=0;
  
  int p2 = gamma.size();
  
  Eigen::VectorXd X = Eigen::VectorXd::Zero(p2);
  Eigen::VectorXd SX = Eigen::VectorXd::Zero(p2);
  Eigen::VectorXd SX2 = Eigen::VectorXd::Zero(p2);
  Eigen::VectorXd SX1 = Eigen::VectorXd::Zero(p2);
  Eigen::VectorXd SRX = Eigen::VectorXd::Zero(p2);
  Eigen::VectorXd SRXX = Eigen::VectorXd::Zero(p2);
  Eigen::MatrixXd SXX1 = Eigen::MatrixXd::Zero(p2, a);
  Eigen::MatrixXd SXX11 = Eigen::MatrixXd::Zero(p2, a);
  Eigen::MatrixXd SRXX1 = Eigen::MatrixXd::Zero(p2, ks);
  Eigen::MatrixXd SRXX2 = Eigen::MatrixXd::Zero(p2, ks);
  
  Eigen::VectorXd N = Eigen::VectorXd::Zero(p1a);
  Eigen::VectorXd TN = Eigen::VectorXd::Zero(p1a);
  Eigen::VectorXd TN2 = Eigen::VectorXd::Zero(p1a);
  Eigen::VectorXd TN1 = Eigen::VectorXd::Zero(p1a);
  Eigen::VectorXd TRN = Eigen::VectorXd::Zero(p1a);
  Eigen::VectorXd TRNN = Eigen::VectorXd::Zero(p1a);
  Eigen::MatrixXd TRNN1 = Eigen::MatrixXd::Zero(p1a,ks);
  Eigen::MatrixXd TRNN2 = Eigen::MatrixXd::Zero(p1a,ks);
  Eigen::MatrixXd TNN1 = Eigen::MatrixXd::Zero(p1a,a);
  Eigen::MatrixXd TNN11 = Eigen::MatrixXd::Zero(p1a,a);
  
  
  for (j=0;j<ks;j++) {
    
    /* calculate score for gamma */
    if (j == 0)
    {
      temp=0;
      risk1_index=risk1_index_temp;
      for (q=j;q<ks;q++)
      {
        
        temp+=exp(MultVV(X2.row(q), gamma))*FUNE(q);
        SX+=exp(MultVV(X2.row(q), gamma))*FUNE(q)*X2.row(q);
        
        if (status(q) == 1)
        {
          if (q == ks-1)
          {
            if (temp != 0) {
              SX2 = H0(risk1_index, 1)/pow(temp, 2)*SX;
              SXX1.col(a-1-risk1_index) = SX2;
              SRXX += SX2;
              SX2 = 1/temp*SX;
              SXX11.col(a-1-risk1_index) = SX2;
            } else continue;
            risk1_index--;
          }
          else if (survtime(q+1) != survtime(q))
          {
            if (temp != 0) {
              SX2 = H0(risk1_index, 1)/pow(temp, 2)*SX;
              SXX1.col(a-1-risk1_index) = SX2;
              SRXX += SX2;
              SX2 = 1/temp*SX;
              SXX11.col(a-1-risk1_index) = SX2;
            } else continue;
            risk1_index--;
          }
          else
          {
            for (q=q+1;q<ks;q++)
            {
              temp+=exp(MultVV(X2.row(q), gamma))*FUNE(q);
              SX+=exp(MultVV(X2.row(q), gamma))*FUNE(q)*X2.row(q);
              
              if (q == ks-1)
              {
                if (temp != 0) {
                  SX2 = H0(risk1_index, 1)/pow(temp, 2)*SX;
                  SXX1.col(a-1-risk1_index) = SX2;
                  SRXX += SX2;
                  SX2 = 1/temp*SX;
                  SXX11.col(a-1-risk1_index) = SX2;
                } else continue;
                risk1_index--;
                break;
              }
              else if (survtime(q+1) != survtime(q))
              {
                if (temp != 0) {
                  SX2 = H0(risk1_index, 1)/pow(temp, 2)*SX;
                  SXX1.col(a-1-risk1_index) = SX2;
                  SRXX += SX2;
                  SX2 = 1/temp*SX;
                  SXX11.col(a-1-risk1_index) = SX2;
                } else continue;
                risk1_index--;
                break;
              }
              else continue;
            }
          }
          
        }
        else continue;
        
      }
      SRXX1.col(j) = SRXX;
    }
    else
    {
      if (risk1_index_temp>=0)
      {
        if (survtime(j) >= H0(risk1_index_temp, 0))
        {
          SRXX1.col(j) = SRXX1.col(j-1);
        }
        else
        {
          risk1_index_temp--;
          if (risk1_index_temp>=0)
          {
            SRXX = SRXX1.col(j-1);
            SRXX -= SXX1.col(a-1-risk1_index_temp-1);
            SRXX1.col(j) = SRXX;
          }
        }
      }
      else
      {
        risk1_index_temp=0;
      }
    }
    SRX = SRXX1.col(j);
    
    if (j==0)
    {
      SRX -= CumuH0(risk1_index_ttemp)*X2.row(j);  
    }
    else if (survtime(j) >= H0(risk1_index_ttemp, 0))
    {
      SRX -= CumuH0(risk1_index_ttemp)*X2.row(j);  
    }
    else
    {
      risk1_index_ttemp--;
      if (risk1_index_ttemp>=0)
      {
        SRX -= CumuH0(risk1_index_ttemp)*X2.row(j);  
      }
      else
      {
        SRX = SRXX1.col(j);
        risk1_index_ttemp=0;
      }
    }
    
    SRX*= exp(MultVV(X2.row(j), gamma))*FUNE(j);
    
    if (survtime(j) >= H0(risk1_index_tttemp, 0))
    {
      if (status(j) == 1)
      {
        X = PslT(j)*X2.row(j);
        X -= PslT(j)*SXX11.col(a-1-risk1_index_tttemp);
        X += SRX;
        for (q=0;q<p2;q++) S(j,q) = X(q);
      }
      else
      {
        for (q=0;q<p2;q++) S(j,q) = SRX(q);
      }
    }
    else
    {
      risk1_index_tttemp--;
      if (risk1_index_tttemp>=0)
      {
        if (status(j) == 1)
        {
          X = PslT(j)*X2.row(j);
          X -= PslT(j)*SXX11.col(a-1-risk1_index_tttemp);
          X += SRX;
          for (q=0;q<p2;q++) S(j,q) = X(q);
        }
        else
        {
          for (q=0;q<p2;q++) S(j,q) = SRX(q);
        }
      }
      else
      {
        risk1_index_tttemp=0;
        for (q=0;q<p2;q++) S(j,q) = SRX(q);
      }
    }
    
    /* calculate score for alpha */
    if (j == 0)
    {
      temp=0;
      
      TN = Eigen::VectorXd::Zero(p1a);
      TRN = Eigen::VectorXd::Zero(p1a);
      
      risk1_index=risk1_index_vtemp;
      for (q=j;q<ks;q++)
      {
        temp += exp(MultVV(X2.row(q), gamma))*FUNE(q);
        TN += exp(MultVV(X2.row(q), gamma))*FUNBWE.row(q);
        if (status(q) == 1)
        {
          if (q == ks-1)
          {
            if (temp != 0) {
              TN2=H0(risk1_index, 1)/pow(temp,2)*TN;
              TNN1.col(a-1-risk1_index) = TN2;
              TRNN += TN2;
              TN2= 1/temp*TN;
              TNN11.col(a-1-risk1_index) = TN2;
            } else continue;
            risk1_index--;
          }
          else if (survtime(q+1) != survtime(q))
          {
            if (temp != 0) {
              TN2=H0(risk1_index, 1)/pow(temp,2)*TN;
              TNN1.col(a-1-risk1_index) = TN2;
              TRNN += TN2;
              TN2= 1/temp*TN;
              TNN11.col(a-1-risk1_index) = TN2;
            } else continue;
            risk1_index--;
          }
          else
          {
            for (q=q+1;q<ks;q++)
            {
              temp += exp(MultVV(X2.row(q), gamma))*FUNE(q);
              TN += exp(MultVV(X2.row(q), gamma))*FUNBWE.row(q);
              if (q == ks-1)
              {
                if (temp != 0) {
                  TN2=H0(risk1_index, 1)/pow(temp,2)*TN;
                  TNN1.col(a-1-risk1_index) = TN2;
                  TRNN += TN2;
                  TN2= 1/temp*TN;
                  TNN11.col(a-1-risk1_index) = TN2;
                } else continue;
                risk1_index--;
                break;
              }
              else if (survtime(q+1) != survtime(q))
              {
                if (temp != 0) {
                  TN2=H0(risk1_index, 1)/pow(temp,2)*TN;
                  TNN1.col(a-1-risk1_index) = TN2;
                  TRNN += TN2;
                  TN2= 1/temp*TN;
                  TNN11.col(a-1-risk1_index) = TN2;
                } else continue;
                risk1_index--;
                break;
              }
              else continue;
            }
          }
          
        }
        else continue;
      }
      TRNN1.col(j) = TRNN;
      
    }
    else
    {
      if (risk1_index_vtemp>=0)
      {
        if (survtime(j) >= H0(risk1_index_vtemp, 0))
        {
          TRNN1.col(j) = TRNN1.col(j-1);
        }
        else
        {
          risk1_index_vtemp--;
          if (risk1_index_vtemp>=0)
          {
            TRNN = TRNN1.col(j-1);
            TRNN -= TNN1.col(a-1-risk1_index_vtemp-1);
            TRNN1.col(j) = TRNN;
          }
        }
      }
      else
      {
        risk1_index_vtemp=0;
      }
    }
    TRN = TRNN1.col(j);
    
    TRN *= exp(MultVV(X2.row(j), gamma))*FUNE(j);
    
    
    if (j==0) {
      N = FUNBWE.row(j);
      N *= CumuH0(risk1_index_vttemp)*exp(MultVV(X2.row(j), gamma));
      TRN -= N;
    } else if (survtime(j) >= H0(risk1_index_vttemp,0))
    {
      N = FUNBWE.row(j);
      N *= CumuH0(risk1_index_vttemp)*exp(MultVV(X2.row(j), gamma));
      TRN -= N;
    } else {
      risk1_index_vttemp--;
      if (risk1_index_vttemp>=0)
      {
        N = FUNBWE.row(j);
        N *= CumuH0(risk1_index_vttemp)*exp(MultVV(X2.row(j), gamma));
        TRN -= N;
      } else {
        risk1_index_vttemp=0;
      }
    }
    
    
    
    if (survtime(j) >= H0(risk1_index_vtttemp,0))
    {
      if (status(j) == 1)
      {
        for (q=0;q<p1a;q++) TN(q) = FUNBW(j,q) - PslT(j)*TNN11(q, a-1-risk1_index_vtttemp);
        TN += TRN;
        for (q=0;q<p1a;q++) S(j,p2+q) = TN(q);
      }
      else
      {
        for (q=0;q<p1a;q++) S(j,p2+q) = TRN(q);
      }
    } else
    {
      risk1_index_vtttemp--;
      if (risk1_index_vtttemp>=0)
      {
        if (status(j) == 1)
        {
          for (q=0;q<p1a;q++) TN(q) = FUNBW(j,q) - PslT(j)*TNN11(q, a-1-risk1_index_vtttemp);
          TN += TRN;
          for (q=0;q<p1a;q++) S(j,p2+q) = TN(q);
        }
        else
        {
          for (q=0;q<p1a;q++) S(j,p2+q) = TRN(q);
        }
      }
      else
      {
        risk1_index_vtttemp=0;
        for (q=0;q<p1a;q++) S(j,p2+q) = TRN(q);
      }
    }
    
    
  }
  
  return S;
  
}

//
// [[Rcpp::export]]
Eigen::MatrixXd getGradY(Eigen::VectorXd & beta, Eigen::VectorXd & tau, 
                         Eigen::MatrixXd & Sig, 
                         const Eigen::MatrixXd & Z, const Eigen::MatrixXd & X1, 
                         const Eigen::MatrixXd & W, const Eigen::VectorXd & Y, 
                         const Eigen::VectorXd & ni, const Eigen::VectorXd & nt, 
                         const Eigen::MatrixXd & Psl,
                         const Eigen::VectorXd & FUNENW, const Eigen::MatrixXd & FUNEBNW, 
                         const Eigen::MatrixXd & FUNEBSNW, const Eigen::MatrixXd & FUNBWS,
                         const double pStol) {
  
  int d = beta.size() + W.cols() + Sig.cols()*(Sig.cols() + 1)/2;
  
  int k = ni.size();
  int p1 = X1.cols();
  int p1a = Z.cols();
  int p1b = W.cols();
  
  int counti=0;
  int countt=0;
  
  int i,p,q,j,t,u;
  
  double temp,temp1;
  
  Eigen::MatrixXd S = Eigen::MatrixXd::Zero(k,d);
  
  
  double epsilon=0;
  double qq=0;
  
  Eigen::VectorXd SZ = Eigen::VectorXd::Zero(p1);
  Eigen::VectorXd SZ1 = Eigen::VectorXd::Zero(p1);
  Eigen::MatrixXd SZZ = Eigen::MatrixXd::Zero(p1, p1a);
  Eigen::VectorXd SZtau = Eigen::VectorXd::Zero(p1b);
  Eigen::MatrixXd bbT  = Eigen::MatrixXd::Zero(p1a, p1a);
  Eigen::MatrixXd bbT2  = Eigen::MatrixXd::Zero(p1a, p1a);
  
  Eigen::MatrixXd bsw = Eigen::MatrixXd::Zero(p1a+1,p1a+1);
  Eigen::MatrixXd bsw2 = Eigen::MatrixXd::Zero(p1a+1,p1a+1);
  
  for (j=0;j<k;j++) {
    /* calculate score for beta and tau*/
    for (t=0;t<nt(j);t++) {
      
      SZ = Eigen::VectorXd::Zero(p1);
      SZ1 = Eigen::VectorXd::Zero(p1);
      
      SZtau = Eigen::VectorXd::Zero(p1b);
      if (Psl(countt, 2) <= pStol) {
        counti += ni(j);
        countt++;
      } else {
        for(u=0;u<p1a;u++)   bbT(u,u) = FUNEBSNW(countt,u);
        if (p1a > 1) {
          for(i=1;i<p1a;i++)
          {
            for(u=0;u<p1a-i;u++) {
              bbT(u,i+u) = FUNEBSNW(countt, p1a+u+(i-1)*(p1a-1));
              bbT(i+u,u) = bbT(u,i+u);
            }
          }
        }
        for (i=0;i<ni(j);i++) {
          SZ1 += (Y(counti) - MultVV(X1.row(counti), beta))*
            exp(-MultVV(W.row(counti),tau))*FUNENW(countt)*X1.row(counti);
          SZ += exp(-MultVV(W.row(counti),tau))*
            MultVV(Z.row(counti), FUNEBNW.row(countt))*X1.row(counti);
          
          epsilon = Y(counti) - MultVV(X1.row(counti), beta);
          bbT2 = MultVVoutprod(Z.row(counti))*bbT;
          qq = pow(epsilon, 2)*FUNENW(countt) - 2*epsilon*MultVV(Z.row(counti), FUNEBNW.row(countt)) + bbT2.trace();
          SZtau += 0.5*(exp(-MultVV(W.row(counti),tau))*qq-Psl(countt, 2))*W.row(counti);
          counti++;
        }
        for (i=0;i<p1;i++) S(j,i) += SZ1(i) - SZ(i);
        for (i=0;i<p1b;i++) S(j,p1+i) += SZtau(i);
        
        /*calculate score for Sig*/
        for(u=0;u<(p1a+1);u++) bsw(u,u) += FUNBWS(countt,u);
        for(i=1;i<(p1a+1);i++)
        {
          for(u=0;u<(p1a+1-i);u++) {
            bsw(u,i+u) += FUNBWS(countt,p1a+1+u+(i-1)*p1a);
            bsw(i+u,u) = bsw(u,i+u);
          }
        }
        
        countt++; 
      }
      
    }
    
    bsw2 = Sig.inverse()*bsw*Sig.inverse() - Sig.inverse();
    
    for (u=0;u<(p1a+1);u++) S(j, p1+p1b+u) = 0.5*bsw2(u,u);
    
    for(q=1;q<(p1a+1);q++)
    {
      for(u=0;u<(p1a+1-q);u++) S(j, p1+p1b+p1a+1+u+(q-1)*(p1a)) = bsw2(u,q+u);
    }
    
    bsw = Eigen::MatrixXd::Zero(p1a+1,p1a+1);
    
  }
  
  
  
  return S;
  
}