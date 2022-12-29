#include <RcppEigen.h>
#include "basics.h"

// [[Rcpp::depends(RcppEigen)]]

//
// [[Rcpp::export]]
Rcpp::List getMC(Eigen::VectorXd & beta, Eigen::VectorXd & tau, 
                 Eigen::VectorXd & gamma, Eigen::VectorXd & alpha,
                 Eigen::MatrixXd & H0, Eigen::MatrixXd & Sig, 
                 const Eigen::MatrixXd & Z, const Eigen::MatrixXd & X1, 
                 const Eigen::MatrixXd & W, const Eigen::VectorXd & Y, 
                 const Eigen::MatrixXd & X2, const Eigen::VectorXd & survtime, 
                 const Eigen::VectorXd & status, const Eigen::VectorXd & ni,
                 const Eigen::VectorXd & nt, const Eigen::MatrixXd & Psl,
                 const Eigen::VectorXd & PslT,
                 const Eigen::VectorXd & FUNENW, const Eigen::MatrixXd & FUNEBNW, 
                 const Eigen::MatrixXd & FUNEBSNW, const Eigen::VectorXd & FUNE, 
                 const Eigen::MatrixXd & FUNBW, const Eigen::MatrixXd & FUNBWE, const Eigen::MatrixXd & FUNBWSE, 
                 const Eigen::MatrixXd & FUNBWS){
  
  int a =H0.rows();
  int k = ni.size();
  int ks = X2.rows();
  int p1 = X1.cols();
  int p1a = Z.cols();
  int p1b = W.cols();
  
  int counti=0;
  int countt=0;
  
  int i,p,q,j,t,u;
  
  double scalef=0;
  
  Eigen::MatrixXd X1X1Ts  = Eigen::MatrixXd::Zero(p1, p1);
  Eigen::MatrixXd X1X1T  = Eigen::MatrixXd::Zero(p1, p1);
  Eigen::VectorXd X11s = Eigen::VectorXd::Zero(p1);
  Eigen::VectorXd X11 = Eigen::VectorXd::Zero(p1);
  Eigen::MatrixXd WWT  = Eigen::MatrixXd::Zero(p1b, p1b);
  Eigen::MatrixXd WWTs  = Eigen::MatrixXd::Zero(p1b, p1b);
  Eigen::VectorXd W11 = Eigen::VectorXd::Zero(p1b);
  Eigen::VectorXd W11s = Eigen::VectorXd::Zero(p1b);
  Eigen::MatrixXd bbT  = Eigen::MatrixXd::Zero(p1a, p1a);
  Eigen::MatrixXd bbT2  = Eigen::MatrixXd::Zero(p1a, p1a);
  
  /* calculate beta */
  for (j=0;j<k;j++) {
    for (t=0;t<nt(j);t++) {
      for (i=0;i<ni(j);i++)
      {
        X1X1Ts += MultVVoutprod(X1.row(counti))/exp(MultVV(W.row(counti),tau));
        X11s += (Y(counti)*FUNENW(j)*X1.row(counti)-
          MultVV(Z.row(counti), FUNEBNW.row(j))*X1.row(counti))/exp(MultVV(W.row(counti),tau)); 
        counti++;
      }
      X1X1Ts *= Psl(countt, 2)*FUNENW(j);
      X11s *= Psl(countt, 2);
      X1X1T += X1X1Ts;
      X11 += X11s;
      X1X1Ts = Eigen::MatrixXd::Zero(p1, p1);
      X11s = Eigen::VectorXd::Zero(p1);
      countt++; 
    }
    
  }
  beta = X1X1T.inverse()*X11;
  
  double epsilon=0;
  double qq=0;

  Eigen::MatrixXd bbT3  = Eigen::MatrixXd::Zero(p1b, p1b);
  Eigen::MatrixXd bbT4  = Eigen::MatrixXd::Zero(p1b, p1b);
  Eigen::VectorXd bT3 = Eigen::VectorXd::Zero(p1b);
  int aa,ss;
  /* calculate tau */
  counti = 0;
  countt = 0;
  for (j=0;j<k;j++) {
    
    for(t=0;t<p1a;t++)   bbT(t,t) = FUNEBSNW(j,t);
    if (p1a > 1) {
      for(i=1;i<p1a;i++)
      {
        for(t=0;t<p1a-i;t++) {
          bbT(t,i+t) = FUNEBSNW(j, p1a+t+(i-1)*(p1a-1));
          bbT(i+t,t) = bbT(t,i+t);
        }
      }
    }

    for (t=0;t<nt(j);t++) {
      for (i=0;i<ni(j);i++)
      {
        epsilon = Y(counti) - MultVV(X1.row(counti), beta);
        bbT2 = MultVVoutprod(Z.row(counti))*bbT;
        qq = pow(epsilon, 2)*FUNENW(j) - 2*epsilon*MultVV(Z.row(counti), FUNEBNW.row(j)) 
          + bbT2.trace();
        WWTs += qq*0.5*exp(-MultVV(W.row(counti),tau))*MultVVoutprod(W.row(counti));
        W11s += 0.5*(exp(-MultVV(W.row(counti),tau))*qq-1)*W.row(counti); 
        counti++;
      }
      WWTs *= Psl(countt, 2);
      W11s *= Psl(countt, 2);
      WWT += WWTs;
      W11 += W11s;
      WWTs = Eigen::MatrixXd::Zero(p1b, p1b);
      W11s = Eigen::VectorXd::Zero(p1b);
      countt++; 
    }
    
  }

  tau+=WWT.inverse()*W11;
  
  /* calculate Sig*/
  Sig = Eigen::MatrixXd::Zero(p1a+1, p1a+1);
  
  for (j=0;j<k;j++) {
    
    for(t=0;t<(p1a+1);t++) Sig(t,t)+=FUNBWS(j,t);
    for(q=1;q<(p1a+1);q++)
    {
      for(t=0;t<(p1a+1-q);t++) {
        Sig(t,q+t)+=FUNBWS(j, p1a+1+t+(q-1)*(p1a));
        Sig(q+t,t)+=FUNBWS(j, p1a+1+t+(q-1)*(p1a));
      }
    }
  }
  Sig/=k;
  
  /* calculate H0*/
  double dem=0;
  int risk1_index=a-1;

  for (j=0;j<ks;j++)
  {

    dem+=FUNE(j)*exp(MultVV(X2.row(j), gamma))*PslT(j);

    if (status(j) == 1)
    {

      if (j == ks-1)
      {
        if (dem != 0) {
          H0(risk1_index, 2)=H0(risk1_index, 1)/dem;
        } else {
          H0(risk1_index, 2)=0;
        }
        risk1_index--;
      }
      else if (survtime(j+1) != survtime(j))
      {
        if (dem != 0) {
          H0(risk1_index, 2)=H0(risk1_index, 1)/dem;
        } else {
          H0(risk1_index, 2)=0;
        }
        risk1_index--;
      }

      else
      {
        for (j=j+1;j<ks;j++)
        {

          dem+=FUNE(j)*exp(MultVV(X2.row(j), gamma))*PslT(j);

          if (j == ks-1)
          {
            if (dem != 0) {
              H0(risk1_index, 2)=H0(risk1_index, 1)/dem;
            } else {
              H0(risk1_index, 2)=0;
            }
            risk1_index--;
            break;
          }
          else if (survtime(j+1) != survtime(j))
          {
            if (dem != 0) {
              H0(risk1_index, 2)=H0(risk1_index, 1)/dem;
            } else {
              H0(risk1_index, 2)=0;
            }
            risk1_index--;
            break;
          }
          else continue;
        }
      }

    }
    else continue;
  }


  /* calculate gamma */
  double scalefH0=0;

  risk1_index=a-1;

  int p2=gamma.size();

  Eigen::VectorXd SX = Eigen::VectorXd::Zero(p2);
  Eigen::VectorXd SX2 = Eigen::VectorXd::Zero(p2);
  Eigen::VectorXd X22 = Eigen::VectorXd::Zero(p2);
  Eigen::MatrixXd XX = Eigen::MatrixXd::Zero(p2, p2);
  Eigen::MatrixXd SXX = Eigen::MatrixXd::Zero(p2, p2);
  Eigen::MatrixXd SXX2 = Eigen::MatrixXd::Zero(p2, p2);
  Eigen::VectorXd SX_new = Eigen::VectorXd::Zero(p2);
  Eigen::VectorXd SX_inter = Eigen::VectorXd::Zero(p2);
  Eigen::MatrixXd SXX_new = Eigen::MatrixXd::Zero(p2, p2);

  for (j=0;j<ks;j++)
  {
    XX = MultVVoutprod(X2.row(j));
    scalef = FUNE(j)*exp(MultVV(X2.row(j), gamma))*PslT(j);
    XX*=scalef;
    SXX+=XX;
    X22=X2.row(j);
    X22*=scalef;
    SX+=X22;

    if (status(j) == 1)
    {
      if (j == ks-1)
      {
        scalefH0 = H0(risk1_index, 2);
        SXX2=SXX*scalefH0;
        SXX_new+=SXX2;
        SX2=SX*scalefH0;
        SX_new+=SX2;
        risk1_index--; 
      }
      else if (survtime(j+1) != survtime(j))
      {
        scalefH0 = H0(risk1_index, 2);
        SXX2=SXX*scalefH0;
        SXX_new+=SXX2;
        SX2=SX*scalefH0;
        SX_new+=SX2;
        risk1_index--;
      }
      else
      {
        for (j=j+1;j<ks;j++)
        {
          XX = MultVVoutprod(X2.row(j));
          scalef = FUNE(j)*exp(MultVV(X2.row(j), gamma))*PslT(j);
          XX*=scalef;
          SXX+=XX;
          X22=X2.row(j);
          X22*=scalef;
          SX+=X22;

          if (j == ks-1)
          {
            scalefH0 = H0(risk1_index, 2);
            SXX2=SXX*scalefH0;
            SXX_new+=SXX2;
            SX2=SX*scalefH0;
            SX_new+=SX2;
            risk1_index--;
            break;
          }
          else if (survtime(j+1) != survtime(j))
          {
            scalefH0 = H0(risk1_index, 2);
            SXX2=SXX*scalefH0;
            SXX_new+=SXX2;
            SX2=SX*scalefH0;
            SX_new+=SX2;
            risk1_index--;
            break;
          }
          else continue;
        }
      }
    }
    else continue;
  }

  for (j=0;j<ks;j++)
  {
    if (status(j) == 1) SX_inter+=PslT(j)*X2.row(j);
  }
  
  gamma+=SXX_new.inverse()*(SX_inter - SX_new);

  /*  calculate alpha*/

  Eigen::MatrixXd bwT = Eigen::MatrixXd::Zero(p1a+1, p1a+1);
  Eigen::MatrixXd TD = Eigen::MatrixXd::Zero(p1a+1, p1a+1);
  Eigen::MatrixXd TD2 = Eigen::MatrixXd::Zero(p1a+1, p1a+1);
  Eigen::MatrixXd TDD = Eigen::MatrixXd::Zero(p1a+1, p1a+1);
  Eigen::VectorXd TN = Eigen::VectorXd::Zero(p1a+1);
  Eigen::VectorXd TN2 = Eigen::VectorXd::Zero(p1a+1);
  Eigen::VectorXd TNN = Eigen::VectorXd::Zero(p1a+1);
  Eigen::VectorXd N = Eigen::VectorXd::Zero(p1a+1);

  risk1_index = a-1;

  for (j=0;j<ks;j++)
  {

    for(t=0;t<(p1a+1);t++) bwT(t,t) = FUNBWSE(j,t);

    for(i=1;i<(p1a+1);i++)
    {
      for(t=0;t<(p1a+1-i);t++) {
        bwT(t,i+t) = FUNBWSE(j,p1a+1+t+(i-1)*(p1a));
        bwT(i+t,t) = bwT(t,i+t);
      }
    }

    N = FUNBWE.row(j);

    bwT*=exp(MultVV(X2.row(j), gamma))*PslT(j);
    TD+=bwT;
    N*=exp(MultVV(X2.row(j), gamma))*PslT(j);
    TN+=N;

    if (status(j) == 1)
    {
      if (j == ks-1)
      {
        TD2=TD*H0(risk1_index,2);
        TDD+=TD2;
        TN2=TN*H0(risk1_index,2);
        TNN+=TN2;
        risk1_index--;

      }
      else if (survtime(j+1) != survtime(j))
      {
        TD2=TD*H0(risk1_index,2);
        TDD+=TD2;
        TN2=TN*H0(risk1_index,2);
        TNN+=TN2;
        risk1_index--;
      }
      else
      {
        for (j=j+1;j<ks;j++)
        {
          
          for(t=0;t<(p1a+1);t++) bwT(t,t) = FUNBWSE(j,t);
          
          for(i=1;i<(p1a+1);i++)
          {
            for(t=0;t<(p1a+1-i);t++) {
              bwT(t,i+t) = FUNBWSE(j,p1a+1+t+(i-1)*(p1a));
              bwT(i+t,t) = bwT(t,i+t);
            }
          }
          
          N = FUNBWE.row(j);
          
          bwT*=exp(MultVV(X2.row(j), gamma))*PslT(j);
          TD+=bwT;
          N*=exp(MultVV(X2.row(j), gamma))*PslT(j);
          TN+=N;

          if (j == ks-1)
          {
            TD2=TD*H0(risk1_index,2);
            TDD+=TD2;
            TN2=TN*H0(risk1_index,2);
            TNN+=TN2;
            risk1_index--;
            break;
          }
          else if (survtime(j+1) != survtime(j))
          {
            TD2=TD*H0(risk1_index,2);
            TDD+=TD2;
            TN2=TN*H0(risk1_index,2);
            TNN+=TN2;
            risk1_index--;
            break;
          }
          else continue;
        }
      }
    }
    else continue;
  }

  TN = Eigen::VectorXd::Zero(p1a+1);

  for (j=0;j<ks;j++)
  {
    if(status(j)==1)
    {
      N = PslT(j)*FUNBW.row(j);
      TN+=N;
    } else continue;
  }

  N=TDD.inverse()*(TN-TNN);
  alpha+=N;
  
 
  
  return Rcpp::List::create(Rcpp::Named("beta")=beta,
                            Rcpp::Named("tau")=tau,
                            Rcpp::Named("Sig")=Sig,
                            Rcpp::Named("H0")=H0,
                            Rcpp::Named("gamma")=gamma,
                            Rcpp::Named("alpha")=alpha);
  
  
}