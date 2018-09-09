#include <RcppArmadillo.h>
#include <functions.hpp>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]] 
using namespace Rcpp;
using namespace arma;


arma::vec dfun(arma::vec par, arma::mat betam, int &n, arma::uvec &nJ, 
               Rcpp::List &matls, Rcpp::List &diagls, 
               arma::vec &indexy, arma::vec &uindexy, arma::vec &y, arma::mat &x)
{
  double sigmav2 = par(0);
  double sigmae2 = par(1);
  
  int n1 = 0;
  arma::mat betami;
  arma::mat M1;
  arma::mat M2;
  arma::uvec indexi;
  arma::vec d11  = zeros<vec>(1);
  arma::vec d12 = zeros<vec>(1);
  for(int i = 0; i < n; i++ )
  {
    betami = trans(betam.row(i));
    arma::mat mat1 =  matls[i];
    arma::mat diag1 = diagls[i];
    n1 = nJ(i);
    
    M1 = - mat1 * trans(mat1)*(1/pow(sigmae2 + n1*sigmav2,2));
    M2 = (sigmav2*(2*sigmae2 + n1*sigmav2)/pow(sigmae2+ n1*sigmav2,2)*mat1*trans(mat1) 
            - diag1)/pow(sigmae2,2);
    
    indexi = find(indexy == uindexy(i));
  
    d11 = d11 -0.5*n1/(sigmae2+ n1*sigmav2) - 0.5*trans(y(indexi) - x.rows(indexi)* betami)* M1* (y(indexi) - x.rows(indexi)*betami);
    d12 = d12 - 0.5*n1*(1- sigmav2/(sigmae2 +n1*sigmav2))/sigmae2 - 0.5*trans(y(indexi) - x.rows(indexi)*betami)*M2*(y(indexi) - x.rows(indexi)*betami);
  }
  
  arma::vec dvalue = zeros<vec>(2);
  dvalue(0) = d11(0);
  dvalue(1) = d12(0);
  
  return(dvalue);
}

arma::mat Ifun(arma::vec par, arma::uvec &nJ)
{
  double sigmav2 = par(0);
  double sigmae2 = par(1);
  
  arma::vec nJ1 = nJ + zeros<vec>(nJ.size());
  
  arma::mat Imat1 = zeros<mat>(2,2);
  double I11 = 0;
  double I12 = 1;
  double I22 = 1;


  I11 = 0.5*sum(square(nJ1)/square(sigmae2 + nJ1*sigmav2));
  I12 = 0.5*sum(nJ1/square(sigmae2 + nJ1*sigmav2));
  I22 = 0.5*sum(nJ1 % (1- sigmav2*(2*sigmae2 + nJ1*sigmav2)/square(sigmae2 + nJ1*sigmav2))/pow(sigmae2,2));
  
  Imat1(0,0) = I11;
  Imat1(0,1) = I12;
  Imat1(1,0) = I12;
  Imat1(1,1) = I22;
  return(Imat1);
}



// [[Rcpp::export]]
Rcpp::List Spgr_bhf2(arma::vec indexy,arma::vec &y,arma::mat &x, 
                    arma::vec &weights, arma::mat &betam0,
                    double nu = 1, double gam =3 , double lam = 0.5 ,
                    int maxiter = 100, double tolabs = 1e-4, double tolrel = 1e-2)
{
  int nt = x.n_rows;
  int p = x.n_cols;
  
  arma::vec uindexy = unique(indexy);
  int n = uindexy.size();
  
  int n0 = n*p;
  int npair = n*(n-1)/2;
  arma::mat Ip =  eye(p,p);
  
  arma::sp_mat Dmat = Dfun(n);
  
  arma::uvec nJ = zeros<uvec>(n);
  arma::uvec indexi;
  int ni;
  
  Rcpp::List matls(n);
  Rcpp::List diagls(n);
  arma::vec yhat = zeros(nt);
  
  for(int i = 0; i <n ; i++)
  {
    indexi = find(indexy == uindexy(i));
    ni = indexi.size();
    nJ(i) = ni;
    matls[i] = ones(ni,1);
    diagls[i] = eye(ni,ni);
    yhat(indexi) = x.rows(indexi)*trans(betam0.row(i));
  }
  
  double sige2 = mean(square(y - yhat));
  double siga2 = 0;
  
  arma::vec sig2est = zeros(2);
  sig2est(0) = siga2;
  sig2est(1) = sige2;
 

  //initial deltam
  arma::mat deltam(p,npair);
  arma::mat deltamold(p, npair);
  arma::mat betadiff(p,npair);
  
  arma::mat vm = zeros(p,npair);
  deltamold = trans(Dmat * betam0);
  
  
  // define some varaibles
  
  arma::vec temp = zeros<vec>(n0);
  arma::vec betanew = zeros<vec>(n0);
  arma::mat betam = zeros(n,p);
  arma::vec normbd(2);
  arma::vec Xty(n0);
  arma::mat XtX = zeros(n0,n0);
  arma::mat xm = zeros(nt,p);
  arma::mat Imat1 = zeros(2,2);
  arma::vec dvalue = zeros<vec>(2);

  int flag = 0;
  double rm  = 1;
  double sm = 1;
  double tolpri;
  double toldual;
  
  
  int m = 0;
  
  for(m = 0; m < maxiter; m++){
    
    for(int i= 0; i < n; i++)
    {
      ni = nJ(i);
      arma::mat mat1 = matls[i];
      arma::mat diag1 = diagls[i];
      arma::mat S1 =  mat1*trans(mat1)*siga2 + diag1*sige2;
      arma::mat S1inv = inv(S1);
      
      indexi = find(indexy == uindexy(i));
      
      xm.rows(indexi) = sqrtmat_sympd(S1inv) * x.rows(indexi);
      Xty(span(i*p,(i+1)*p - 1)) = trans(x.rows(indexi))*S1inv*y(indexi);
    }
    
    temp =  reshape((deltamold - 1/nu * vm)*Dmat,n0,1);
    betanew = inverserx(indexy,xm,nu)*(Xty + nu * temp);
    betam = trans(reshape(betanew, p, n));
    betadiff = trans(Dmat * betam);
    

    for(int j =0; j <10; j++)
    {
      Imat1 = Ifun(sig2est, nJ);
      dvalue = dfun(sig2est, betam, n, nJ, matls, diagls, indexy, uindexy, y, x);
      sig2est = sig2est + inv(Imat1)* dvalue;
      sig2est.elem(find(sig2est <0)).zeros();
    }

    //sig2est.elem(find(sig2est <0)).zeros();
    
    sige2 = sig2est(1);
    siga2 = sig2est(0); 
    
    deltam = betadiff + (1/nu) * vm;
    // update deltam
    for(int i = 0; i < npair; i++)
    {
      deltam.col(i) = scad(deltam.col(i),weights(i)*lam,nu,gam);
    }
    
    vm =  vm + nu * (betadiff - deltam);
    
    normbd(0) = norm(betadiff,"fro");
    normbd(1) = norm(deltam,"fro");
    
    tolpri = tolabs*sqrt(npair*p) + tolrel*max(normbd);
    toldual = tolabs*sqrt(n * p) + tolrel * norm(vm * Dmat, "fro");
    
    rm = norm(betadiff - deltam, "fro");
    sm = nu * norm((deltam - deltamold)*Dmat, "fro");
    
    deltamold = deltam;
    
    if(rm <= tolpri & sm <= toldual)
      break;
  }
  
  if(m == maxiter) {flag = 1;}
  
  arma::uvec group = getgroup(deltam,n,tolabs);
  arma::uvec ugroup = unique(group);
  int ng = ugroup.size();
  int nj;

  arma::uvec indexj;
  arma::mat betaest(ng,p);

  for(int j = 0; j < ng; j++)
  {
    indexj = find(group == ugroup(j));
    nj = indexj.size();
    betaest.row(j) = sum(betam.rows(indexj),0)/nj;
  }
  
  arma::mat betami;
  arma::vec loglikvalue = zeros(1);
  for(int i = 0; i < n; i++)
  {
    ni = nJ(i);
    arma::mat mat1 = matls[i];
    arma::mat diag1 = diagls[i];
    arma::mat S1 =  mat1*trans(mat1)*siga2 + diag1*sige2;
    arma::mat S1inv = inv(S1);

    indexi = find(indexy == uindexy(i));
    betami = trans(betaest.row(group(i) - 1));
    
    double val = 0;
    double sign = 0;

    log_det(val, sign, S1);

    loglikvalue = loglikvalue - 0.5*val -
     0.5*trans(y(indexi) - x.rows(indexi)*betami)*S1inv*(y(indexi) - x.rows(indexi)*betami);
  }
  
  
  return Rcpp::List::create(Named("beta") = betam,
                            Named("betaest") = betaest,
                            Named("index") = uindexy,
                            Named("sig2est") = sig2est,
                            Named("group") = group,
                            Named("loglikvalue") = loglikvalue,
                            Named("deltam") = deltam,
                            Named("flag") = flag,
                            Named("rm") = rm,
                            Named("sm") = sm,
                            Named("tolpri") = tolpri,
                            Named("toldual") = toldual,
                            Named("niteration") = m);
  
}


// [[Rcpp::export]]
double BICc_bhf2(Rcpp::List obj, double c0=0.2)
{
  arma::mat beta = obj("beta");
  arma::uvec group= obj("group");
  arma::uvec ugroup = unique(group);
  double loglikvalue = obj("loglikvalue");
  int ncx =  beta.n_cols;
  int nobs = beta.n_rows;
  double ngest = ugroup.size();
  double Cn = c0*log(log(nobs*ncx));
  
  double bicvalue = -2*loglikvalue + Cn*log(nobs)*(ngest*ncx);
  return(bicvalue);
}