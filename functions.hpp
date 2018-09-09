#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <RcppArmadillo.h>
#define ARMA_64BIT_WORD 1
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]] 
using namespace Rcpp;
using namespace arma;


// sfun
arma::vec sfun(arma::vec x, double th)
{
  double xn = norm(x,2);
  double thval = 1 - th/xn;
  x = thval*(thval >0)*x;
  return(x);
}


//scad 
// [[Rcpp::export]]
arma::vec scad(arma::vec v, double lam,  double nu , 
               double gam ) 
{
  
  double temp1 = lam/nu;
  double temp2 = gam * lam;
  double xn = norm(v,2);
  arma::vec vnew;
  
  if(xn <= lam + temp1)
  {
    vnew = sfun(v, temp1);
  }else if(xn <= temp2 & xn >= lam + temp1)
  {
    vnew = sfun(v, temp2/((gam-1)*nu))/(1- 1/((gam - 1 )*nu));
  }else{
    vnew = v;
  }
  return(vnew);
}

//inverse 
arma::mat inverse( arma::mat &x, arma::mat &z, 
                  double nu)
{
  int n = x.n_rows;
  int p = x.n_cols;
  int q = z.n_cols;
  
  int n0 = n*p;
  arma::mat Ip = 1/nu * eye(p,p);
  arma::mat nIp = nu* n * eye(p,p);
  arma::mat xt = x.t();
  arma::mat zt = z.t();
  arma::mat matinv = zeros(n0,n0); 
  
  
  arma::mat DB = zeros(p,p);
  arma::mat zAz = zeros(q, q);
  arma::mat Az = zeros(p,q);
  arma::mat AB = zeros(n0, p);
  arma::mat AZ = zeros(n0,q);
  
  arma::mat mati = zeros(p,p);
  
  int idp1 = 0;
  int idp2 = p - 1;
  for(int i=0; i < n; i++){
    
    mati = inv(xt.col(i) * x.row(i) + nIp);
    DB = DB + mati;
    zAz = zAz + zt.col(i)*z.row(i) -  zt.col(i)*x.row(i) *mati * xt.col(i)*z.row(i);
    Az = Az + mati * xt.col(i)*z.row(i);
    AB.rows(idp1,idp2) = mati;
    AZ.rows(idp1,idp2) = mati*xt.col(i)* z.row(i);
    matinv.submat(idp1,idp1,idp2,idp2) = mati; // output
    
    idp1 = idp1 + p;
    idp2 = idp2 + p;
  }
  
  arma::mat zAzinv = zAz.i();
  arma::mat Azt = Az.t();
  arma::mat temp1 = inv(Ip - DB - Az * zAzinv* Azt);
  
  arma::mat A1 = AB * temp1 * AB.t();
  arma::mat A2 = AZ * zAzinv * Az.t();
  arma::mat A3 = A2 * temp1;
  arma::mat A22 = A3 * A2.t();
  arma::mat A12 = A3 * AB.t();
  
  matinv = matinv + AZ * zAzinv * AZ.t() + A1 + A22 + A12 + A12.t();
  return(matinv);
  
}

arma::mat inversex( arma::mat &x, double nu)
{
  int n = x.n_rows;
  int p = x.n_cols;
  
  int n0 = n*p;
  arma::mat Ip = 1/nu * eye(p,p);
  arma::mat nIp = nu* n * eye(p,p);
  arma::mat xt = x.t();
  arma::mat matinv = zeros(n0,n0); 
  
  
  arma::mat DB = zeros(p,p);
  arma::mat AB = zeros(n0, p);
  arma::mat mati = zeros(p,p);
  
  int idp1 = 0;
  int idp2 = p - 1;
  for(int i=0; i < n; i++){
    
    mati = inv(xt.col(i) * x.row(i) + nIp);
    DB = DB + mati;
    AB.rows(idp1,idp2) = mati;
    matinv.submat(idp1,idp1,idp2,idp2) = mati; // output
    
    idp1 = idp1 + p;
    idp2 = idp2 + p;
  }
  
  arma::mat IB = inv(Ip - DB);
  
  matinv = matinv + AB * IB * AB.t();
  return(matinv);
  
}




// sparse matrix
arma::sp_mat Dfun(int n)
{
  int npair = n*(n-1)/2;
  int index1;
  arma::mat loc = zeros(2,npair*2);
  
  
  for(int j = 0; j < n-1; j++)
  {
    index1 = (2*n -j- 1)*j/2;
    loc(0,span(index1,index1 + n - j - 2)) = linspace<rowvec>(index1,index1+n-j-2,n-j-1);
    loc(1,span(index1,index1 + n - j - 2)) = j*ones<rowvec>(n-j-1);
    
    loc(0,span(index1 + npair,index1 + n - j - 2 + npair)) = linspace<rowvec>(index1,index1+n-j-2,n-j-1);
    loc(1,span(index1 + npair,index1 + n - j - 2 + npair)) = linspace<rowvec>(j+1,n -1,n-j-1);
  }
  
  arma::umat loc1 =  arma::conv_to<arma::umat>::from(loc);
  
  arma::vec values = zeros<vec>(npair*2);
  values.head(npair) = ones(npair);
  values.tail(npair) = (-1)*ones(npair);
  
  arma::sp_mat res(loc1,values,npair,n);
  
  return(res);
}

// a function to implement in for uvec index of a%in%b
arma::uvec findin(arma::uvec &a, arma::uvec &b) 
{
  int asize = a.size();
  arma::uvec indexa = zeros<uvec>(asize);
  arma::uvec bi;
  
  for(int i=0; i<asize; i++)
  {
    bi = find( b == a(i));
    if(bi.size() >0){indexa(i) = 1;}
  }
  return(indexa);
}

arma::uvec getgroup(arma::mat &b2pair, 
                    int n, double tol = 1e-4)
{
  int p = b2pair.n_rows;
  arma::vec b2value  =  arma::conv_to<arma::vec>::from(sqrt(sum(square(b2pair),0)/p));
  arma::uvec indexb = find(abs(b2value) <= tol);
  b2value(indexb) = zeros(size(indexb));
  
  arma::mat d2 = zeros(n,n);
  int indexj1;
  int indexj2;
  
  for(int j = 0; j < n - 1; j++)
  {
    indexj1 = (2*n -j- 1)*j/2;
    indexj2 = indexj1 + n - j - 2;
    d2(span(n - indexj2 + indexj1 - 1, n - 1),j) = b2value(span(indexj1,indexj2));
  }
  
  d2 = d2.t() + d2;
  
  int j = 0;
  arma::uvec ngj = linspace<uvec>(0, n - 1, n);
  arma::uvec gj;
  arma::uvec groupest = zeros<uvec>(n);
  
  do{
    j = j + 1;
    gj = find(d2.row(ngj(0))==0);
    gj = ngj(find(findin(ngj, gj)==1));
    ngj = ngj(find(findin(ngj,gj)==0));
    groupest(gj) = j*ones<uvec>(gj.size());
  } while (ngj.size()>0);
  
  return(groupest);
}

double BIClog(List &obj, arma::vec &y,  arma::mat &z,  arma::mat &x, 
              arma::uvec &group, double tol )
{
  int n = x.n_rows;
  int p = x.n_cols;
  int q = z.n_cols;
  arma::uvec ugroup = unique(group);
  int ng = ugroup.size();
  int nj;
  
  arma::mat beta = obj("beta");
  arma::vec eta = obj("eta");
  arma::vec estexp(n);
  arma::uvec indexj;
  arma::rowvec betaj(p);
  
  for(int j = 0; j < ng; j++)
  {
    indexj = find(group == ugroup(j));
    nj = indexj.size();
    betaj = sum(beta.rows(indexj),0)/nj;
    estexp(indexj) = x.rows(indexj) * betaj.t();
  }
  
  estexp = estexp + z* eta;
  
  double Cn = log(n * p + q);
  
  double bicvalue =  log(sum(square(y - estexp))/n) + Cn * log(n) *(ng * p + q)/n;
  
  return(bicvalue);
}
double BIClogx(List &obj, arma::vec &y,arma::mat &x,
               arma::uvec &group, double tol = 1e-4)
{
  int n = x.n_rows;
  int p = x.n_cols;
  arma::uvec ugroup = unique(group);
  int ng = ugroup.size();
  int nj;
  
  arma::mat beta = obj("beta");
  arma::vec estexp(n);
  arma::uvec indexj;
  arma::rowvec betaj(p);
  
  for(int j = 0; j < ng; j++)
  {
    indexj = find(group == ugroup(j));
    nj = indexj.size();
    betaj = sum(beta.rows(indexj),0)/nj;
    estexp(indexj) = x.rows(indexj) * betaj.t();
  }
  
  
  double Cn = log(n * p);
  
  double bicvalue =  log(sum(square(y - estexp))/n) + Cn * log(n) *(ng * p)/n;
  
  return(bicvalue);
}


arma::mat inverser( arma::vec indexy, arma::mat &x, arma::mat &z, 
                    double nu)
{
  int p = x.n_cols;
  int q = z.n_cols;
  
  arma::vec uindexy = unique(indexy);
  int n = uindexy.size();
  
  int n0 = n*p;
  arma::mat Ip = 1/nu * eye(p,p);
  arma::mat nIp = nu* n * eye(p,p);
  arma::mat xt = x.t();
  arma::mat zt = z.t();
  arma::mat matinv = zeros(n0,n0);
  
  
  arma::mat DB = zeros(p,p);
  arma::mat zAz = zeros(q, q);
  arma::mat Az = zeros(p,q);
  arma::mat AB = zeros(n0, p);
  arma::mat AZ = zeros(n0,q);
  
  arma::mat mati = zeros(p,p);
  
  arma::uvec indexi;
  int idp1 = 0;
  int idp2 = p - 1;
  for(int i=0; i < n; i++){
    
    indexi = find(indexy == uindexy(i));
    
    mati = inv(xt.cols(indexi) * x.rows(indexi) + nIp);
    DB = DB + mati;
    zAz = zAz + zt.cols(indexi)*z.rows(indexi) -  
      zt.cols(indexi)*x.rows(indexi) *mati * xt.cols(indexi)*z.rows(indexi);
    
    Az = Az + mati * xt.cols(indexi)*z.rows(indexi);
    AB.rows(idp1,idp2) = mati;
    AZ.rows(idp1,idp2) = mati*xt.cols(indexi)* z.rows(indexi);
    matinv.submat(idp1,idp1,idp2,idp2) = mati; // output
    
    idp1 = idp1 + p;
    idp2 = idp2 + p;
  }
  
  arma::mat zAzinv = zAz.i();
  arma::mat Azt = Az.t();
  arma::mat temp1 = inv(Ip - DB - Az * zAzinv* Azt);
  
  arma::mat A1 = AB * temp1 * AB.t();
  arma::mat A2 = AZ * zAzinv * Az.t();
  arma::mat A3 = A2 * temp1;
  arma::mat A22 = A3 * A2.t();
  arma::mat A12 = A3 * AB.t();
  
  matinv = matinv + AZ * zAzinv * AZ.t() + A1 + A22 + A12 + A12.t();
  return(matinv);
  
}

arma::mat inverserx( arma::vec indexy, arma::mat &x,double nu)
{
  int p = x.n_cols;
  
  arma::vec uindexy = unique(indexy);
  int n = uindexy.size();
  
  int n0 = n*p;
  arma::mat Ip = 1/nu * eye(p,p);
  arma::mat nIp = nu* n * eye(p,p);
  arma::mat xt = x.t();
  arma::mat matinv = zeros(n0,n0);
  
  
  arma::mat DB = zeros(p,p);
  arma::mat AB = zeros(n0, p);
  arma::mat mati = zeros(p,p);
  
  arma::uvec indexi;
  int idp1 = 0;
  int idp2 = p - 1;
  for(int i=0; i < n; i++){
    
    indexi = find(indexy == uindexy(i));
    
    mati = inv(xt.cols(indexi) * x.rows(indexi) + nIp);
    DB = DB + mati;
    AB.rows(idp1,idp2) = mati;
    matinv.submat(idp1,idp1,idp2,idp2) = mati; // output
    
    idp1 = idp1 + p;
    idp2 = idp2 + p;
  }
  
  
  arma::mat IB = inv(Ip - DB);
  
  matinv = matinv + AB * IB * AB.t();
  return(matinv);
  
}




#endif