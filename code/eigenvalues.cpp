#include "eigenvalues.hpp"
#include <iostream>
#include <fstream>
#include "time.h"
#include <fstream>
#define TRUE 1
#define FALSE 0

void eigenvalues::solve(double tolerance, int maxiter){
offdiag();
Jacobi_rotate();
iterations = 1;

runtime = 0;
start = clock();
for(int i = 0;i<2;i++){
  for(int j=0;j<2;j++){
    if(j==i){
      R_eigen(i,j)=1;
    }
    else{
      R_eigen(i,j)=0;
    }
  }
}
finish = clock();
runtime = ( (finish - start)*1./CLOCKS_PER_SEC  );
//R_eigen.print();


while ( m_max > tolerance && iterations <= maxiter)
{
   offdiag();
   Jacobi_rotate();
   iterations++;
}
//R_eigen = 1*R_eigen;
//cout << "iterations:" <<  iterations << endl;
}


void eigenvalues::offdiag(){
   double max;
   for (int i = 0; i < n; ++i)
   {
       for ( int j = i+1; j < n; ++j)
       {
           double aij = fabs(A(i,j));
           if ( aij > max)
           {
              max = aij;  p = i; q = j;
           }
       }
   }
   m_max = max;
   }

void eigenvalues::Jacobi_rotate()
{
  int k = p;
  int l = q;

  if ( A(k,l) != 0.0 ) {
    tau = (A(l,l) - A(k,k))/(2*A(k,l));

    if ( tau >= 0 ) {
      t = 1.0/(tau + sqrt(1.0 + tau*tau));
    }
    else {
      t = -1.0/(-tau +sqrt(1.0 + tau*tau));
    }

    c = 1/sqrt(1+t*t);
    s = c*t;
  }
  else {
    c = 1.0;
    s = 0.0;
  }

  a_kk = A(k,k);
  a_ll = A(l,l);
  A(k,k) = c*c*a_kk - 2.0*c*s*A(k,l) + s*s*a_ll;
  A(l,l) = s*s*a_kk + 2.0*c*s*A(k,l) + c*c*a_ll;
  A(k,l) = 0.0;  // hard-coding non-diagonal elements by hand
  A(l,k) = 0.0;  // same here

  for ( int i = 0; i < n; i++ ) {
    if ( i != k && i != l ) {
      a_ik = A(i,k);
      a_il = A(i,l);
      A(i,k) = c*a_ik - s*a_il;
      A(k,i) = A(i,k);
      A(i,l) = c*a_il + s*a_ik;
      A(l,i) = A(i,l);
    }
//  And finally the new eigenvectors
    r_ik = R_eigen(i,k);
    r_il = R_eigen(i,l);

    R_eigen(i,k) = c*r_ik - s*r_il;
    R_eigen(i,l) = c*r_il + s*r_ik;
  }
  }


void eigenvalues::order_eigenvalues(){
  running = TRUE;
  col_swap = TRUE;
  while(running==TRUE){
    for(int i=0;i<n-1;i++){
      if(A(i+1,i+1)<A(i,i)){
        float a_ip1ip1 = A(i+1,i+1);
        float a_ii = A(i,i);
        A(i,i) = a_ip1ip1;
        A(i+1,i+1) = a_ii;
        R_eigen.swap_cols(i+1,i);
        col_swap = TRUE;
      }
      }
    if(col_swap == FALSE){
      running = FALSE;
    }
    col_swap = FALSE;
    }
}


vec eigenvalues::get_eigenvectors(int n_){
  vec v = zeros(n);
  for(int i=0;i<n;i++){
    v(i) = R_eigen(i,n_);
  }
  return v;
}

float eigenvalues::get_eigenvalues(int n_){
  return A(n_,n_);
}

void eigenvalues::QR_GS(){
  //Gram-Scmidt method
  mat Q = zeros(n,n);
  Q.col(0) = A.col(0)*1./norm(A.col(0));

  for(int k = 1;k<n;k++){
    vec proj = zeros(n);
    vec v_k = A.col(k);

    for(int j=0;j<k-1;j++){
      vec u_j = Q.col(j);
      proj = proj + ( dot(v_k,u_j)*1./dot(u_j,u_j) ) * u_j;
    }
    Q.col(k) = v_k - proj;
    Q.col(k) = Q.col(k)*1./norm(Q.col(k),2);


  }
  R_eigen = Q.i()*A;
  mat D = Q.t()*A*Q;

  D.print();

}
/*
void eigenvalues::Lanczos(){
  //NEED NON ZERO GUESS FOR r_0
  mat Q = zeros(n,n);
  //mat R = zeros(n,n);

  int k=0;
  mat I = eye(n,n);

  vec beta = zeros(n);
  beta(0) = 1;
  while(beta(k) != 0){
    Q.col(k+1) = R.col(k)*1./beta(k);
    k +=1;

    float alpha_k =  dot(Q.col(k).t(),A*Q.col(k));

    R.col(k) = (A-alpha_k*I)*Q.col(k)-beta(k-1)*Q.col(k-1);
    beta(k) = norm(R.col(k),2);
  }
  //Q.print();

  return;
}
*/
