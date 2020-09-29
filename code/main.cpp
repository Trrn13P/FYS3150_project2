#include <iostream>
#include <armadillo>
#include "eigenvalues.hpp"
#include <fstream>
#include <cmath>
#include "time.h"

using namespace arma;
using namespace std;

float V_one_electron(float rho_){
  return rho_ * rho_;
}
float V_two_electron(float omega_, float rho_){
  return omega_ * omega_ * rho_ * rho_ + 1./rho_;
}

mat A_notQ(int n){
  int N = n+1;
  float h = 1./N;
  float a = -1*1./pow(h,2);
  float d = 2*1./pow(h,2);

  mat A = zeros(n,n);
  for(int i=0;i<n;i++){
    if (i!=0){
        A(i,i-1) = a;
      }
    A(i,i) = d;
    if(i!=n-1){
        A(i,i+1) = a;
      }
    }
    return A;
}


mat A_quantum(int n, float rho_0, float rho_N, int electron_number, float omega_r){
  mat A = zeros(n,n);
  float h = (rho_N-rho_0)*1./(n-1);

  vec d = zeros(n);
  float e = -1*1./(h*h);
  vec rho = zeros(n);

  for(int i=0;i<n;i++){
    rho(i) = rho_0+i*h;
    if(electron_number==1){
      d(i) = 2*1./(h*h) + V_one_electron(rho(i));
    }
    if(electron_number==2){
    d(i) = 2*1./(h*h) + V_two_electron(omega_r,rho(i));
  }

    if (i!=0){
        A(i,i-1) = e;
      }
    A(i,i) = d(i);
    if(i!=n-1){
        A(i,i+1) = e;
      }
    }
  return A;
}


void one_electron_print(int n,float rho_0,float rho_N,float omega_r, double tolerance, int maxiter){
  int electron_number = 1;
  mat A = A_quantum(n,rho_0,rho_N,electron_number, omega_r);
  eigenvalues test(A,n);
  test.solve(tolerance,maxiter);
  test.order_eigenvalues();
  /*
  mat B = test.get_solution(0,rho_0,rho_N);

  string filename = "../data/one_electron/N"+to_string(n)+"_jacobi.txt";
  ofstream outfile(filename);
  outfile << "tolerance=" << tolerance << " ";
  outfile << "iterations="<< test.iterations << " n=" << n << " n_electrons=" << electron_number<< endl;
  outfile << "rho:\n";
  for(int i=0;i<n+2;i++){
    outfile << B(i,0)<<" ";
  }
  outfile << endl;
  outfile << "eigenvalue: " << "eigenvector: \n";
  for(int i=0;i<10;i++){
    outfile << test.get_eigenvalues(i);
    B = test.get_solution(i,rho_0,rho_N);
    outfile << B.col(1).t();
  }
  outfile.close();
  */
}

void two_electron_print(int n, double tolerance, int maxiter){
  return;
}

void buck_beam_print(double tolerance, int maxiter){
  int n = 100;
  mat A = A_notQ(n);
  eigenvalues test(A,n);
  test.solve(tolerance,maxiter);
  test.order_eigenvalues();
  string filename = "../data/buck_beam/N"+to_string(n)+"_jacobi.txt";
  ofstream outfile(filename);
  outfile << "tolerance=" << tolerance << " ";
  outfile << "iterations="<< test.iterations << " n=" << n << endl;

  vec analytic = zeros(n);
  int N = n+1;
  for(int i=0;i<n-1;i++){
    analytic(i) = sin((i)*M_PI/N);
  }

  outfile << "Eigenvalue_0: " << test.get_eigenvalues(0) << endl;
  outfile << "rho: u_a: u_n:" << endl;
  vec rho = zeros(n);
  for(int i=0;i<n-1;i++){
    rho(i+1)=1./(n-1)+rho(i);
  }
  outfile << "rho: "<< rho.t();
  outfile << "u_a: " << analytic.t();
  outfile << "u_n: " << test.get_eigenvectors(0).t();
outfile.close();
}

int main(int argc, char const *argv[]) {
  float start1, finish1, runtime1;
  float start2, finish2, runtime2;


  double tolerance = 1.0E-6;
  int maxiter = 1000000;
  int n_;
  //buck_beam_print(tolerance, maxiter);
  int arr[7] = {10,50,100,200,300,400,500};


  string filename = "../data/CPU_TIMES.txt";
  ofstream outfile(filename);

  for(int i=0;i<7;i++){
    n_ = arr[i];
    mat A = A_notQ(n_);
    eigenvalues test(A,n_);

    runtime1 = 0;
    start1 = clock();

    test.solve(tolerance,maxiter);
    finish1 = clock();
    runtime1 = ( (finish1 - start1)*1./CLOCKS_PER_SEC );

    test.order_eigenvalues();


    runtime2 = 0;
    start2 = clock();
    vec eigval =  eig_sym(A);
    float lambda_0 = eigval(0);

    finish2 = clock();
    runtime2 = ( (finish2 - start2)*1./CLOCKS_PER_SEC );
    outfile << "n=" << n_ << " Runtime_arma=" << runtime2 << " runtime jacobi=" << runtime1
    << " eigenval_arma=" << lambda_0 << " eigenval_jacobi=" << test.get_eigenvalues(0) << endl;
  }
  outfile.close();



/*
plot computing times for different n's for Jacobi_rotate and the other method
Og se paa mean error mellom eigenvalues man finner paa forskjellige
metoder og analytiske
*/


  /*
  //Printer jacobi for 1 elektron
  int n = 10;
  double tolerance = 1.0E-6;
  int maxiter = 1000000;
  float rho_0 = 0;
  float rho_N = 4;
  int electron_number = 1;
  float omega_r = 1;
  one_electron_print(n,rho_0,rho_N,omega_r, tolerance, maxiter);
  */

  /*
  int n = 100;
  double tolerance = 1.0E-6;
  int maxiter = 1000000;
  double rho_0 = 1E-8;
  float rho_N = 4;
  int electron_number = 2;
  float omega_r = 1;
*/

  /*
  mat A = A_quantum(n,rho_0,rho_N,electron_number,omega_r);
  eigenvalues test(A,n);
  //A.print();
  test.solve(tolerance,maxiter);
  test.order_eigenvalues();
  cout << test.get_eigenvalues(0) << endl;
  */
  //mat A1 = A_notQ(n);
  //eigenvalues test(A1,n);
  //float rho_0 = 0;
  //float rho_N = 1;
  //mat A2 = A_quantum(n,rho_0,rho_N);




  //test.solve(tolerance,maxiter);
  //test.QR_GS();
  //test.Lanczos();





}
