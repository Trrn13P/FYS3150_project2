#include <iostream>
#include <armadillo>
#include "eigenvalues.hpp"
#include <fstream>

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
  float a = -1*1./(h*h);
  float d = 2*1./(h*h);

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
  mat B = test.get_solution(0,rho_0,rho_N);

  string filename = "../data/one_electron/N"+to_string(n)+"_jacobi.txt";
  ofstream outfile(filename);
  outfile << "tolerance=" << tolerance << " ";
  outfile << "iterations="<< test.iterations << " n=" << n << "n_electrons=" << electron_number<< endl;
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

}

void two_electron_print(int n, double tolerance, int maxiter){
  return;
}



int main(int argc, char const *argv[]) {
/*
plot computing times for different n's for Jacobi_rotate and the other method
Og se paa mean error mellom eigenvalues man finner paa forskjellige
metoder og analytiske
*/

  int n = 10;
  double tolerance = 1.0E-6;
  int maxiter = 1000000;
  float rho_0 = 0;
  float rho_N = 4;
  int electron_number = 1;
  float omega_r = 1;
  one_electron_print(n,rho_0,rho_N,omega_r, tolerance, maxiter);

  //mat A1 = A_notQ(n);
  //eigenvalues test(A1,n);
  //float rho_0 = 0;
  //float rho_N = 1;
  //mat A2 = A_quantum(n,rho_0,rho_N);




  //test.solve(tolerance,maxiter);
  //test.QR_GS();
  //test.Lanczos();





}
