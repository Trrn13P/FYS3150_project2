#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "catch.hpp"
#include "eigenvalues.hpp"
#include <armadillo>
#include <cmath>
#include <iostream>

using namespace std;
using namespace arma;

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


TEST_CASE("Frobenius norm"){
    int n = 4;
    mat A = A_notQ(n);


    double tolerance = 1.0E-10;
    int maxiter = 1000;

  float norm_before = 0;
  for(int i=0;i<n;i++){
    for(int j=0;j<n;j++){
      norm_before+=pow(A(i,j),2);
    }
  }
  norm_before = pow(norm_before,0.5);

  eigenvalues test(A,n);
  test.solve(tolerance,maxiter);

  float norm_after = 0;
  for(int i=0;i<n;i++){
    norm_after+= pow(test.get_eigenvalues(i),2);
  }
  norm_after = pow(norm_after,0.5);

  REQUIRE(norm_before==Approx(norm_after).epsilon(0.0001));
}


TEST_CASE("Testing eigenvalues"){
  int n = 4;
  int N = n+1;
  float h = 1./N;
  float a = -1*1./pow(h,2);
  float d = 2*1./pow(h,2);
  mat A = A_notQ(n);

    double tolerance = 1.0E-10;
    int maxiter = 1000;
      eigenvalues test(A,n);
      test.solve(tolerance,maxiter);
      test.order_eigenvalues();

      for(int j=0;j<n;j++){
        float eig1 = test.get_eigenvalues(j);
        float eig2 = d+2*a*cos((j+1)*M_PI*1./N);
        REQUIRE(eig1==Approx(eig2).epsilon(0.0001));
      }
}

TEST_CASE("Testing eigenvector orthogonality"){
  int n = 4;
  mat A = A_notQ(n);

    double tolerance = 1.0E-10;
    int maxiter = 1000;
    eigenvalues test(A,n);
    test.solve(tolerance,maxiter);

    float x;
    vec vec1,vec2;
    for(int i=0;i<n;i++){
      vec1 = test.get_eigenvectors(i);
      for(int j=0;j<n;j++){
        if(j!=i){
          vec2 = test.get_eigenvectors(j);
          x = dot(vec1,vec2);
          REQUIRE(x==Approx(0).epsilon(0.0001));
        }
      }
    }
  }
