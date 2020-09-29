#ifndef EIGENVALUES_HPP
#define EIGENVALUES_HPP

#include <armadillo>
using namespace arma;

class eigenvalues {
  private:
    int n;
    int N;
    double m_max;


    int p, q;

    double s, c;
    double t, tau;
    double a_kk, a_ll, a_ik, a_il, r_ik, r_il;

    Mat<double> R_eigen;
    Mat<double> A;
    mat Q;
    mat R_qr;

    bool running;
    bool col_swap;

    void Initialize(mat A_, int n_){
      n = n_;
      A = A_;
      R_eigen = zeros(n,n);
      for(int i=0;i<n;i++){
        R_eigen(i,i) = 1;
      }
      Q = zeros(n,n);
      R_qr = zeros(n,n);
    }
      //overload function
      void Initialize(){
        n = 4;
        A = zeros(n);
        for(int i=0;i<n;i++){
          if (i!=0){
              A(i,i-1) = -1;
            }
          A(i,i) = 2;
          if(i!=n-1){
              A(i,i+1) = -1;
            }
          }

        std::cout << "Running on overload function\n" << std::endl;
        Initialize(A,n);
      }


  public:
    int iterations;
    void offdiag();
    void Jacobi_rotate();
    void solve(double tolerance,int maxiter);

    vec get_eigenvectors(int n_);
    float get_eigenvalues(int n_);

    void order_eigenvalues();
    mat get_solution(int n_, float rho_0, float rho_N);

    void Lanczos();
    void QR_GS();

    //setting up the overload
    eigenvalues(mat A, int n_){
      Initialize(A,n_);
    }
    eigenvalues(){
      Initialize();
    }
};
#endif
