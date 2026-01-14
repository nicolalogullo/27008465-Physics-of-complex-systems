// qkr_timeevolve.cpp
// Build the Floquet operator by explicitly time-evolving the single-period
// Schrödinger equation U(T,0) using small time steps.
// Then compute its eigenvalues (quasienergies).
//
// Requires Eigen (header-only): https://eigen.tuxfamily.org
//
// Compile:
//   g++ -O2 -std=c++17 -I /path/to/eigen qkr_timeevolve.cpp -o qkr_timeevolve
//
// Run example:
//   ./qkr_timeevolve 40 10.0 1.0 200
//
// Model:
//   H(t) = p^2/2I + K cos(theta) * sum_n delta(t - nT)
//   We take I = 1, T = 1. Free evolution for t in (0,1), delta kick at t=0.
//   Basis: momentum |m>, m = -M,...,M (dimension N=2M+1)

#include <Eigen/Dense>
#include <complex>
#include <iostream>
#include <vector>
#include <fstream>

using namespace Eigen;
using namespace std;

int main(int argc, char** argv) {
    if(argc < 5){
        cerr << "Usage: ./qkr_timeevolve M K tau Nt\n";
        cerr << "M: momentum cutoff => N = 2M+1\n";
        cerr << "K: kick strength\n";
        cerr << "tau: free propagation time (period = 1)\n";
        cerr << "Nt: number of time slices for free evolution\n";
        return 1;
    }

    int M = atoi(argv[1]);
    double K = atof(argv[2]);
    double tau = atof(argv[3]);
    int Nt = atoi(argv[4]);

    int N = 2*M + 1;
    using cplx = complex<double>;
    const cplx I(0.0,1.0);

    // Momentum indices m = -M,...,+M
    vector<int> m_vals(N);
    for(int i=0;i<N;i++) m_vals[i] = i - M;

    // Position grid for finite Fourier transform
    int Nx = N;
    vector<double> xgrid(Nx);
    for(int j=0;j<Nx;j++)
        xgrid[j] = 2.0*M_PI * double(j)/Nx;

    // Precompute Fourier transform matrices: momentum <-> position
    MatrixXcd F(N, N), Finv(N, N);
    for(int j=0;j<N;j++){
        for(int k=0;k<N;k++){
            double phase = -2.0*M_PI * double(m_vals[k]) * double(j) / double(N);
            F(j,k) = exp(I * phase)/sqrt(double(N));
            // Inverse transform
            double phase2 = +2.0*M_PI * double(m_vals[j]) * double(k) / double(N);
            Finv(j,k) = exp(I * phase2)/sqrt(double(N));
        }
    }

    // Build the kick operator U_kick in position basis: exp(-i K cos(x))
    VectorXcd phaseKick_x(N);
    for(int j=0;j<N;j++){
        phaseKick_x[j] = exp(-I * K * cos(xgrid[j]));
    }
    // Convert to matrix in momentum basis: U_k = F^-1 diag * F
    MatrixXcd U_kick(N,N);
    {
        MatrixXcd D = phaseKick_x.asDiagonal();
        U_kick = Finv * D * F;
    }

    // Build free evolution operator over small dt:
    double dt = tau / double(Nt);
    VectorXcd freePhase(N);
    for(int i=0;i<N;i++){
        double m = m_vals[i];
        double E = 0.5*m*m; // I=1
        freePhase[i] = exp(-I * E * dt);
    }
    MatrixXcd U_free_dt = freePhase.asDiagonal();
    MatrixXcd U_free = MatrixXcd::Identity(N,N);
    for(int n=0;n<Nt;n++){
        U_free = U_free_dt * U_free;
    }

    // Full Floquet operator: first kick at t=0, then free evolution
    MatrixXcd U_F = U_free * U_kick;

    // Diagonalize U_F
    ComplexEigenSolver<MatrixXcd> ces;
    ces.compute(U_F, true);
    VectorXcd evals = ces.eigenvalues();

    ofstream out("quasienergies.csv");
    out << "index,theta\n";
    for(int i=0;i<N;i++){
        cplx l = evals[i];
        double theta = atan2(imag(l), real(l));
        out << i << "," << theta << "\n";
    }
    out.close();

    return 0;
}