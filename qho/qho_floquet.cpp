#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <eigen-5.0.0/Eigen/Dense>
#include <eigen-5.0.0/unsupported/Eigen/MatrixFunctions>
#include <complex>

using namespace Eigen;
using namespace std;

int main()
{
    // -----------------------------
    // Parameters
    // -----------------------------
    const int N = 80;          // HO basis size
    const int N_low = N / 2;   // keep low-energy part
    const double lambda = 0.01;
    const double F = 0.2;      // drive amplitude
    const double Omega = 1.2;  // drive frequency
    const int Nt = 500;        // time slices per period

    const double T = 2.0 * M_PI / Omega;
    const double dt = T / Nt;


    // -----------------------------
    // Annihilation operators
    // -----------------------------
    MatrixXd a = MatrixXd::Zero(N, N);
    MatrixXd adag = MatrixXd::Zero(N, N);

    for (int n = 0; n < N-1; ++n) {
        a(n, n+1) = sqrt(n + 1); //Annihilation operator
    }

    // -----------------------------
    // Number operator
    // -----------------------------
    MatrixXd num = a.transpose() * a;

    // -----------------------------
    // Position operator x
    // x = (a + a†)/sqrt(2)
    // -----------------------------
    MatrixXd x = (a + a.transpose()) / sqrt(2.0);

    // -----------------------------
    // Nonlinear term x^4
    // -----------------------------
    MatrixXd x2 = x * x;
    MatrixXd x4 = x2 * x2;

    // -----------------------------
    // Hamiltonian
    // -----------------------------
    MatrixXd H0 = num + 0.5 * MatrixXd::Identity(N, N)
               + lambda * x4;

    // -----------------------------
    // Floquet operator initialization
    // -----------------------------
    MatrixXcd U = MatrixXcd::Identity(N, N);

    // -----------------------------
    // Time evolution over one period
    // -----------------------------
    for (int k = 0; k < Nt; ++k) {
        double t = k * dt;
        MatrixXcd Ht = -complex<double>(0,1) * dt* (H0 + F * cos(Omega * t) * x);
        MatrixXcd step = Ht.exp();
        U = step * U;
    }

    // -----------------------------
    // Diagonalize Floquet operator
    // -----------------------------
    ComplexEigenSolver<MatrixXcd> solver(U);
    VectorXcd eigvals = solver.eigenvalues();

    // -----------------------------
    // Extract quasienergies
    // -----------------------------
    vector<double> quasienergy;

    for (int i = 0; i < N; ++i) {
        double phase = arg(eigvals(i));     // in (-π, π]
        double eps = phase / T;             // quasienergy
        quasienergy.push_back(eps);
    }

    // Sort quasienergies
    sort(quasienergy.begin(), quasienergy.end());

    // Keep central part (avoid edge effects)
    vector<double> levels;
    for (int i = 0; i < N_low; ++i)
        levels.push_back(quasienergy[i]);

    // -----------------------------
    // Save quasienergy spectrum
    // -----------------------------
    ofstream spec("quasienergy.dat");
    for (double e : levels) spec << e << "\n";
    spec.close();

    // -----------------------------
    // Level spacing statistics
    // -----------------------------
    ofstream sp("floquet_spacings.dat");

    vector<double> spacings;
    for (size_t i = 0; i < levels.size() - 1; ++i) {
        double s = levels[i+1] - levels[i];
        spacings.push_back(s);
        sp << s << "\n";
    }
    sp.close();

    // -----------------------------
    // Output
    // -----------------------------
    cout << "Floquet spectrum computed\n";
    cout << "Drive amplitude F = " << F << "\n";
    cout << "Drive frequency Omega = " << Omega << "\n";
    cout << "Quasienergy spacings computed: " << spacings.size() << "\n";

    return 0;
}

