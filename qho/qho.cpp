#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <eigen-5.0.0/Eigen/Dense>

using namespace Eigen;
using namespace std;

int main()
{
    // -----------------------------
    // Parameters
    // -----------------------------
    const int N = 80;          // HO basis size
    const int N_low = N / 2;   // keep low-energy part
    const double lambda = 0.01; //x^4 coefficient

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
    MatrixXd H = num + 0.5 * MatrixXd::Identity(N, N)
               + lambda * x4;

    // -----------------------------
    // Diagonalization
    // -----------------------------
    SelfAdjointEigenSolver<MatrixXd> solver(H);
    VectorXd E = solver.eigenvalues();

    // -----------------------------
    // Save spectrum
    // -----------------------------
    vector<double> levels(N_low);
    ofstream spec("spectrum.dat");

    for (int i = 0; i < N_low; ++i) {
        levels[i] = E(i);
        spec << levels[i] << "\n";
    }
    spec.close();

    // -----------------------------
    // Energy level spacing
    // -----------------------------
    vector<double> spacings;
    ofstream sp("spacings.dat");

    for (int i = 1; i < N_low-1; ++i) {
        double s = (levels[i+1] - levels[i]);
        spacings.push_back(s);
        sp << s << "\n";
    }
    sp.close();

    // -----------------------------
    // Output
    // -----------------------------
    cout << "Anharmonic oscillator spectrum computed\n";
    cout << "Basis size N = " << N << "\n";
    cout << "lambda = " << lambda << "\n";
    cout << "Unfolded spacings computed: " << spacings.size() << "\n";

    return 0;
}

