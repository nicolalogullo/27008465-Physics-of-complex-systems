// qkr_floquet.cpp
// Build: g++ -O2 -std=c++17 -I /path/to/eigen qkr_floquet.cpp -o qkr_floquet
//
// Example run: ./qkr_floquet 50 10.0 1.0 42
//   M = 50  -> matrix size N = 101
//   K = 10.0
//   kbar = 1.0
//   seed = 42 (0 for none)
//
// Produces CSV: qkr_output_M50_K10000_kbar1000.csv

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <complex>
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <random>

using namespace Eigen;
using namespace std;
using dcomplex = complex<double>;

// Note: uses std::cyl_bessel_j (C++17). If your standard library doesn't provide it,
// you can replace with a Bessel implementation or use std::boost::math::cyl_bessel_j.

int main(int argc, char** argv){
    if(argc<5){
        cout << "Usage: "<< argv[0] << " M K kbar seed\n";
        cout << " M: momentum cutoff (half-size), matrix size N=2M+1\n";
        cout << " K: kick strength\n";
        cout << " kbar: effective hbar (kbar)\n";
        cout << " seed: integer seed for small phase averaging (0 for none)\n";
        return 1;
    }
    int M = atoi(argv[1]);
    double K = atof(argv[2]);
    double kbar = atof(argv[3]);
    int seed = atoi(argv[4]);

    int N = 2*M+1;
    MatrixXcd U(N,N);
    U.setZero();

    // momentum indices m = -M..M mapped to 0..N-1
    auto idx = [&](int m){ return m + M; };

    // optionally add a small twist phase to break accidental symmetries; controlled by seed
    double twist = 0.0;
    if(seed!=0){
        std::mt19937 rng(seed);
        std::uniform_real_distribution<double> unif(0.0, 2.0*M_PI);
        twist = unif(rng) * 1e-3; // small
    }

    // Precompute Bessel coefficients J_{n}(K/kbar)
    int maxj = 2*M; // needed range for m'-m
    double arg = K / kbar;

    // Construct U_{m',m} = e^{-i * m^2 * kbar / 2 + i * twist * m} * (-i)^{m'-m} * J_{m'-m}(K/kbar)
    const double PI = acos(-1.0);
    for(int mp=-M; mp<=M; ++mp){
        for(int m=-M; m<=M; ++m){
            int n = mp - m; // order of Bessel
            // compute J_n(arg)
            double Jn;
            if(n >= 0){
                Jn = std::cyl_bessel_j(n, arg);
            } else {
                int nn = -n;
                double Jnn = std::cyl_bessel_j(nn, arg);
                // J_{-n}(x) = (-1)^n J_n(x) for integer n
                Jn = ((nn % 2) == 0) ? Jnn : -Jnn;
            }

            // factor (-i)^{n} = exp(-i * pi/2 * n)
            dcomplex pref = std::exp(dcomplex(0.0, -PI/2.0 * n));
            double phase_free = -0.5 * m * m * kbar;
            dcomplex phase = std::exp(dcomplex(0.0, phase_free + twist * m));
            U(idx(mp), idx(m)) = phase * pref * Jn;
        }
    }

    // diagonalize U (unitary) using Eigen's ComplexEigenSolver
    ComplexEigenSolver<MatrixXcd> ces;
    ces.compute(U, /* computeEigenvectors = */ true);
    if(ces.info() != Success){
        cerr << "Eigen decomposition failed\n";
        return 2;
    }

    VectorXcd evals = ces.eigenvalues();
    // extract phases in [0,2pi)
    vector<double> phases; phases.reserve(N);
    for(int i=0;i<evals.size();++i){
        double th = std::arg(evals(i));
        if(th < 0) th += 2.0*PI;
        phases.push_back(th);
    }

    sort(phases.begin(), phases.end());

    // compute spacings (including wraparound)
    vector<double> spacings;
    spacings.reserve(N);
    for(int i=0;i<N-1;++i) spacings.push_back(phases[i+1]-phases[i]);
    spacings.push_back(2.0*PI - phases.back() + phases.front());

    // mean spacing
    double mean_s = 0.0;
    for(double s: spacings) mean_s += s;
    mean_s /= spacings.size();

    // normalized spacings
    vector<double> s_norm; s_norm.reserve(spacings.size());
    for(double s: spacings) s_norm.push_back(s / mean_s);

    // r-statistics (consecutive spacing ratios)
    vector<double> rvals;
    for(int i=0;i<(int)s_norm.size();++i){
        double s_i = s_norm[i];
        double s_prev = s_norm[(i-1 + s_norm.size()) % s_norm.size()];
        double r = min(s_i, s_prev) / max(s_i, s_prev);
        rvals.push_back(r);
    }

    double rmean = 0.0; for(double r: rvals) rmean += r; rmean /= rvals.size();

    // output to CSV: phases_unfolded (normalized), s_norm, rvals
    // create a filename encoding parameters (rounded)
    auto round_int = [](double x){ return static_cast<int>(std::round(x*1000.0)); };
    string fname = "qkr_output_M" + to_string(M) + "_K" + to_string(round_int(K)) + "_kbar" + to_string(round_int(kbar)) + ".csv";
    ofstream fout(fname);
    fout << "phase,spacing_norm,r\n";
    for(int i=0;i<N;++i){
        fout << phases[i]/(2.0*PI) << "," << s_norm[i] << "," << rvals[i] << "\n";
    }
    fout.close();

    cout << "Wrote "<<fname<<" (N="<<N<<"). Mean normalized spacing=1 (by construction).\n";
    cout << "Mean r = " << rmean << "\n";
    cout << "To visualize, plot spacing histogram of column spacing_norm and compare with Poisson/Wigner surmise.\n";

    return 0;
}
