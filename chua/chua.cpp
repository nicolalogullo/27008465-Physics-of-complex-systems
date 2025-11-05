#include <iostream>
#include <fstream>
#include <cmath>

struct params {
    double alpha, beta, gamma, a, b;
};

// Non-linear resistor potential drop
double nr_pot(double v, const params& p) {
    //if (v > 1.0)  return p.b * v + (p.a - p.b);
    //if (v < -1.0) return p.b * v - (p.a - p.b);
    return p.b * v + 0.5*( p.a - p.b ) * ( abs( v + 1.0 ) - abs( v - 1.0) );
}

// Right-hand side of the system of differential equations
void rhs(const double X[3], double dX[3], const params& p) {
    double x = X[0], y = X[1], z = X[2];
    dX[0] = p.alpha * (y - x - nr_pot(x, p));
    dX[1] = x - y + z;
    dX[2] = -p.beta * y - p.gamma * z;
}

// RK4 integrator
void rk4(double X[3], double dt, const params& p) {
    double k1[3], k2[3], k3[3], k4[3], tmp[3];
    rhs(X, k1, p);
    for (int i=0;i<3;i++) tmp[i]=X[i]+0.5*dt*k1[i];
    rhs(tmp, k2, p);
    for (int i=0;i<3;i++) tmp[i]=X[i]+0.5*dt*k2[i];
    rhs(tmp, k3, p);
    for (int i=0;i<3;i++) tmp[i]=X[i]+dt*k3[i];
    rhs(tmp, k4, p);
    for (int i=0;i<3;i++)
        X[i]+=dt/6.0*(k1[i]+2*k2[i]+2*k3[i]+k4[i]);
}

int main() {
    //params p = {C2/C1,R*R*C2/L ,R*R0*C2/L,R*Ga,R*Gb};
    params p = {1.0, 10.0, 0.05,-1.0,-1.0};
    //params p = {9.0, 14.286, 0.05,-1.143, -0.714};
    double X[3] = {0.1, 0.0, 0.0};
    double dt = 0.001;
    int N = 200000;

    std::ofstream traj("trajectory.txt");
    std::ofstream poinc("poincare.txt");

    double y_prev = X[1];
    for (int i = 0; i < N; ++i) {
        rk4(X, dt, p);
        double t = i * dt;

        if (i % 100 == 0)
            traj << t << " " << X[0] << " " << X[1] << " " << X[2] << "\n";

        // Poincaré section y=0 upward
        if (y_prev < 0.0 && X[1] >= 0.0)
            poinc << X[0] << " " << X[2] << "\n";
        y_prev = X[1];
    }

    std::cout << "Data saved: trajectory.txt, poincare.txt\n";
    return 0;
}

