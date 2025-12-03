#include <iostream>
#include <cmath>
#include <fstream>
using std::cout;
using std::cos, std::sin;
using std::ofstream;

struct state {
    double th1, th2;     // angles
    double w1, w2;       // angular velocities
};

const double g = 9.81;
double m1 = 1.0, m2 = 1.0;
double l1 = 1.0, l2 = 1.0;

// Compute derivatives: dy/dt = f(y)
state rhs(const state& s) {
    double th1 = s.th1, th2 = s.th2;
    double w1 = s.w1, w2 = s.w2;

    double delta = th1 - th2;
    double den1 = l1*((m1 + m2)-m2*cos(delta)*cos(delta));
    double den2 = (l2/l1)*den1;

    state ds;
    ds.th1 = w1;
    ds.th2 = w2;

    ds.w1 =
        (m2*l1*w1*w1*sin(2*delta)
         +m2*l2*w2*w2*sin(delta)
         +m2*g*sin(th2)*cos(delta)
         -(m1+m2)*g*sin(th1))
        /den1;

    ds.w2 =
        (-ds.w1*cos(delta)
         -l1*w1*w1*sin(delta)
         -g*sin(th2)) 
         /den2;

    return ds;
}

// RK4 step
state rk4_step(const state& s, double dt) {
    state k1 = rhs(s);

    state s2 = { s.th1 + 0.5*dt*k1.th1,
                 s.th2 + 0.5*dt*k1.th2,
                 s.w1  + 0.5*dt*k1.w1,
                 s.w2  + 0.5*dt*k1.w2 };

    state k2 = rhs(s2);

    state s3 = { s.th1 + 0.5*dt*k2.th1,
                 s.th2 + 0.5*dt*k2.th2,
                 s.w1  + 0.5*dt*k2.w1,
                 s.w2  + 0.5*dt*k2.w2 };

    state k3 = rhs(s3);

    state s4 = { s.th1 + dt*k3.th1,
                 s.th2 + dt*k3.th2,
                 s.w1  + dt*k3.w1,
                 s.w2  + dt*k3.w2 };

    state k4 = rhs(s4);

    state out;
    out.th1 = s.th1 + dt*(k1.th1 + 2*k2.th1 + 2*k3.th1 + k4.th1)/6.0;
    out.th2 = s.th2 + dt*(k1.th2 + 2*k2.th2 + 2*k3.th2 + k4.th2)/6.0;
    out.w1  = s.w1  + dt*(k1.w1  + 2*k2.w1  + 2*k3.w1  + k4.w1 )/6.0;
    out.w2  = s.w2  + dt*(k1.w2  + 2*k2.w2  + 2*k3.w2  + k4.w2 )/6.0;

    return out;
}

int main() {
    state s = { 0.002, 0.001, 0.0, 0.0 }; // initial conditions
    double dt = 0.001;
    int steps = 50000;

    ofstream file("results.dat");

    for(int i = 0; i < steps; i++) {
        file << i*dt << " "
             << s.th1 << " " << s.th2 << " "
             << s.w1  << " " << s.w2  << "\n";

        s = rk4_step(s, dt); 
    }

    file.close();
    cout << "DONE!\n";
}

