#include <iostream>
#include <fstream>
#include <cmath>

// Constants
constexpr double mu = 398600.4418; // km^3/s^2
constexpr int K = 4;               // AM-4
constexpr double dt = 10.0;        // time step (s)
constexpr int steps = 10000;       // number of steps to simulate

// AM-4 coefficients (predictor = x + dt * sum(ak[i] * dx[i]))
constexpr double ak[3] = {-1.0/12, 2.0/3, 5.0/12};
constexpr double aK = 3.0/8;

// State history (t, x, dx): x = {x, y, z, vx, vy, vz}
double t[K];
double x[K][6];   // position + velocity
double dx[K][6];  // time derivatives

// Compute gravitational acceleration at a point in space
void compute_derivative(double t_now, const double* x_now, double* dx_out) {
    double rx = x_now[0], ry = x_now[1], rz = x_now[2];
    double vx = x_now[3], vy = x_now[4], vz = x_now[5];

    double r3 = std::pow(rx*rx + ry*ry + rz*rz, 1.5);

    dx_out[0] = vx;
    dx_out[1] = vy;
    dx_out[2] = vz;
    dx_out[3] = -mu * rx / r3;
    dx_out[4] = -mu * ry / r3;
    dx_out[5] = -mu * rz / r3;
}

// Initialize first K points using RK4 (simpler than recursive AM bootstrapping)
void rk4_step(double t_now, const double* x_now, double* x_out) {
    double k1[6], k2[6], k3[6], k4[6], temp[6];

    compute_derivative(t_now, x_now, k1);
    for (int i = 0; i < 6; ++i) temp[i] = x_now[i] + 0.5 * dt * k1[i];
    compute_derivative(t_now + 0.5*dt, temp, k2);
    for (int i = 0; i < 6; ++i) temp[i] = x_now[i] + 0.5 * dt * k2[i];
    compute_derivative(t_now + 0.5*dt, temp, k3);
    for (int i = 0; i < 6; ++i) temp[i] = x_now[i] + dt * k3[i];
    compute_derivative(t_now + dt, temp, k4);

    for (int i = 0; i < 6; ++i)
        x_out[i] = x_now[i] + (dt/6.0)*(k1[i] + 2*k2[i] + 2*k3[i] + k4[i]);
}

// Perform one Adams–Moulton step
void drive() {
    double t_next = t[K-1] + dt;
    double pred[6] = {};

    // Predictor step: x_new = x + dt * sum(ak * dx)
    for (int j = 0; j < 6; ++j) {
        pred[j] = x[K-1][j];
        for (int i = 0; i < K - 1; ++i)
            pred[j] += dt * ak[i] * dx[i][j];
    }

    // Corrector step: add dt * aK * dx_new
    double d_next[6];
    compute_derivative(t_next, pred, d_next);
    for (int j = 0; j < 6; ++j)
        pred[j] += dt * aK * d_next[j];

    // Shift state history
    for (int i = 0; i < K - 1; ++i) {
        t[i] = t[i+1];
        for (int j = 0; j < 6; ++j) {
            x[i][j] = x[i+1][j];
            dx[i][j] = dx[i+1][j];
        }
    }
    t[K-1] = t_next;
    for (int j = 0; j < 6; ++j) {
        x[K-1][j] = pred[j];
        dx[K-1][j] = d_next[j];
    }
}

int main() {
    // Initial circular orbit at 7000 km altitude
    double r0 = 7000.0;                      // km
    double v0 = std::sqrt(mu / r0);         // km/s
    double state0[6] = {r0, 0, 0, 0, v0, 0}; // [x,y,z,vx,vy,vz]

    // Initialize first point
    t[0] = 0.0;
    for (int j = 0; j < 6; ++j) x[0][j] = state0[j];
    compute_derivative(t[0], x[0], dx[0]);

    // Fill initial K using RK4
    for (int i = 1; i < K; ++i) {
        t[i] = t[i-1] + dt;
        rk4_step(t[i-1], x[i-1], x[i]);
        compute_derivative(t[i], x[i], dx[i]);
    }

    // Open CSV output
    std::ofstream out("eci_orbit.csv");
    out << "t,x,y,z,vx,vy,vz\n";
    for (int i = 0; i < K; ++i)
        out << t[i] << "," << x[i][0] << "," << x[i][1] << "," << x[i][2] << ","
            << x[i][3] << "," << x[i][4] << "," << x[i][5] << "\n";

    // Main propagation loop
    for (int step = 0; step < steps; ++step) {
        drive();
        out << t[K-1] << "," << x[K-1][0] << "," << x[K-1][1] << "," << x[K-1][2] << ","
            << x[K-1][3] << "," << x[K-1][4] << "," << x[K-1][5] << "\n";
    }

    std::cout << "Orbit propagation complete. Output: eci_orbit.csv\n";
    return 0;
}
