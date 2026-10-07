#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <random>

/* Particle state:
x, y   - coordinates [fm]
vx, vy - velocity [c], c = 1 */

struct State {
    double x;
    double y;
    double vx;
    double vy;
};

/* State derivative:
dx/dt  = vx
dy/dt  = vy
dvx/dt = ax
dvy/dt = ay */

struct Derivative {
    double dx;
    double dy;
    double dvx;
    double dvy;
};


/* ==========================================================
Physical constants
============================================================*/

// e^2 / (4*pi*epsilon_0), Coulomb conversion factor [MeV*fm]
constexpr double eSquared = 1.439964;

// Charge numbers
constexpr double Z_ALPHA = 2.0;
constexpr double Z_GOLD  = 79.0;

// Alpha particle mass [MeV/c^2]
constexpr double MASS_ALPHA = 3727.38;

// k = Z1*Z2*e^2/(4*pi*epsilon_0) [MeV*fm]
constexpr double K =
    eSquared * Z_ALPHA * Z_GOLD;


/* ==========================================================
Equations of motion
============================================================*/

Derivative equations(const State& s)
{
    double r2 = s.x * s.x + s.y * s.y;
    double r  = std::sqrt(r2);
    double r3 = r2 * r;

    /*a = F/m
    F_x = K*x/r^3
    F_y = K*y/r^3 */

    double ax = (K / MASS_ALPHA) * s.x / r3;
    double ay = (K / MASS_ALPHA) * s.y / r3;

    return {
        s.vx,
        s.vy,
        ax,
        ay
    };
}


/* ==========================================================
Function for RK4
============================================================*/

State add(const State& s,
          const Derivative& d,
          double factor)
{
    return {
        s.x  + factor * d.dx,
        s.y  + factor * d.dy,
        s.vx + factor * d.dvx,
        s.vy + factor * d.dvy
    };
}


/* ==========================================================
One Runge-Kutta-4 step
============================================================*/

State rk4Step(const State& s, double dt)
{
    // K1 = dt * f(Y_i)
    Derivative f1 = equations(s);

    Derivative K1 {
        dt * f1.dx,
        dt * f1.dy,
        dt * f1.dvx,
        dt * f1.dvy
    };


    // K2 = dt * f(Y_i + K1/2)
    State s2 = add(s, K1, 0.5);

    Derivative f2 = equations(s2);

    Derivative K2 {
        dt * f2.dx,
        dt * f2.dy,
        dt * f2.dvx,
        dt * f2.dvy
    };


    // K3 = dt * f(Y_i + K2/2)
    State s3 = add(s, K2, 0.5);

    Derivative f3 = equations(s3);

    Derivative K3 {
        dt * f3.dx,
        dt * f3.dy,
        dt * f3.dvx,
        dt * f3.dvy
    };


    // K4 = dt * f(Y_i + K3)
    State s4 = add(s, K3, 1.0);

    Derivative f4 = equations(s4);

    Derivative K4 {
        dt * f4.dx,
        dt * f4.dy,
        dt * f4.dvx,
        dt * f4.dvy
    };


    // Delta Y = (K1 + 2K2 + 2K3 + K4) / 6
    State result;

    result.x =
        s.x +
        (K1.dx + 2.0*K2.dx + 2.0*K3.dx + K4.dx) / 6.0;

    result.y =
        s.y +
        (K1.dy + 2.0*K2.dy + 2.0*K3.dy + K4.dy) / 6.0;

    result.vx =
        s.vx +
        (K1.dvx + 2.0*K2.dvx + 2.0*K3.dvx + K4.dvx) / 6.0;

    result.vy =
        s.vy +
        (K1.dvy + 2.0*K2.dvy + 2.0*K3.dvy + K4.dvy) / 6.0;

    return result;
}


/* ==========================================================
Total mechanical energy
============================================================*/

double totalEnergy(const State& s)
{
    double r = std::sqrt(s.x*s.x + s.y*s.y);

    double v2 =
        s.vx*s.vx + s.vy*s.vy;

    double kinetic =
        0.5 * MASS_ALPHA * v2;

    double potential =
        K / r;

    return kinetic + potential;
}


/* ==========================================================
Angular momentum L_z
============================================================*/

double angularMomentum(const State& s)
{
    return MASS_ALPHA *
           (s.x*s.vy - s.y*s.vx);
}

/* ==========================================================
Scattering simulation
============================================================*/

double simulateScattering(double b)
{
    const double E = 5.0;         // MeV
    const double x0 = 10000.0;    // fm
    const double dt = 1.0;        // fm/c

    // Start to the left of the nucleus

    double initialX = -x0;
    double initialY = b;

    double r0 =
        std::sqrt(initialX * initialX +
                  initialY * initialY);

    /* Initial velocity

    E = 1/2 m v^2 + K/r0
    v = sqrt(2(E - K/r0)/m) */

    double potential0 = K / r0;

    double kinetic0 = E - potential0;

    double v0 =
        std::sqrt(2.0 * kinetic0 / MASS_ALPHA);

    State state {
        initialX,
        initialY,
        v0,
        0.0
    };


    // Simulation loop

    const int MAX_STEPS = 1000000;

    for (int i = 0; i < MAX_STEPS; ++i)
    {
        double r =
            std::sqrt(state.x * state.x +
                      state.y * state.y);

        double radialVelocity =
            state.x * state.vx +
            state.y * state.vy;

        /* The particle starts at approximately r = x0.
        Stop when it reaches the same distance again after scattering. */

        if (r >= x0 &&
            radialVelocity > 0.0)
        {
            break;
        }

        // One RK4 step

        state = rk4Step(state, dt);
    }


    // Numerical scattering angle

    double theta =
        std::atan2(state.vy, state.vx);

    return theta;
}


// MAIN

int main()
{
    // Monte Carlo parameters

    const int N = 5000;
    const double bMax = 100.0;


    /* ========================================================
    Random number generator
    ========================================================*/

    std::random_device rd;

    std::mt19937 generator(rd());

    std::uniform_real_distribution<double> uniform(0.0, 1.0);


    // Output file

    std::ofstream file("data/monte_carlo_scattering.csv");

    file << "event,U,b,theta_numerical,theta_theory\n";


    /* ========================================================
    Monte Carlo simulation
    ========================================================*/

    const double E = 5.0;
    const double pi = std::acos(-1.0);

    for (int i = 0; i < N; ++i)
    {
        // Generate U uniformly in [0, 1)

        double U = uniform(generator);


        // Sample impact parameter

        double b = bMax * std::sqrt(U);


        // Numerical scattering angle from RK4

        double thetaNumerical =
            simulateScattering(b);


        // Analytical Rutherford angle

        double thetaTheory =
            2.0 * std::atan(
                K / (2.0 * E * b)
            );


        // Convert radians to degrees

        double thetaNumericalDeg =
            thetaNumerical * 180.0 / pi;

        double thetaTheoryDeg =
            thetaTheory * 180.0 / pi;


        // Save event

        file
            << i << ","
            << U << ","
            << b << ","
            << thetaNumericalDeg << ","
            << thetaTheoryDeg << "\n";


        // Print progress every 100 events

        if ((i + 1) % 100 == 0)
        {
            std::cout
                << "Completed "
                << i + 1
                << " / "
                << N
                << " events\n";
        }
    }


    file.close();


    std::cout
        << "\nMonte Carlo scattering simulation completed.\n";

    std::cout
        << "Results saved to data/monte_carlo_scattering.csv\n";


    return 0;
}
