#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include <algorithm>


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
constexpr double Z_GOLD = 79.0;

// Alpha particle mass [MeV/c^2]
constexpr double MASS_ALPHA = 3727.38;

// k = Z1*Z2*e^2/(4*pi*epsilon_0) [MeV*fm]
constexpr double K =
    eSquared * Z_ALPHA * Z_GOLD;

// Gold nucleus mass number
constexpr double A_GOLD = 197.0;

// Nuclear radius coefficient [fm]
constexpr double r0 = 1.2;

// Nuclear radius R = r0 * A^(1/3) [fm]
const double R = r0 * std::cbrt(A_GOLD);

const double PI = std::acos(-1.0);


/* ==========================================================
Nuclear models
============================================================*/

enum class NuclearModel {
    PointNucleus,
    UniformSphere
};


/* ==========================================================
Potential energy U(r)
============================================================*/

double potentialEnergy(
    double r,
    NuclearModel model
)
{
    // Point nucleus: U(r) = K/r

    if (model == NuclearModel::PointNucleus) {
        return K / r;
    }

    // Outside uniform sphere: U(r) = K/r

    if (r >= R) {
        return K / r;
    }

    /* Inside uniform sphere:

    U(r) = K/(2R) * (3 - r^2/R^2) */

    return (K / (2.0 * R)) *
           (3.0 - r * r / (R * R));
}


/* ==========================================================
Equations of motion
============================================================*/

Derivative equations(
    const State& s,
    NuclearModel model
)
{
    double r2 = s.x * s.x + s.y * s.y;
    double r = std::sqrt(r2);

    double factor;

    if (model == NuclearModel::UniformSphere && r < R)
    {
        /* Inside uniform charged sphere:

        F_x = K*x/R^3
        F_y = K*y/R^3

        a = F/m */

        factor = K / (MASS_ALPHA * R * R * R);
    }
    else
    {
        /* Outside nucleus or point nucleus:

        F_x = K*x/r^3
        F_y = K*y/r^3

        a = F/m */

        if (r == 0.0) {
            throw std::runtime_error(
                "Point Coulomb singularity at r = 0"
            );
        }

        factor = K / (MASS_ALPHA * r2 * r);
    }

    double ax = factor * s.x;
    double ay = factor * s.y;

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

State add(
    const State& s,
    const Derivative& d,
    double factor
)
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

State rk4Step(
    const State& s,
    double dt,
    NuclearModel model
)
{
    // K1 = dt * f(Y_i)

    Derivative f1 = equations(s, model);

    Derivative K1 {
        dt * f1.dx,
        dt * f1.dy,
        dt * f1.dvx,
        dt * f1.dvy
    };


    // K2 = dt * f(Y_i + K1/2)

    State s2 = add(s, K1, 0.5);

    Derivative f2 = equations(s2, model);

    Derivative K2 {
        dt * f2.dx,
        dt * f2.dy,
        dt * f2.dvx,
        dt * f2.dvy
    };


    // K3 = dt * f(Y_i + K2/2)

    State s3 = add(s, K2, 0.5);

    Derivative f3 = equations(s3, model);

    Derivative K3 {
        dt * f3.dx,
        dt * f3.dy,
        dt * f3.dvx,
        dt * f3.dvy
    };


    // K4 = dt * f(Y_i + K3)

    State s4 = add(s, K3, 1.0);

    Derivative f4 = equations(s4, model);

    Derivative K4 {
        dt * f4.dx,
        dt * f4.dy,
        dt * f4.dvx,
        dt * f4.dvy
    };


    // Delta Y = (K1 + 2K2 + 2K3 + K4) / 6

    State result;

    result.x =
        s.x + (K1.dx + 2.0*K2.dx +
               2.0*K3.dx + K4.dx) / 6.0;

    result.y =
        s.y + (K1.dy + 2.0*K2.dy +
               2.0*K3.dy + K4.dy) / 6.0;

    result.vx =
        s.vx + (K1.dvx + 2.0*K2.dvx +
                2.0*K3.dvx + K4.dvx) / 6.0;

    result.vy =
        s.vy + (K1.dvy + 2.0*K2.dvy +
                2.0*K3.dvy + K4.dvy) / 6.0;

    return result;
}


/* ==========================================================
Total mechanical energy
============================================================*/

double totalEnergy(
    const State& s,
    NuclearModel model
)
{
    double r = std::hypot(s.x, s.y);

    double v2 =
        s.vx * s.vx + s.vy * s.vy;

    double kinetic =
        0.5 * MASS_ALPHA * v2;

    double potential =
        potentialEnergy(r, model);

    return kinetic + potential;
}


/* ==========================================================
Angular momentum L_z
============================================================*/

double angularMomentum(const State& s)
{
    return MASS_ALPHA *
           (s.x * s.vy - s.y * s.vx);
}


/* ==========================================================
Simulation result
============================================================*/

struct SimulationResult {
    double theta;
    double rMin;
    double energyError;
    double angularMomentumError;
    bool enteredNucleus;
};


/* ==========================================================
Scattering simulation
============================================================*/

SimulationResult simulateScattering(
    double E,
    double b,
    double x0,
    double dt,
    NuclearModel model
)
{
    // Initial position

    double initialX = -x0;
    double initialY = b;

    double rInitial =
        std::hypot(initialX, initialY);


    /* Initial velocity

    E = 1/2 m v^2 + U(r0)
    v = sqrt(2(E - U(r0))/m) */

    double potential0 =
        potentialEnergy(rInitial, model);

    double kinetic0 = E - potential0;

    if (kinetic0 <= 0.0) {
        throw std::runtime_error(
            "Initial kinetic energy must be positive"
        );
    }

    double initialSpeed =
        std::sqrt(2.0 * kinetic0 / MASS_ALPHA);

    State state {
        initialX,
        initialY,
        initialSpeed,
        0.0
    };


    // Initial energy and angular momentum

    double energyInitial =
        totalEnergy(state, model);

    double LInitial =
        angularMomentum(state);


    // Minimum distance from nucleus center

    double rMin = rInitial;

    bool finished = false;


    // Simulation loop

    const int MAX_STEPS = 2000000;

    for (int i = 0; i < MAX_STEPS; ++i)
    {
        double r =
            std::hypot(state.x, state.y);

        rMin = std::min(rMin, r);

        double radialVelocity =
            state.x * state.vx +
            state.y * state.vy;

        if (r >= rInitial && radialVelocity > 0.0)
        {
            finished = true;
            break;
        }

        // One RK4 step

        state = rk4Step(state, dt, model);
    }

    if (!finished) {
        throw std::runtime_error(
            "Simulation did not reach the exit boundary"
        );
    }


    // Final energy and angular momentum

    double energyFinal =
        totalEnergy(state, model);

    double LFinal =
        angularMomentum(state);


    // Relative energy error

    double energyError =
        std::abs(energyFinal - energyInitial)
        / std::abs(energyInitial);


    // Angular momentum error

    double angularMomentumError;

    if (std::abs(LInitial) > 1e-12)
    {
        angularMomentumError =
            std::abs(LFinal - LInitial)
            / std::abs(LInitial);
    }
    else
    {
        // For b = 0, use absolute error

        angularMomentumError =
            std::abs(LFinal - LInitial);
    }


    // Numerical scattering angle

    double theta =
        std::abs(std::atan2(state.vy, state.vx));


    // Return simulation results

    return {
        theta,
        rMin,
        energyError,
        angularMomentumError,
        rMin < R
    };
}


// MAIN


int main()
{
    // Simulation parameters

    const double x0 = 10000.0;  // fm
    const double dt = 1.0;      // fm/c


    // Alpha particle energies [MeV]

    std::vector<double> energies;

    for (double E = 5.0; E < 30.0; E += 5.0)
    {
        energies.push_back(E);
    }

    for (double E = 30.0; E <= 55.0; E += 1.0)
    {
        energies.push_back(E);
    }

    for (double E = 60.0; E <= 100.0; E += 5.0)
    {
        energies.push_back(E);
    }


    // Impact parameters [fm]

    std::vector<double> impactParameters;

    for (double b = 0.0; b < 10.0; b += 0.25)
    {
        impactParameters.push_back(b);
    }

    for (double b = 10.0; b < 30.0; b += 1.0)
    {
        impactParameters.push_back(b);
    }

    for (double b = 30.0; b <= 100.0; b += 5.0)
    {
        impactParameters.push_back(b);
    }


    // Output

    std::ofstream file(
        "data/finite_size_scattering.csv"
    );

    if (!file.is_open()) {
        std::cerr << "Cannot open output CSV file.\n";
        return 1;
    }

    file << std::setprecision(15);

    file
        << "E,b,"
        << "r_min_point,r_min_sphere,"
        << "theta_point,theta_sphere,"
        << "delta_theta,"
        << "entered_nucleus,"
        << "energy_error_point,"
        << "energy_error_sphere,"
        << "angular_momentum_error_point,"
        << "angular_momentum_error_sphere\n";


    // Print simulation parameters

    std::cout << std::setprecision(10);

    std::cout
        << "Finite-size nucleus scattering simulation\n\n";

    std::cout
        << "Gold nucleus radius R = "
        << R << " fm\n";

    std::cout
        << "Surface Coulomb energy K/R = "
        << K / R << " MeV\n\n";


    /* ========================================================
    ENERGY AND IMPACT PARAMETER SCAN
    ======================================================*/

    const int total =
        static_cast<int>(
            energies.size() * impactParameters.size()
        );

    int completed = 0;

    for (double E : energies)
    {
        std::cout
            << "\nEnergy E = " << E << " MeV\n";

        for (double b : impactParameters)
        {
            // Point nucleus simulation

            SimulationResult point =
                simulateScattering(
                    E, b, x0, dt,
                    NuclearModel::PointNucleus
                );


            // Uniform charged sphere simulation

            SimulationResult sphere =
                simulateScattering(
                    E, b, x0, dt,
                    NuclearModel::UniformSphere
                );


            // Convert radians to degrees

            double thetaPointDeg =
                point.theta * 180.0 / PI;

            double thetaSphereDeg =
                sphere.theta * 180.0 / PI;


            // Difference between nuclear models

            double deltaTheta =
                thetaSphereDeg - thetaPointDeg;


            // Save results

            file
                << E << ","
                << b << ","
                << point.rMin << ","
                << sphere.rMin << ","
                << thetaPointDeg << ","
                << thetaSphereDeg << ","
                << deltaTheta << ","
                << (sphere.enteredNucleus ? 1 : 0) << ","
                << point.energyError << ","
                << sphere.energyError << ","
                << point.angularMomentumError << ","
                << sphere.angularMomentumError
                << "\n";


            // Print progress every 100 simulations

            ++completed;

            if (completed % 100 == 0 || completed == total)
            {
                std::cout
                    << "Completed "
                    << completed
                    << " / "
                    << total
                    << " parameter combinations\n";
            }
        }
    }


    file.close();


    std::cout
        << "\nSimulation completed.\n";

    std::cout
        << "Results saved to "
        << "data/finite_size_scattering.csv\n";


    return 0;
}


