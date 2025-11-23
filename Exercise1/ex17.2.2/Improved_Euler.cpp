#include <iostream>
#include <fstream>
#include <memory>
#include <vector>
#include <cmath>

#include "vector.hpp"
#include "timestepper.hpp"
#include "nonlinfunc.hpp"

using namespace ASC_ode;

// def mass spring system
class MassSpring : public NonlinearFunction
{
    double m_mass;
    double m_stiffness;

public:
    MassSpring(double m, double k)
        : m_mass(m), m_stiffness(k) {}

    size_t dimX() const override { return 2; }
    size_t dimF() const override { return 2; }

    void evaluate(VectorView<double> x, VectorView<double> f) const override
    {
        f(0) = x(1);                         // dx/dt = v
        f(1) = -m_stiffness/m_mass * x(0);   // dv/dt = -k/m x
    }

    void evaluateDeriv(VectorView<double> x, MatrixView<double> df) const override
    {
        df = 0.0;
        df(0,1) = 1.0;
        df(1,0) = -m_stiffness/m_mass;
    }
};


// Improved Euler (Heun)

class ImprovedEuler : public TimeStepper
{
    Vector<double> m_vecf;
    Vector<double> m_vecf_tilde;
    Vector<double> m_ytilde;

public:
    ImprovedEuler(std::shared_ptr<NonlinearFunction> rhs)
        : TimeStepper(rhs),
          m_vecf(rhs->dimF()),
          m_vecf_tilde(rhs->dimF()),
          m_ytilde(rhs->dimF())
    {}

    void DoStep(double tau, VectorView<double> y) override
    {
        // Predictor stage
        m_rhs->evaluate(y, m_vecf);
        m_ytilde = y + 0.5 * tau * m_vecf;

        // Corrector stage
        m_rhs->evaluate(m_ytilde, m_vecf_tilde);
        y += tau * m_vecf_tilde;
    }
};

// main
int main()
{
    // time step sizes
    std::vector<double> taus = {0.1, 0.05, 0.01};

    // end time T = 8π
    double T = 8.0 * M_PI;

    // define ODE system
    std::shared_ptr<NonlinearFunction> rhs =
        std::make_shared<MassSpring>(1.0, 1.0);

    ImprovedEuler stepper(rhs);

    // CSV output
    std::ofstream fout("improved_euler.csv");
    fout << "tau,t,x,v\n";

    // run simulation for each tau
    for(double tau : taus)
    {
        int steps = int(T / tau);

        Vector<double> y(2);
        y(0) = 1.0;   // x(0)
        y(1) = 0.0;   // v(0)

        double t = 0.0;

        // write initial condition
        fout << tau << "," << t << "," << y(0) << "," << y(1) << "\n";

        for(int i = 0; i < steps; i++)
        {
            stepper.DoStep(tau, y);
            t += tau;

            fout << tau << "," << t << "," << y(0) << "," << y(1) << "\n";
        }
    }

    fout.close();
    std::cout << "Output written to improved_euler.csv\n";
    return 0;
}
