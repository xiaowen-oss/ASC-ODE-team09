#include <iostream>
#include <fstream>
#include <memory>

#include "nonlinfunc.hpp"
#include "Newton.hpp"
#include "vector.hpp"
#include "matrix.hpp"

using namespace ASC_ode;

// Mass–spring system
class MassSpringRHS : public NonlinearFunction
{
public:
    size_t dimX() const override { return 2; }
    size_t dimF() const override { return 2; }

    void evaluate(VectorView<double> u, VectorView<double> f) const override
    {
        f(0) = u(1);
        f(1) = -u(0);
    }

    void evaluateDeriv(VectorView<double> u, MatrixView<double> df) const override
    {
        df = 0.0;
        df(0,1) = 1.0;
        df(1,0) = -1.0;
    }
};


// Crank–Nicolson time stepper

class CrankNicolsonMS
{
    std::shared_ptr<NonlinearFunction> m_rhs;

    std::shared_ptr<NonlinearFunction> m_equ;
    std::shared_ptr<ConstantFunction> m_uold;
    std::shared_ptr<Parameter> m_tau_half;

public:
    CrankNicolsonMS()
    {
        m_rhs = std::make_shared<MassSpringRHS>();
        m_uold = std::make_shared<ConstantFunction>(2);

        m_tau_half = std::make_shared<Parameter>(0.0);

        auto unew = std::make_shared<IdentityFunction>(2);

        // Crank–Nicolson:
        // u_new - u_old - (τ/2) f(u_new) - (τ/2) f(u_old) = 0
        auto f_old = std::make_shared<ComposeFunction>(m_rhs, m_uold);

        m_equ = unew - m_uold - (m_tau_half * m_rhs) - (m_tau_half * f_old);
    }

    void DoStep(double tau, VectorView<double> u)
    {
        m_uold->set(u);
        m_tau_half->set(0.5 * tau);

        NewtonSolver(m_equ, u, 1e-12, 20);
    }
};


// main
int main()
{
    std::vector<double> taus = {0.1, 0.05, 0.01};
    double T = 8 * M_PI;

    std::ofstream fout("crank_nicolson.csv");
    fout << "tau,t,y,v\n";

    for(double tau : taus)
    {
        CrankNicolsonMS stepper;

        Vector<> u = {1.0, 0.0};
        double t = 0.0;

        fout << tau << "," << t << "," << u(0) << "," << u(1) << "\n";

        int steps = int(T / tau);
        for(int i = 0; i < steps; i++)
        {
            stepper.DoStep(tau, u);
            t += tau;

            fout << tau << "," << t << "," << u(0) << "," << u(1) << "\n";
        }
    }

    fout.close();
    std::cout << "Output written to crank_nicolson.csv\n";
    return 0;
}
