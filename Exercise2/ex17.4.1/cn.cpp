#include <iostream>
#include <fstream>
#include <memory>
#include <cmath>

#include "nonlinfunc.hpp"
#include "Newton.hpp"
#include "vector.hpp"
#include "matrix.hpp"

#include "RCRHS.hpp"

using namespace ASC_ode;


// Crank–Nicolson time stepper for RC circuit
class CrankNicolsonRC
{
    std::shared_ptr<NonlinearFunction> m_rhs;

    std::shared_ptr<NonlinearFunction> m_equ;
    std::shared_ptr<ConstantFunction> m_uold;
    std::shared_ptr<Parameter> m_tau_half;

public:
    CrankNicolsonRC(double R, double C)
    {
        m_rhs     = std::make_shared<RCRHS>(R, C);
        m_uold    = std::make_shared<ConstantFunction>(2);
        m_tau_half = std::make_shared<Parameter>(0.0);

        auto unew = std::make_shared<IdentityFunction>(2);

        // Crank–Nicolson equation:
        // u_new - u_old - (tau/2)*f(u_new) - (tau/2)*f(u_old) = 0
        // f(u_old) = rhs(u_old)
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



// --------------------------- main ----------------------------------

int main()
{
    std::vector<double> taus = {0.004, 0.002, 0.001};
    double T = 0.2;     // RC time scale small

    double R = 100.0;
    double C = 1e-6;

    std::ofstream fout("cn.csv");
    fout << "tau,t,UC,timevar\n";

    for(double tau : taus)
    {
        CrankNicolsonRC solver(R, C);

        Vector<> u = {1.0, 0.0};   // UC(0) = 1, t(0)=0
        double t = 0.0;

        fout << tau << "," << t << "," << u(0) << "," << u(1) << "\n";

        int steps = int(T / tau);
        for(int i = 0; i < steps; i++)
        {
            solver.DoStep(tau, u);
            t += tau;

            fout << tau << "," << t << "," << u(0) << "," << u(1) << "\n";
        }
    }

    fout.close();
    std::cout << "Output written to cn.csv\n";
    return 0;
}
