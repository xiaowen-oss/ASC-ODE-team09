#include <iostream>
#include <fstream>
#include <vector>
#include <memory>
#include <cmath>

#include "nonlinfunc.hpp"
#include "Newton.hpp"
#include "vector.hpp"
#include "matrix.hpp"
#include "RCCircuit.hpp"

using namespace ASC_ode;

// ---------------- Crank–Nicolson stepper --------------------

class CrankNicolsonRC
{
    std::shared_ptr<NonlinearFunction> m_rhs;      // f(u)
    std::shared_ptr<ConstantFunction> m_uold;      // u_old
    std::shared_ptr<Parameter> m_tau_half;         // tau/2
    std::shared_ptr<NonlinearFunction> m_equ;      // equation for Newton

public:

    CrankNicolsonRC(double R, double C)
    {
        m_rhs = std::make_shared<RCCircuit>(R, C);
        m_uold = std::make_shared<ConstantFunction>(2);
        m_tau_half = std::make_shared<Parameter>(0.0);

        auto unew = std::make_shared<IdentityFunction>(2);

        // f(u_old)
        auto f_old = std::make_shared<ComposeFunction>(m_rhs, m_uold);

        // CN equation:
        //  u_new - u_old - τ/2 f(u_new) - τ/2 f(u_old) = 0
        m_equ = unew - m_uold - m_tau_half * m_rhs - m_tau_half * f_old;
    }

    void DoStep(double tau, VectorView<double> u)
    {
        m_uold->set(u);
        m_tau_half->set(0.5 * tau);

        NewtonSolver(m_equ, u, 1e-12, 20);
    }
};


// --------------------- main --------------------------

int main()
{
    std::vector<int> Ns = {50, 100, 200, 300};
    double T = 0.2;

    double R = 100.0;
    double C = 1e-6;

    std::ofstream fout("cn.csv");
    fout << "N,t,UC,timevar\n";

    for(int N : Ns)
    {
        double tau = T / N;

        CrankNicolsonRC solver(R, C);

        Vector<> u(2);
        u(0) = 1.0;    // UC(0)
        u(1) = 0.0;    // t(0)

        fout << N << "," << u(1) << "," << u(0) << "," << u(1) << "\n";

        for(int i = 0; i < N; i++)
        {
            solver.DoStep(tau, u);

            fout << N << "," << u(1) << "," << u(0) << "," << u(1) << "\n";
        }
    }

    fout.close();
    std::cout << "cn.\n";
    return 0;
}
