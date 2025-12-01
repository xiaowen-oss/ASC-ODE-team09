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

// ------------------ Implicit Euler Stepper (UC, t) --------------------
class ImplicitEulerRC
{
    std::shared_ptr<NonlinearFunction> m_rhs;      // f(u)
    std::shared_ptr<ConstantFunction> m_uold;      // u_old
    std::shared_ptr<Parameter> m_tau;              // tau
    std::shared_ptr<NonlinearFunction> m_equ;      // equation to solve

public:

    ImplicitEulerRC(double R, double C)
    {
        m_rhs  = std::make_shared<RCCircuit>(R, C);
        m_uold = std::make_shared<ConstantFunction>(2);
        m_tau  = std::make_shared<Parameter>(0.0);

        auto unew = std::make_shared<IdentityFunction>(2);

        // Implicit Euler:
        // u_new - u_old - tau * f(u_new) = 0
        m_equ = unew - m_uold - m_tau * m_rhs;
    }

    void DoStep(double tau, VectorView<double> u)
    {
        m_uold->set(u);
        m_tau->set(tau);

        // solve F(u_new) = 0
        NewtonSolver(m_equ, u, 1e-12, 20);
    }
};


// -------------------------- main ------------------------------

int main()
{
    // Use N steps instead of tau directly
    std::vector<int> Ns = {50, 100, 200, 300};
    double T = 0.2;

    double R = 100.0;
    double C = 1e-6;

    std::ofstream fout("yin.csv");
    fout << "N,t,UC,timevar\n";

    for(int N : Ns)
    {
        double tau = T / N;

        ImplicitEulerRC solver(R, C);

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
    std::cout << "yin.csv written.\n";
    return 0;
}
