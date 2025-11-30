#include <iostream>
#include <fstream>
#include <vector>
#include <memory>
#include <cmath>

#include "nonlinfunc.hpp"
#include "Newton.hpp"
#include "vector.hpp"
#include "matrix.hpp"
#include "RCRHS.hpp"

using namespace std;
using namespace ASC_ode;

// implicit Euler stepper for RC circuit
class ImplicitEulerRC
{
    shared_ptr<NonlinearFunction> m_rhs;
    shared_ptr<NonlinearFunction> m_equ;
    shared_ptr<Parameter> m_tau;
    shared_ptr<ConstantFunction> m_uold;

public:
    ImplicitEulerRC(double R, double C)
    {
        m_rhs  = make_shared<RCRHS>(R, C);
        m_tau  = make_shared<Parameter>(0.0);
        m_uold = make_shared<ConstantFunction>(2);

        auto unew = make_shared<IdentityFunction>(2);

        // equation: u - u_old - tau * f(u) = 0
        m_equ = unew - m_uold - m_tau * m_rhs;
    }

    void DoStep(double tau, VectorView<double> u)
    {
        m_uold->set(u);
        m_tau->set(tau);

        NewtonSolver(m_equ, u, 1e-12, 20);
    }
};


// main
int main()
{
    std::vector<double> taus = {0.004, 0.002, 0.001};
    double T = 0.2;   


    double R = 100.0;
    double C = 1e-6;

    ofstream fout("yin.csv");
    fout << "tau,t,UC,timevar\n";

    for(double tau : taus)
    {
        ImplicitEulerRC solver(R, C);

        // state u = (UC, t)
        Vector<> u(2);
        u(0) = 1.0;  // UC(0) = cos(0)
        u(1) = 0.0;  // t(0)

        double t = 0.0;

        fout << tau << "," << t << "," << u(0) << "," << u(1) << "\n";

        int steps = int(T / tau);

        for(int k = 0; k < steps; k++)
        {
            solver.DoStep(tau, u);
            t += tau;

            fout << tau << "," << t << "," << u(0) << "," << u(1) << "\n";
        }
    }

    fout.close();
    cout << "Output written to yin.csv\n";
    return 0;
}
