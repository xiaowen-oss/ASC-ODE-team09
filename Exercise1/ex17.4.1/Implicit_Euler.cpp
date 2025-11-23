#include <iostream>
#include <fstream>
#include <cmath>
#include "nonlinfunc.hpp"
#include "Newton.hpp"
#include "vector.hpp"
#include "matrix.hpp"

using namespace std;
using namespace ASC_ode;

// system
class MassSpringRHS : public NonlinearFunction
{
public:
    MassSpringRHS() {}

    size_t dimX() const override { return 2; }
    size_t dimF() const override { return 2; }

    void evaluate(VectorView<double> u, VectorView<double> f) const override
    {
        double y = u(0);
        double v = u(1);

        f(0) = v;
        f(1) = -y;
    }

    void evaluateDeriv(VectorView<double> u, MatrixView<double> df) const override
    {
        df = 0.0;
        df(0,1) = 1.0;
        df(1,0) = -1.0;
    }
};


// Implicit Euler time-stepper:  u - u_old - tau f(u) = 0

class ImplicitEulerMS
{
    shared_ptr<NonlinearFunction> m_rhs;
    shared_ptr<NonlinearFunction> m_equ;
    shared_ptr<Parameter> m_tau;
    shared_ptr<ConstantFunction> m_uold;

public:
    ImplicitEulerMS()
    {
        m_rhs  = make_shared<MassSpringRHS>();
        m_tau  = make_shared<Parameter>(0.0);
        m_uold = make_shared<ConstantFunction>(2);

        auto unew = make_shared<IdentityFunction>(2);

        // u - u_old - tau * f(u)
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
    vector<double> taus = {0.1, 0.05, 0.01};
    double T = 8.0 * M_PI;

    ofstream fout("implicit_euler.csv");
    fout << "tau,t,y,v\n";

    for(double tau : taus)
    {
        ImplicitEulerMS solver;

        Vector<> u(2);
        u(0) = 1.0;   // initial position
        u(1) = 0.0;   // initial velocity

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
    cout << "Output written to implicit_euler.csv\n";

    return 0;
}
