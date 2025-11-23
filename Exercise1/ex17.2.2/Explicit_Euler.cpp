#include <iostream>
#include <fstream>
#include <memory>
#include <vector>

#include "vector.hpp"
#include "timestepper.hpp"
#include "nonlinfunc.hpp"

using namespace ASC_ode;

//define mass-spring system
class MassSpring : public NonlinearFunction
{
    double m_mass, m_stiffness;

public:
    MassSpring(double m, double k)
        : m_mass(m), m_stiffness(k) {}

    size_t dimX()  const override { return 2; }
    size_t dimF()  const override { return 2; }

    void evaluate(VectorView<double> x, VectorView<double> f) const override
    {
        f(0) = x(1);
        f(1) = -m_stiffness/m_mass * x(0);
    }

    void evaluateDeriv(VectorView<double> x, MatrixView<double> df) const override
    {
        df = 0.0;
        df(0,1) = 1.0;
        df(1,0) = -m_stiffness/m_mass;
    }
};

//use explicit method with different time steps
int main()
{
    std::vector<double> taus = {0.1, 0.05, 0.01};
    double T = 8*M_PI;

    std::shared_ptr<NonlinearFunction> rhs = std::make_shared<MassSpring>(1.0,1.0);

    ExplicitEuler stepper(rhs);   

    std::ofstream fout("explicit_euler.csv");
    fout << "tau,t,x,v\n";

    for(double tau : taus)
    {
        int steps = T / tau;

        Vector<double> y(2);
        y(0) = 1.0;
        y(1) = 0.0;

        double t = 0.0;
        fout << tau << "," << t << "," << y(0) << "," << y(1) << "\n";

        for(int i=0;i<steps;i++)
        {
            stepper.DoStep(tau, y);   
            t += tau;

            fout << tau << "," << t << "," << y(0) << "," << y(1) << "\n";
        }
    }

    fout.close();
    std::cout << "Output written to explicit_euler.csv\n";
    return 0;
}
