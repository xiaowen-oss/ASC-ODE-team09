#include <iostream>
#include <fstream>
#include <memory>
#include <vector>
#include <cmath>

#include "vector.hpp"
#include "timestepper.hpp"
#include "nonlinfunc.hpp"
#include "RCCircuit.hpp"

using namespace ASC_ode;

int main()
{
    std::vector<double> taus = {0.004, 0.002, 0.001};
    double T = 0.2;   

    double R = 100.0;
    double C = 1e-6;

    std::shared_ptr<NonlinearFunction> rhs = std::make_shared<RCCircuit>(R, C);
    ExplicitEuler stepper(rhs);

    std::ofstream fout("xian.csv");
    fout << "tau,t,UC,timevar\n";

    for(double tau : taus)
    {
        int steps = T / tau;

        Vector<double> x(2);
        x(0) = 1.0;   // UC(0)
        x(1) = 0.0;   // t(0)

        fout << tau << "," << x(1) << "," << x(0) << "," << x(1) << "\n";

        for(int i=0; i<steps; i++)
        {
            stepper.DoStep(tau, x);

            fout << tau << "," << x(1) << "," << x(0) << "," << x(1) << "\n";
        }
    }

    fout.close();
    return 0;
}
