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
    // n = number of steps
    std::vector<int> Ns = {50, 100, 200, 300};
    double T = 0.2;

    double R = 100.0;
    double C = 1e-6;

    std::shared_ptr<NonlinearFunction> rhs = std::make_shared<RCCircuit>(R, C);
    ExplicitEuler stepper(rhs);

    std::ofstream fout("xian.csv");
    fout << "N,t,UC,timevar\n";

    for(int N : Ns)
    {
        double tau = T / N;   //  compute tau automatically

        Vector<double> x(2);
        x(0) = 1.0;   // UC(0)
        x(1) = 0.0;   // t(0) — time is state variable

        // output initial condition
        fout << N << "," << x(1) << "," << x(0) << "," << x(1) << "\n";

        for(int i = 0; i < N; i++)
        {
            stepper.DoStep(tau, x);   //  update both UC and t

            fout << N << "," << x(1) << "," << x(0) << "," << x(1) << "\n";
        }
    }

    fout.close();
    std::cout << "xian.\n";
    return 0;
}
