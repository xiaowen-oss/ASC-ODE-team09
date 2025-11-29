#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

#include "vector.hpp"
#include "matrix.hpp"
#include "nonlinfunc.hpp"
#include "Newton.hpp"
#include "timestepper.hpp"
#include "implicitRK.hpp"
#include "explicitRK.hpp"

using namespace std;
using namespace ASC_ode;
using namespace nanoblas;

// mass-spring system: y = (position, velocity)
class MassSpringRHS : public NonlinearFunction
{
public:
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

// 2-stage explicit RK 
void SetupRK2(Matrix<> &A, Vector<> &b, Vector<> &c)
{
    A = 0.0;
    A(1,0) = 0.5;

    b(0) = 0.0;
    b(1) = 1.0;

    c(0) = 0.0;
    c(1) = 0.5;
}

// classical 4-stage RK4
void SetupRK4(Matrix<> &A, Vector<> &b, Vector<> &c)
{
    A = 0.0;
    A(1,0) = 0.5;
    A(2,1) = 0.5;
    A(3,2) = 1.0;

    b(0) = 1.0/6.0;
    b(1) = 1.0/3.0;
    b(2) = 1.0/3.0;
    b(3) = 1.0/6.0;

    c(0) = 0.0;
    c(1) = 0.5;
    c(2) = 0.5;
    c(3) = 1.0;
}

// run the method for several step sizes and write CSV
void SolveAndWrite(TimeStepper &stepper,
                   const vector<double> &taus,
                   double T,
                   const string &filename)
{
    ofstream fout(filename);
    fout << "tau,t,y,v\n";

    for(double tau : taus)
    {
        Vector<> u(2);
        u(0) = 1.0;   // initial position
        u(1) = 0.0;   // initial velocity

        double t = 0.0;
        fout << tau << "," << t << "," << u(0) << "," << u(1) << "\n";

        int steps = int(T / tau);
        for(int k = 0; k < steps; k++)
        {
            stepper.DoStep(tau, u);
            t += tau;
            fout << tau << "," << t << "," << u(0) << "," << u(1) << "\n";
        }
    }

    fout.close();
    cout << "Wrote " << filename << endl;
}

int main()
{
    vector<double> taus = {0.1, 0.05, 0.01};
    double T = 8.0 * M_PI;

    auto rhs = make_shared<MassSpringRHS>();

    // RK2 
    {
        int s = 2;
        Matrix<> A(s,s);
        Vector<> b(s), c(s);
        SetupRK2(A, b, c);

        ExplicitRungeKutta rk2(rhs, A, b, c);
        SolveAndWrite(rk2, taus, T, "ex19_4_rk2.csv");
    }

    // RK4 
    {
        int s = 4;
        Matrix<> A(s,s);
        Vector<> b(s), c(s);
        SetupRK4(A, b, c);

        ExplicitRungeKutta rk4(rhs, A, b, c);
        SolveAndWrite(rk4, taus, T, "ex19_4_rk4.csv");
    }



    // implicit RK with Gauss-Legendre 2-stage
    {
        ImplicitRungeKutta irk(rhs, Gauss2a, Gauss2b, Gauss2c);
        SolveAndWrite(irk, taus, T, "ex19_4_irk_gauss2.csv");
    }

    return 0;
}
