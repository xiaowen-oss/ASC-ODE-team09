#include <iostream>
#include <fstream>
#include "../src/autodiff.hpp"
#include "../src/nonlinfunc.hpp"

using namespace ASC_ode;

// ---------------------------------------------------------
// PendulumAD: f(x) for the pendulum using AutoDiff
// State x = [theta, omega]
// f(x) = [ omega,
//         -(g/L) * sin(theta) ]
// ---------------------------------------------------------
class PendulumAD : public NonlinearFunction
{
private:
    double m_length;
    double m_gravity;

public:
    PendulumAD(double length, double gravity = 9.81)
        : m_length(length), m_gravity(gravity) {}

    size_t dimX() const override { return 2; }
    size_t dimF() const override { return 2; }

    // -------------------------------------------------------
    // REQUIRED BY BASE CLASS:
    // f(x)
    // -------------------------------------------------------
    void evaluate(VectorView<double> x,
                  VectorView<double> f) const override
    {
        double theta = x(0);
        double omega = x(1);

        f(0) = omega;
        f(1) = -(m_gravity / m_length) * std::sin(theta);
    }

    // -------------------------------------------------------
    // df/dx using AutoDiff
    // -------------------------------------------------------
    void evaluateDeriv(VectorView<double> x,
                       MatrixView<double> df) const override
    {
        using AD = AutoDiff<2,double>;

        AD th = Variable<0,double>( x(0) );
        AD om = Variable<1,double>( x(1) );

        AD f0 = om;

        AD C(-(m_gravity/m_length));
        AD f1 = C * sin(th);

        df(0,0) = f0.deriv()[0];
        df(0,1) = f0.deriv()[1];
        df(1,0) = f1.deriv()[0];
        df(1,1) = f1.deriv()[1];
    }
};



// ---------------------------------------------------------
// Simple Explicit Euler to simulate the pendulum
// ---------------------------------------------------------
int main()
{
    PendulumAD pend(1.0);      // L = 1 meter
    double dt = 0.01;
    double T = 10.0;

    double theta0 = 0.5;        // initial angle
    double omega0 = 0.0;        // initial velocity

    double t = 0.0;
    double x0[2] = {theta0, omega0};
    double x1[2];

    std::ofstream fout("ex18_5.csv");
    fout << "t,theta,omega\n";

    Vector<double> xv(2), fv(2);

    while (t <= T)
    {
        xv(0) = x0[0];
        xv(1) = x0[1];

        pend.evaluate(xv, fv);

        // Explicit Euler
        x1[0] = x0[0] + dt * fv(0);
        x1[1] = x0[1] + dt * fv(1);

        // write CSV
        fout << t << "," << x0[0] << "," << x0[1] << "\n";

        // update
        x0[0] = x1[0];
        x0[1] = x1[1];
        t += dt;
    }

    fout.close();
    std::cout << "Simulation complete. CSV written: ex18_5.csv\n";

    return 0;
}
