#pragma once

#include <cmath>
#include "nonlinfunc.hpp"

using namespace ASC_ode;

class RCCircuit : public NonlinearFunction
{
    double m_R, m_C;

public:
    RCCircuit(double R, double C)
        : m_R(R), m_C(C) {}

    size_t dimX() const override { return 2; }
    size_t dimF() const override { return 2; }

    // f(x)
    void evaluate(VectorView<double> x, VectorView<double> f) const override
    {
        double UC = x(0);   // capacitor voltage
        double t  = x(1);   // time state variable
        double omega = 100.0 * M_PI;

        f(0) = (1.0 / (m_R * m_C)) * ( std::cos(omega * t) - UC );
        f(1) = 1.0;   // dx2/dt = 1
    }

    // df/dx (Jacobian) — for later implicit Euler & CN
    void evaluateDeriv(VectorView<double> x, MatrixView<double> df) const override
    {
        df = 0.0;

        double RC = m_R * m_C;
        double t  = x(1);
        double omega = 100.0 * M_PI;

        df(0,0) = -(1.0 / RC);
        df(0,1) = -(omega / RC) * std::sin(omega * t);

        df(1,0) = 0.0;
        df(1,1) = 0.0;
    }
};
