#pragma once
#include <cmath>
#include "nonlinfunc.hpp"


class RCRHS : public NonlinearFunction
{
    double m_R, m_C, m_omega;

public:
    RCRHS(double R, double C)
        : m_R(R), m_C(C), m_omega(100.0 * M_PI) {}

    size_t dimX() const override { return 2; }
    size_t dimF() const override { return 2; }

    void evaluate(VectorView<double> u, VectorView<double> f) const override
    {
        double UC = u(0);
        double t  = u(1);

        f(0) = (1.0 / (m_R * m_C)) * (std::cos(m_omega * t) - UC);
        f(1) = 1.0;
    }

    void evaluateDeriv(VectorView<double> u, MatrixView<double> df) const override
    {
        df = 0.0;

        double RC = m_R * m_C;
        double t = u(1);

        df(0,0) = -(1.0 / RC);
        df(0,1) = -(m_omega / RC) * std::sin(m_omega * t);

        df(1,0) = 0.0;
        df(1,1) = 0.0;
    }
};
