#ifndef EXPLICITRK_HPP
#define EXPLICITRK_HPP

#include <vector.hpp>
#include <matrix.hpp>

#include "timestepper.hpp"
#include "nonlinfunc.hpp"

namespace ASC_ode
{
  using namespace nanoblas;

  // Explicit Runge–Kutta method 
  class ExplicitRungeKutta : public TimeStepper
  {
    Matrix<> m_a;
    Vector<> m_b, m_c;
    int m_stages;
    int m_n;

    Vector<> m_k;   // all stage derivatives
    Vector<> m_y;   // all stage states

  public:
    ExplicitRungeKutta(std::shared_ptr<NonlinearFunction> rhs,
                       const Matrix<> &a,
                       const Vector<> &b,
                       const Vector<> &c)
      : TimeStepper(rhs),
        m_a(a), m_b(b), m_c(c),
        m_stages(int(c.size())),
        m_n(int(rhs->dimX())),
        m_k(m_stages * m_n),
        m_y(m_stages * m_n)
    {
      if (m_a.rows() != m_stages || m_a.cols() != m_stages)
        throw std::runtime_error("ExplicitRungeKutta: A must be s x s");
      if (m_b.size() != size_t(m_stages))
        throw std::runtime_error("ExplicitRungeKutta: b must have size s");
    }

    void DoStep(double tau, VectorView<double> y) override
    {
      // compute stages
      for (int j = 0; j < m_stages; j++)
      {
        auto yj = m_y.range(j*m_n, (j+1)*m_n);
        yj = y;

        for (int ell = 0; ell < j; ell++)
          yj += tau * m_a(j, ell) * m_k.range(ell*m_n, (ell+1)*m_n);

        m_rhs->evaluate(yj, m_k.range(j*m_n, (j+1)*m_n));
      }

      // update solution
      for (int j = 0; j < m_stages; j++)
        y += tau * m_b(j) * m_k.range(j*m_n, (j+1)*m_n);
    }
  };

} // namespace ASC_ode

#endif // EXPLICITRK_HPP
