#ifndef TIMERSTEPPER_HPP
#define TIMERSTEPPER_HPP

#include <functional>
#include <exception>

#include "Newton.hpp"


namespace ASC_ode
{
  
  class TimeStepper
  { 
  protected:
    std::shared_ptr<NonlinearFunction> m_rhs;
  public:
    TimeStepper(std::shared_ptr<NonlinearFunction> rhs) : m_rhs(rhs) {}
    virtual ~TimeStepper() = default;
    virtual void DoStep(double tau, VectorView<double> y) = 0;
  };

  class ExplicitEuler : public TimeStepper
  {
    Vector<> m_vecf;
  public:
    ExplicitEuler(std::shared_ptr<NonlinearFunction> rhs) 
    : TimeStepper(rhs), m_vecf(rhs->dimF()) {}
    void DoStep(double tau, VectorView<double> y) override
    {
      this->m_rhs->evaluate(y, m_vecf);
      y += tau * m_vecf;
    }
  };

  class ImplicitEuler : public TimeStepper
  {
    std::shared_ptr<NonlinearFunction> m_equ;
    std::shared_ptr<Parameter> m_tau;
    std::shared_ptr<ConstantFunction> m_yold;
  public:
    ImplicitEuler(std::shared_ptr<NonlinearFunction> rhs) 
    : TimeStepper(rhs), m_tau(std::make_shared<Parameter>(0.0)) 
    {
      m_yold = std::make_shared<ConstantFunction>(rhs->dimX());
      auto ynew = std::make_shared<IdentityFunction>(rhs->dimX());
      m_equ = ynew - m_yold - m_tau * m_rhs;
    }

    void DoStep(double tau, VectorView<double> y) override
    {
      m_yold->set(y);
      m_tau->set(tau);
      NewtonSolver(m_equ, y);
    }
  };

class ImprovedEuler : public TimeStepper
{
    Vector<> m_vecf;
    Vector<> m_y_hat;

public:
    ImprovedEuler(std::shared_ptr<NonlinearFunction> rhs)
      : TimeStepper(rhs),
        m_vecf(rhs->dimF()),
        m_y_hat(rhs->dimX())
    {}

    void DoStep(double tau, VectorView<double> y) override
    {
        m_rhs->evaluate(y, m_vecf);          // k1 = f(y_n)
        m_y_hat = y;                         // y_hat = y_n
        m_y_hat += (tau/2.0) * m_vecf;       // y_hat += τ/2 * k1
        m_rhs->evaluate(m_y_hat, m_vecf);    // k2 = f(y_hat)
        y += tau * m_vecf;                   // y_{n+1} = y_n + τ * k2
    }
};

  

}


#endif
