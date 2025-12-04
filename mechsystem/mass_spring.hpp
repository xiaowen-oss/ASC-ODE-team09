#ifndef MASS_SPRING_HPP
#define MASS_SPRING_HPP

#include <nonlinfunc.hpp>
#include <timestepper.hpp>

using namespace ASC_ode;

#include <vector.hpp>

using namespace nanoblas;

// operator a/b
template <size_t D>
Vec<D> operator/(const Vec<D> &v, double s)
{
    Vec<D> r;
    for (size_t i = 0; i < D; i++)
        r(i) = v(i) / s;
    return r;
}



template <int D>
class Mass
{
public:
  double mass;
  Vec<D> pos;
  Vec<D> vel = 0.0;
  Vec<D> acc = 0.0;
};


template <int D>
class Fix
{
public:
  Vec<D> pos;
};


class Connector
{
public:
  enum CONTYPE { FIX=1, MASS=2 };
  CONTYPE type;
  size_t nr;
};

std::ostream & operator<< (std::ostream & ost, const Connector & con)
{
  ost << "type = " << int(con.type) << ", nr = " << con.nr;
  return ost;
}

class Spring
{
public:
  double length;  
  double stiffness;
  std::array<Connector,2> connectors;
};


class DistanceConstraint
{
public:
    Connector c1;
    Connector c2;
    double rest_length;
};


template <int D>
class MassSpringSystem
{
  std::vector<Fix<D>> m_fixes;
  std::vector<Mass<D>> m_masses;
  std::vector<Spring> m_springs;
  std::vector<DistanceConstraint> m_constraints;
  Vec<D> m_gravity=0.0;
public:
  void setGravity (Vec<D> gravity) { m_gravity = gravity; }
  Vec<D> getGravity() const { return m_gravity; }

  Connector addFix (Fix<D> p)
  {
    m_fixes.push_back(p);
    return { Connector::FIX, m_fixes.size()-1 };
  }

  Connector addMass (Mass<D> m)
  {
    m_masses.push_back (m);
    return { Connector::MASS, m_masses.size()-1 };
  }
  
  size_t addSpring (Spring s) 
  {
    m_springs.push_back (s); 
    return m_springs.size()-1;
  }

  auto & fixes() { return m_fixes; } 
  auto & masses() { return m_masses; } 
  auto & springs() { return m_springs; }
  auto & constraints() { return m_constraints;}

  void addDistanceConstraint(const DistanceConstraint &dc)
  {
    m_constraints.push_back(dc);
  }

  void getState (VectorView<> values, VectorView<> dvalues, VectorView<> ddvalues)
  {
    auto valmat = values.asMatrix(m_masses.size(), D);
    auto dvalmat = dvalues.asMatrix(m_masses.size(), D);
    auto ddvalmat = ddvalues.asMatrix(m_masses.size(), D);

    for (size_t i = 0; i < m_masses.size(); i++)
      {
        valmat.row(i) = m_masses[i].pos;
        dvalmat.row(i) = m_masses[i].vel;
        ddvalmat.row(i) = m_masses[i].acc;
      }
  }

  void setState (VectorView<> values, VectorView<> dvalues, VectorView<> ddvalues)
  {
    auto valmat = values.asMatrix(m_masses.size(), D);
    auto dvalmat = dvalues.asMatrix(m_masses.size(), D);
    auto ddvalmat = ddvalues.asMatrix(m_masses.size(), D);

    for (size_t i = 0; i < m_masses.size(); i++)
      {
        m_masses[i].pos = valmat.row(i);
        m_masses[i].vel = dvalmat.row(i);
        m_masses[i].acc = ddvalmat.row(i);
      }
  }
};

template <int D>
std::ostream & operator<< (std::ostream & ost, MassSpringSystem<D> & mss)
{
  ost << "fixes:" << std::endl;
  for (auto f : mss.fixes())
    ost << f.pos << std::endl;

  ost << "masses: " << std::endl;
  for (auto m : mss.masses())
    ost << "m = " << m.mass << ", pos = " << m.pos << std::endl;

  ost << "springs: " << std::endl;
  for (auto sp : mss.springs())
    ost << "length = " << sp.length << ", stiffness = " << sp.stiffness
        << ", C1 = " << sp.connectors[0] << ", C2 = " << sp.connectors[1] << std::endl;
  return ost;
}


template <int D>
class MSS_Function : public NonlinearFunction
{
  MassSpringSystem<D> & mss;
public:
  MSS_Function (MassSpringSystem<D> & _mss)
    : mss(_mss) { }

  virtual size_t dimX() const override { return D*mss.masses().size(); }
  virtual size_t dimF() const override{ return D*mss.masses().size(); }

  virtual void evaluate (VectorView<double> x, VectorView<double> f) const override
  {
    f = 0.0;

    auto xmat = x.asMatrix(mss.masses().size(), D);
    auto fmat = f.asMatrix(mss.masses().size(), D);

    for (size_t i = 0; i < mss.masses().size(); i++)
      fmat.row(i) = mss.masses()[i].mass*mss.getGravity();

    for (auto spring : mss.springs())
      {
        auto [c1,c2] = spring.connectors;
        Vec<D> p1, p2;
        if (c1.type == Connector::FIX)
          p1 = mss.fixes()[c1.nr].pos;
        else
          p1 = xmat.row(c1.nr);
        if (c2.type == Connector::FIX)
          p2 = mss.fixes()[c2.nr].pos;
        else
          p2 = xmat.row(c2.nr);

        double force = spring.stiffness * (norm(p1-p2)-spring.length);
        Vec<D> dir12 = 1.0/norm(p1-p2) * (p2-p1);
        if (c1.type == Connector::MASS)
          fmat.row(c1.nr) += force*dir12;
        if (c2.type == Connector::MASS)
          fmat.row(c2.nr) -= force*dir12;
      }

    //distance constraints
    double kp = 10000.0;

    for (auto &dc : mss.constraints())
    {
      auto c1 = dc.c1;
      auto c2 = dc.c2;

        // positions
        Vec<D> p1 = (c1.type == Connector::FIX ? mss.fixes()[c1.nr].pos : xmat.row(c1.nr));
        Vec<D> p2 = (c2.type == Connector::FIX ? mss.fixes()[c2.nr].pos : xmat.row(c2.nr));

        Vec<D> d = p1 - p2;
        double L = norm(d);
        if (L < 1e-12) continue;

        Vec<D> dir = d / L;
        double phi = L - dc.rest_length;   // constraint violation

        Vec<D> f12 = -kp * phi * dir;      // penalty force

        if (c1.type == Connector::MASS)
            fmat.row(c1.nr) += f12;

        if (c2.type == Connector::MASS)
            fmat.row(c2.nr) -= f12;
    }



    for (size_t i = 0; i < mss.masses().size(); i++)
      fmat.row(i) *= 1.0/mss.masses()[i].mass;
  }
  
  virtual void evaluateDeriv (VectorView<double> x, MatrixView<double> df) const override
  {
    // TODO: exact differentiation
    df = 0.0;

    size_t N = mss.masses().size();
    auto X = x.asMatrix(N, D);

    // Loop over all springs
    for (auto &spring : mss.springs())
    {
        auto [c1, c2] = spring.connectors;

        //get positions
        Vec<D> p1 = (c1.type == Connector::FIX ? mss.fixes()[c1.nr].pos : X.row(c1.nr));
        Vec<D> p2 = (c2.type == Connector::FIX ? mss.fixes()[c2.nr].pos : X.row(c2.nr));

        Vec<D> d = p2 - p1;
        double L = norm(d);

        if (L < 1e-12) continue; // avoid 0 division

        Vec<D> n = d / L;               // unit direction
        double k = spring.stiffness;
        double L0 = spring.length;
        double fmag = k * (L - L0);

        // For Jacobian:  df/dp = k * ( n*n^T + ((L-L0)/L)*(I - n*n^T) )
        double a = fmag / L;            // (L-L0)/L * k

        // Build 3x3 or 2x2 block
        for (size_t i = 0; i < D; i++)
        for (size_t j = 0; j < D; j++)
        {
            double dij = n(i)*n(j);                // n n^T
            double Kij = k * dij + a * ((i==j?1.0:0.0) - dij);

            // F1 = +K * dX
            // F2 = -K * dX

            if (c1.type == Connector::MASS)
            {
                df(c1.nr*D + i, c1.nr*D + j) -= Kij;
            }

            if (c2.type == Connector::MASS)
            {
                df(c2.nr*D + i, c2.nr*D + j) -= Kij;
            }

            // off-diagonal blocks
            if (c1.type == Connector::MASS && c2.type == Connector::MASS)
            {
                df(c1.nr*D + i, c2.nr*D + j) += Kij;
                df(c2.nr*D + i, c1.nr*D + j) += Kij;
            }
        }
      }
  }
  
};

#endif
