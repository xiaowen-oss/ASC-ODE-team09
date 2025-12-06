#ifndef MASS_SPRING_HPP
#define MASS_SPRING_HPP

#include <nonlinfunc.hpp>
#include <timestepper.hpp>
#include <vector>
#include <array>
#include <iostream>

using namespace ASC_ode;

#include <vector.hpp>

using namespace nanoblas;

// Helper operator to divide a vector by a scalar
template <size_t D>
Vec<D> operator/(const Vec<D> &v, double s)
{
    Vec<D> r;
    for (size_t i = 0; i < D; i++)
        r(i) = v(i) / s;
    return r;
}

// --- BASIC CLASSES ---

template <int D>
class Mass
{
public:
  double mass;
  Vec<D> pos;
  Vec<D> vel = Vec<D>({0, 0, -1});  // Default initial velocity (change as needed)

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

// Helper for printing connector info
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
  std::array<Connector, 2> connectors;
};

class DistanceConstraint
{
public:
    Connector c1;
    Connector c2;
    double rest_length;
    
    // Constructors for convenience
    DistanceConstraint() = default;
    DistanceConstraint(Connector _c1, Connector _c2, double _len) 
      : c1(_c1), c2(_c2), rest_length(_len) {}
};

// --- MASS-SPRING SYSTEM CLASS ---

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

  void addDistanceConstraint(const DistanceConstraint &dc)
  {
    m_constraints.push_back(dc);
  }

  auto & fixes() { return m_fixes; } 
  auto & masses() { return m_masses; } 
  auto & springs() { return m_springs; }
  auto & constraints() { return m_constraints; }

  // Get physical state (positions, velocities, accelerations)
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

  // Set physical state
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
  for (auto f : mss.fixes()) ost << f.pos << std::endl;
  ost << "masses: " << std::endl;
  for (auto m : mss.masses()) ost << "m = " << m.mass << ", pos = " << m.pos << std::endl;
  ost << "springs: " << std::endl;
  for (auto sp : mss.springs()) ost << "L=" << sp.length << ", k=" << sp.stiffness << std::endl;
  return ost;
}

// --- CLASS 1: SYSTEM MASS MATRIX (DAE) ---
// This part remains the same

// --- CLASS 2: FORCES AND CONSTRAINTS (Exact Derivative) ---
// This part remains the same
#endif
