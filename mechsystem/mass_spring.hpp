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
  std::array<Connector,2> connectors;
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
  auto & constraints() { return m_constraints;}

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
// Returns the mass matrix for the augmented system:
// M = diag(m1, m1, ..., 0, 0)
// The zeros on the diagonal correspond to the Lagrange multipliers (constraints)

template <int D>
class SystemMassFunction : public NonlinearFunction
{
  MassSpringSystem<D> & mss;
public:
  SystemMassFunction (MassSpringSystem<D> & _mss) : mss(_mss) { }

  virtual size_t dimX() const override { 
      return D * mss.masses().size() + mss.constraints().size(); 
  }
  virtual size_t dimF() const override { 
      return D * mss.masses().size() + mss.constraints().size(); 
  }

  // f = M * x (where 'x' is the acceleration vector in the solver context)
  virtual void evaluate (VectorView<double> x, VectorView<double> f) const override
  {
    f = 0.0;
    // Mass rows (physical DOFs)
    for (size_t i = 0; i < mss.masses().size(); i++)
      for (size_t d = 0; d < D; d++)
        f(i*D + d) = mss.masses()[i].mass * x(i*D + d);
    
    // Constraint rows (Lagrange multipliers): 0 * lambda = 0
  }

  // Returns the constant mass matrix (Deriv)
  virtual void evaluateDeriv (VectorView<double> x, MatrixView<double> df) const override
  {
    df = 0.0;
    for (size_t i = 0; i < mss.masses().size(); i++)
      for (size_t d = 0; d < D; d++)
        df(i*D + d, i*D + d) = mss.masses()[i].mass;
  }
};


// --- CLASS 2: FORCES AND CONSTRAINTS (Exact Derivative) ---
// This function evaluates the RHS of the equation: F(x, lambda)
// F = [ Forces_physical + Gradient^T * lambda ]
//     [ Constraint_Equation C(x)              ]

template <int D>
class MSS_Function : public NonlinearFunction
{
  MassSpringSystem<D> & mss;
public:
  MSS_Function (MassSpringSystem<D> & _mss)
    : mss(_mss) { }

  virtual size_t dimX() const override { return D*mss.masses().size() + mss.constraints().size(); }
  virtual size_t dimF() const override{ return D*mss.masses().size() + mss.constraints().size(); }

  virtual void evaluate (VectorView<double> x, VectorView<double> f) const override
  {
    f = 0.0;
    size_t n_masses = mss.masses().size();
    size_t n_constraints = mss.constraints().size();

    auto xmat = x.asMatrix(n_masses, D);
    auto fmat = f.asMatrix(n_masses, D); 

    // 1. Gravity (External Force)
    for (size_t i = 0; i < n_masses; i++)
      fmat.row(i) = mss.masses()[i].mass*mss.getGravity();

    // 2. Springs (Elastic Force)
    for (auto spring : mss.springs())
      {
        auto [c1,c2] = spring.connectors;
        Vec<D> p1, p2;
        if (c1.type == Connector::FIX) p1 = mss.fixes()[c1.nr].pos;
        else p1 = xmat.row(c1.nr);
        
        if (c2.type == Connector::FIX) p2 = mss.fixes()[c2.nr].pos;
        else p2 = xmat.row(c2.nr);

        double L = norm(p1-p2);
        if (L < 1e-12) continue;

        double force = spring.stiffness * (L-spring.length);
        Vec<D> dir12 = (1.0/L) * (p2-p1); // Direction p1 -> p2
        
        if (c1.type == Connector::MASS) fmat.row(c1.nr) += force*dir12;
        if (c2.type == Connector::MASS) fmat.row(c2.nr) -= force*dir12;
      }

    // 3. Constraints (Lagrange Forces and Constraint Equations)
    for (size_t i = 0; i < n_constraints; i++)
    {
        auto& dc = mss.constraints()[i];
        double lambda = x(D * n_masses + i); // The multiplier is at the end of vector x

        Vec<D> p1 = (dc.c1.type == Connector::FIX) ? mss.fixes()[dc.c1.nr].pos : xmat.row(dc.c1.nr);
        Vec<D> p2 = (dc.c2.type == Connector::FIX) ? mss.fixes()[dc.c2.nr].pos : xmat.row(dc.c2.nr);
        
        Vec<D> d = p1 - p2; // Vector p2 -> p1
        double L = norm(d);
        if (L < 1e-12) continue;
        
        Vec<D> dir = d / L; // Gradient of C with respect to p1

        // Constraint Force (-lambda * gradient) moved to the RHS (force side)
        // Since we solve Ma = F, and Lagrange term is usually on LHS (Ma + G^T lambda = F_ext),
        // we move it to RHS: Ma = F_ext - G^T lambda.
        if (dc.c1.type == Connector::MASS) fmat.row(dc.c1.nr) -= lambda * dir;
        if (dc.c2.type == Connector::MASS) fmat.row(dc.c2.nr) += lambda * dir;

        // Constraint Equation: C(x) = L - L0 = 0
        // Since the solver solves M*a - F = 0, and for the lambda row M=0, we need -F_lambda = 0.
        // We set F_lambda = L - L0. So 0 - (L - L0) = 0  => L = L0.
        f(D * n_masses + i) = (L - dc.rest_length);
    }
  }
  
  // Exact Derivative (Jacobian Matrix)
  virtual void evaluateDeriv (VectorView<double> x, MatrixView<double> df) const override
  {
    df = 0.0;
    size_t n_masses = mss.masses().size();
    auto X = x.asMatrix(n_masses, D);

    // --- Part A: Spring Stiffness ---
    for (auto &spring : mss.springs())
    {
        auto [c1, c2] = spring.connectors;
        Vec<D> p1 = (c1.type == Connector::FIX) ? mss.fixes()[c1.nr].pos : X.row(c1.nr);
        Vec<D> p2 = (c2.type == Connector::FIX) ? mss.fixes()[c2.nr].pos : X.row(c2.nr);

        Vec<D> d = p2 - p1;
        double L = norm(d);
        if (L < 1e-12) continue; 

        Vec<D> n = d / L;
        double k = spring.stiffness;
        double L0 = spring.length;
        
        // Geometric stiffness term due to spring tension
        double force_over_L = k * (L - L0) / L;

        for (size_t i = 0; i < D; i++)
        for (size_t j = 0; j < D; j++)
        {
            double ninj = n(i)*n(j);
            // K_ij = k * n_i*n_j + (f/L) * (delta_ij - n_i*n_j)
            double Kij = k * ninj + force_over_L * ((i==j?1.0:0.0) - ninj);
            
            // We want dF/dx. Spring force pulls towards the other point.
            if (c1.type == Connector::MASS) 
                df(c1.nr*D + i, c1.nr*D + j) -= Kij; 
            
            if (c2.type == Connector::MASS) 
                df(c2.nr*D + i, c2.nr*D + j) -= Kij; 
            
            if (c1.type == Connector::MASS && c2.type == Connector::MASS) {
                df(c1.nr*D + i, c2.nr*D + j) += Kij; 
                df(c2.nr*D + i, c1.nr*D + j) += Kij; 
            }
        }
    }

    // --- Part B: Constraints and Multipliers ---
    for (size_t k = 0; k < mss.constraints().size(); k++)
    {
        auto &dc = mss.constraints()[k];
        double lambda = x(D * n_masses + k);
        size_t idx_lambda = D * n_masses + k;

        Vec<D> p1 = (dc.c1.type == Connector::FIX) ? mss.fixes()[dc.c1.nr].pos : X.row(dc.c1.nr);
        Vec<D> p2 = (dc.c2.type == Connector::FIX) ? mss.fixes()[dc.c2.nr].pos : X.row(dc.c2.nr);
        
        Vec<D> d = p1 - p2;
        double L = norm(d);
        if (L < 1e-12) continue;
        
        Vec<D> n = d / L; 

        // Geometric stiffness due to constraint tension (lambda)
        double lambda_over_L = lambda / L;

        for (size_t i = 0; i < D; i++)
        for (size_t j = 0; j < D; j++)
        {
            // Hessian of the constraint: (lambda/L) * (I - n*n^T)
            // It contributes to dF/dx (Linearization of the constraint force)
            double Hij = lambda_over_L * ((i==j?1.0:0.0) - n(i)*n(j));

            if (dc.c1.type == Connector::MASS) df(dc.c1.nr*D + i, dc.c1.nr*D + j) -= Hij;
            if (dc.c2.type == Connector::MASS) df(dc.c2.nr*D + i, dc.c2.nr*D + j) -= Hij;
            if (dc.c1.type == Connector::MASS && dc.c2.type == Connector::MASS) {
                df(dc.c1.nr*D + i, dc.c2.nr*D + j) += Hij;
                df(dc.c2.nr*D + i, dc.c1.nr*D + j) += Hij;
            }
        }

        // Cross-blocks: dF/dlambda (Gradient^T) and dConstraint/dx (Gradient)
        for (size_t i = 0; i < D; i++) {
            if (dc.c1.type == Connector::MASS) {
                df(dc.c1.nr*D + i, idx_lambda) -= n(i); // dF1 / dLambda
                df(idx_lambda, dc.c1.nr*D + i) += n(i); // dConstraint / dp1
            }
            if (dc.c2.type == Connector::MASS) {
                df(dc.c2.nr*D + i, idx_lambda) += n(i); // dF2 / dLambda
                df(idx_lambda, dc.c2.nr*D + i) -= n(i); // dConstraint / dp2
            }
        }
    }
  }
};

#endif