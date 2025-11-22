#include <iostream>
#include <fstream> 

#include <nonlinfunc.hpp>
#include <timestepper.hpp>

using namespace ASC_ode;


class MassSpring : public NonlinearFunction
{
private:
  double mass;
  double stiffness;

public:
  MassSpring(double m, double k) : mass(m), stiffness(k) {}

  size_t dimX() const override { return 2; }
  size_t dimF() const override { return 2; }
  
  void evaluate (VectorView<double> x, VectorView<double> f) const override
  {
    f(0) = x(1);
    f(1) = -stiffness/mass*x(0);
  }
  
  void evaluateDeriv (VectorView<double> x, MatrixView<double> df) const override
  {
    df = 0.0;
    df(0,1) = 1;
    df(1,0) = -stiffness/mass;
  }
};


int main(int argc, char* argv[])
{
    if (argc != 2) return 1;

    int steps = std::stoi(argv[1]);
    double T = 4.0 * M_PI;
    double tau = T / steps;

    Vector<> y_exp = {1.0, 0.0};
    Vector<> y_imp = {1.0, 0.0};
    Vector<> y_impv = {1.0, 0.0};

    auto rhs = std::make_shared<MassSpring>(1.0, 1.0);

    ExplicitEuler  e_exp(rhs);
    ImplicitEuler  e_imp(rhs);
    ImprovedEuler  e_impv(rhs);

    std::string fe = "compare_results/explicit_" + std::to_string(steps) + ".txt";
    std::string fi = "compare_results/implicit_" + std::to_string(steps) + ".txt";
    std::string fm = "compare_results/improved_" + std::to_string(steps) + ".txt";

    std::ofstream out_e(fe);
    std::ofstream out_i(fi);
    std::ofstream out_m(fm);

    out_e << 0.0 << " " << y_exp(0) << " " << y_exp(1) << " " << steps << "\n";
    out_i << 0.0 << " " << y_imp(0) << " " << y_imp(1) << " " << steps << "\n";
    out_m << 0.0 << " " << y_impv(0) << " " << y_impv(1) << " " << steps << "\n";

    for (int j = 0; j < steps; j++)
    {
        double t = (j + 1) * tau;

        e_exp.DoStep(tau, y_exp);
        e_imp.DoStep(tau, y_imp);
        e_impv.DoStep(tau, y_impv);

        out_e << t << " " << y_exp(0) << " " << y_exp(1) << " " << steps << "\n";
        out_i << t << " " << y_imp(0) << " " << y_imp(1) << " " << steps << "\n";
        out_m << t << " " << y_impv(0) << " " << y_impv(1) << " " << steps << "\n";
    }

    return 0;
}
