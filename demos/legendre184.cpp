#include <iostream>
#include <vector>
#include <iomanip>
#include "autodiff.hpp"

using namespace ASC_ode;

template <typename T>
void LegendrePolynomials(int n, T x, std::vector<T>& P) {
    // make sur we are introducing a valid number to don't get an error
    if (n < 0) {
        P.clear();
        return;
    }
    P.resize(n + 1); // because we want up to P5 and need to savve all from 0: so 6=5+1

    P[0] = T(1); //known
    if (n == 0) return;

    P[1] = x; // known

    // now we use the given Legendre formula
    for (int k = 2; k <= n; ++k) {
        P[k] = ((T(2 * k - 1) * x * P[k - 1]) - T(k - 1) * P[k - 2]) / T(k);
    }
    // we use T to make sure we don't loose decimals
}

int main() {
    int N_order = 5; // we want P5
    std::vector< AutoDiff<1> > P; // to save the polynomials

    std::cout << std::setw(10) << "x" 
              << std::setw(15) << "P5(x)" 
              << std::setw(15) << "P5'(x)" << std::endl;

    // we use 1.05 to make sure the last value (ain approx 1) is considered, just in case
    // we have something like .000000001 becuase of the decimals and just putig 1.00 the code
    // will not consider the value
    for (double val = -1.0; val <= 1.05; val += 0.2) {
        AutoDiff<1> x = Variable<0>(val);
        
        LegendrePolynomials(N_order, x, P); // since x is Autodiff the function can compute
        // the value and the derivative at the same time

        std::cout << std::fixed << std::setprecision(4) 
                  << std::setw(10) << val 
                  << std::setw(15) << P[5].value()
                  << std::setw(15) << P[5].deriv()[0]
                  << std::endl;
    }

    return 0;
}