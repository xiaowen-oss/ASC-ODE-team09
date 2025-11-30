#include <iostream>
#include <iomanip>
#include "autodiff.hpp"

using namespace ASC_ode;

int main()
{
    std::cout << std::fixed << std::setprecision(4);

    double val_x = 2.0;
    double val_y = 4.0;

    AutoDiff<2> x = Variable<0>(val_x);
    AutoDiff<2> y = Variable<1>(val_y);

    std::cout << "--- INITIAL VARIABLES ---" << std::endl;
    std::cout << "x = " << x << std::endl;
    std::cout << "y = " << y << std::endl << std::endl;

    std::cout << "--- TEST 1: x - y ---" << std::endl;
    AutoDiff<2> sub = x - y;
    std::cout << "Calculated: " << sub << std::endl;
    std::cout << "Expected:   Value: -2.0000, Deriv: [1.0000, -1.0000]" << std::endl;
    
    std::cout << "--- TEST 2: -x ---" << std::endl;
    AutoDiff<2> neg = -x;
    std::cout << "Calculated: " << neg << std::endl;
    std::cout << "Expected:   Value: -2.0000, Deriv: [-1.0000, 0.0000]" << std::endl << std::endl;

    std::cout << "--- TEST 3: DIVISION x / y ---" << std::endl;
    AutoDiff<2> div = x / y;
    std::cout << "Calculated: " << div << std::endl;
    std::cout << "Expected:   Value: 0.5000, Deriv: [0.2500, -0.1250]" << std::endl << std::endl;

    std::cout << "--- TEST 4 COS(x) ---" << std::endl;
    AutoDiff<2> c = cos(x); 
    std::cout << "Calculated: " << c << std::endl;
    std::cout << "Expected:   Value: -0.4161, Deriv: [-0.9093, 0.0000]" << std::endl;

    std::cout << "--- TEST 5: EXP(x) ---" << std::endl;
    AutoDiff<2> e = exp(x);
    std::cout << "Calculated: " << e << std::endl;
    std::cout << "Expected:   Value: 7.3891, Deriv: [7.3891, 0.0000]" << std::endl;

    std::cout << "--- TEST 6: LOG(x) ---" << std::endl;
    AutoDiff<2> l = log(x);
    std::cout << "Calculated: " << l << std::endl;
    std::cout << "Expected:   Value: 0.6931, Deriv: [0.5000, 0.0000]" << std::endl;

    return 0;
}