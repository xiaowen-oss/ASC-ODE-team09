#include "mass_spring.hpp"
#include "Newmark.hpp"
#include <iomanip>

int main()
{
    std::cout << "==== Testing Mass-Spring System WITH Distance Constraints ====\n";

    MassSpringSystem<3> mss;
    mss.setGravity({0, 0, -9.81});

    // --- Create masses ---
    auto mA = mss.addMass( Mass<3>{ 1.0, {1.0, 0.0, 0.0} } );
    auto mB = mss.addMass( Mass<3>{ 1.0, {2.2, 0.0, 0.0} } ); // intentionally not at rest length

    // --- Add a spring ---
    mss.addSpring( {1.0, 10000.0, {mA, mB}} );

    // --- Add a distance constraint between A and B ---
    DistanceConstraint dc;
    dc.c1 = mA;
    dc.c2 = mB;
    dc.rest_length = 1.0;   // desired distance
    mss.addDistanceConstraint(dc);

    std::cout << "System initialized.\n";
    std::cout << "Initial positions:\n";
    std::cout << " A = " << mss.masses()[mA.nr].pos << "\n";
    std::cout << " B = " << mss.masses()[mB.nr].pos << "\n\n";

    double tend = 5.0;
    int steps = 200;
    double dt = tend / steps;

    Vector<> x(3 * mss.masses().size());
    Vector<> dx(3 * mss.masses().size());
    Vector<> ddx(3 * mss.masses().size());

    mss.getState(x, dx, ddx);

    auto f = std::make_shared<MSS_Function<3>>(mss);
    auto M = std::make_shared<IdentityFunction>(x.size());

    // --- Time stepping ---
    for (int i = 0; i <= steps; i++)
    {
        double t = i * dt;

        // Print positions
        auto posA = mss.masses()[mA.nr].pos;
        auto posB = mss.masses()[mB.nr].pos;
        double dist = norm(posA - posB);
        double phi = dist - dc.rest_length;

        std::cout << std::fixed << std::setprecision(4);
        std::cout << "t = " << t
                  << " | dist = " << dist
                  << " | phi = " << phi << "\n";

        // simulate
        SolveODE_Newmark(dt, 1, x, dx, f, M);

        mss.setState(x, dx, ddx);
    }

    std::cout << "Simulation finished.\n";
    return 0;
}
