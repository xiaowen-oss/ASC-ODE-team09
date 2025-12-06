#include <sstream>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include "mass_spring.hpp"
#include "Newmark.hpp"

namespace py = pybind11;

PYBIND11_MAKE_OPAQUE(std::vector<Mass<3>>);
PYBIND11_MAKE_OPAQUE(std::vector<Fix<3>>);
PYBIND11_MAKE_OPAQUE(std::vector<Spring>);
PYBIND11_MAKE_OPAQUE(std::vector<DistanceConstraint>);

PYBIND11_MODULE(mass_spring, m) {
    m.doc() = "mass-spring-system simulator"; 

    // --- 2D CLASSES ---
    py::class_<Mass<2>> (m, "Mass2d")
          .def_property("mass",
                    [](Mass<2> & m) { return m.mass; },
                    [](Mass<2> & m, double mass) { m.mass = mass; })
      .def_property_readonly("pos",
                             [](Mass<2> & m) { return m.pos.data(); });
    ;
    m.def("Mass", [](double m, std::array<double,2> p) {
      return Mass<2>{m, { p[0], p[1] }};
    });
    
    // --- 3D CLASSES ---
    py::class_<Mass<3>> (m, "Mass3d")
      .def_property("mass",
                    [](Mass<3> & m) { return m.mass; },
                    [](Mass<3> & m, double mass) { m.mass = mass; })
      .def_property_readonly("pos",
                             [](Mass<3> & m) { return m.pos.data(); });
    ;
    m.def("Mass", [](double m, std::array<double,3> p) {
      return Mass<3>{m, { p[0], p[1], p[2] }};
    });

    // --- FIXES ---
    py::class_<Fix<2>> (m, "Fix2d")
      .def_property_readonly("pos", [](Fix<2> & f) { return f.pos.data(); });
    m.def("Fix", [](std::array<double,2> p) { return Fix<2>{ { p[0], p[1] } }; });

    py::class_<Fix<3>> (m, "Fix3d")
      .def_property_readonly("pos", [](Fix<3> & f) { return f.pos.data(); });
    m.def("Fix", [](std::array<double,3> p) { return Fix<3>{ { p[0], p[1], p[2] } }; });

    // --- CONNECTORS (CRITICAL FIX: exposing properties) ---
    py::class_<Connector> (m, "Connector")
      .def_readwrite("nr", &Connector::nr)
      // Expose 'type' as an integer so Python can read/compare it easily
      .def_property("type", 
            [](const Connector &c) { return (int)c.type; },
            [](Connector &c, int t) { c.type = (Connector::CONTYPE)t; }
      );

    // --- SPRING ---
    py::class_<Spring> (m, "Spring")
      .def(py::init<double, double, std::array<Connector,2>>())
      .def_property_readonly("connectors", [](Spring & s) { return s.connectors; });

    // --- CONSTRAINTS ---
    py::class_<DistanceConstraint>(m, "DistanceConstraint")
    .def(py::init<Connector, Connector, double>())
    .def_readwrite("c1", &DistanceConstraint::c1)
    .def_readwrite("c2", &DistanceConstraint::c2)
    .def_readwrite("rest_length", &DistanceConstraint::rest_length);
    
    py::bind_vector<std::vector<DistanceConstraint>>(m, "DistanceConstraints");
    py::bind_vector<std::vector<Mass<3>>>(m, "Masses3d");
    py::bind_vector<std::vector<Fix<3>>>(m, "Fixes3d");
    py::bind_vector<std::vector<Spring>>(m, "Springs");        
    
    // --- 2D SYSTEM ---
    py::class_<MassSpringSystem<2>> (m, "MassSpringSystem2d")
      .def(py::init<>())

      .def("add", [](MassSpringSystem<2> & mss, Mass<2> m) { return mss.addMass(m); });
      
    // --- 3D SYSTEM ---
    py::class_<MassSpringSystem<3>> (m, "MassSpringSystem3d")
      .def(py::init<>())

      .def("__str__", [](MassSpringSystem<3> & mss) {
        std::stringstream sstr;
        sstr << mss;
        return sstr.str();
      })

      .def_property("gravity", [](MassSpringSystem<3> & mss) { return mss.getGravity(); },
                    [](MassSpringSystem<3> & mss, std::array<double,3> g) { mss.setGravity(Vec<3>{g[0],g[1],g[2]}); })

      .def("add", [](MassSpringSystem<3> & mss, Mass<3> m) { return mss.addMass(m); })
      .def("add", [](MassSpringSystem<3> & mss, Fix<3> f) { return mss.addFix(f); })
      .def("add", [](MassSpringSystem<3> & mss, Spring s) { return mss.addSpring(s); })

      .def_property_readonly("masses", [](MassSpringSystem<3> & mss) -> auto& { return mss.masses(); })
      .def_property_readonly("fixes", [](MassSpringSystem<3> & mss) -> auto& { return mss.fixes(); })
      .def_property_readonly("springs", [](MassSpringSystem<3> & mss) -> auto& { return mss.springs(); })

      .def("__getitem__", [](MassSpringSystem<3> mss, Connector & c) {
        if (c.type==Connector::FIX) return py::cast(mss.fixes()[c.nr]);
        else return py::cast(mss.masses()[c.nr]);
      })
      
      .def("getState", [] (MassSpringSystem<3> & mss) {
        Vector<> x(3*mss.masses().size());
        Vector<> dx(3*mss.masses().size());
        Vector<> ddx(3*mss.masses().size());
        mss.getState (x, dx, ddx);
        return std::vector<double>(x);
      })
      
      .def("addDistanceConstraint",
           [](MassSpringSystem<3>& mss, DistanceConstraint dc){
                mss.addDistanceConstraint(dc);
           })
      .def_property_readonly("constraints",
           [](MassSpringSystem<3>& mss) -> auto& { return mss.constraints(); })


      // --- SIMULATION (MODIFIED FOR DAE/LAGRANGE) ---
      .def("simulate", [](MassSpringSystem<3> & mss, double tend, size_t steps) {
        // Augmented dimensions (Mass DOFs + Constraint DOFs)
        size_t n_mass_dofs = 3 * mss.masses().size();
        size_t n_constraints = mss.constraints().size();
        size_t n_total = n_mass_dofs + n_constraints;

        Vector<> x(n_total);
        Vector<> dx(n_total);
        Vector<> ddx(n_total);
        x = 0.0; dx = 0.0; ddx = 0.0;

        // Retrieve current physical state
        Vector<> x_mass(n_mass_dofs);
        Vector<> dx_mass(n_mass_dofs);
        Vector<> ddx_mass(n_mass_dofs);
        mss.getState(x_mass, dx_mass, ddx_mass);

        // Copy state to augmented vectors
        for(size_t i=0; i<n_mass_dofs; i++) {
            x(i) = x_mass(i);
            dx(i) = dx_mass(i);
            ddx(i) = ddx_mass(i);
        }
        // Multipliers (at the end of x) are initialized to 0

        auto mss_func = std::make_shared<MSS_Function<3>> (mss);
        // Use SystemMassFunction (which puts zeros on diagonal for multipliers)
        auto mass = std::make_shared<SystemMassFunction<3>> (mss);

        // Solve using Generalized Alpha (rho_inf=0.8 damps high freq noise)
        SolveODE_Alpha(tend, steps, 0.8, x, dx, ddx, mss_func, mass);

        // Copy physical state back to mss
        for(size_t i=0; i<n_mass_dofs; i++) {
            x_mass(i) = x(i);
            dx_mass(i) = dx(i);
            ddx_mass(i) = ddx(i);
        }
        mss.setState (x_mass, dx_mass, ddx_mass);  
    });
}
