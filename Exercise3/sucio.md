# Exercise 3 (Due 09.12.2025)

*All documents we mentioned can be found in branch "Exercise3".

## Exercise 20.4

### Q1: Test the examples

1. Update the upsteam and get folder `mechsystem`
2. Modify Cmakelist:

 ① First, add cmake version, project name and standard

```
cmake_minimum_required(VERSION 3.20)
project(mechsystem)
set (CMAKE_CXX_STANDARD 20)
```
 ② Then, include the head file in the root 'src'

```
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_executable (test_mass_spring mass_spring.cpp ${CMAKE_SOURCE_DIR}/src)
```
③ In the end, use pybind to create python module
```
pybind11_add_module(mass_spring
    bind_mass_spring.cpp
)

target_include_directories(mass_spring PRIVATE 
    ${CMAKE_CURRENT_SOURCE_DIR}
)
```

3. Build and compile in `.venv` enviroment. (Since we install python in .venv enviroment in case of extra conflicts)

```
# activate the enviroment
source .venv/bin/activate

# build file
mkdir build
cd build

# use venv Python to open CMake
cmake -DPython_EXECUTABLE=$(which python) ..

# Compile 
cmake ..
make -j
```

4. We can test the exmaple files:

```
./test_mass_spring > test.txt
```
[Test1: The result link](data_ex3/test.txt)

```
cd ~/ASC-ODE-team09/mechsystem

PYTHONPATH=../build/mechsystem python3 test_mass_spring.py >test2.txt
```

[Test2 : The result link](data_ex3/test2.txt)


5. We can also test the 3D version according to `mass_spring.ipynb`

① Install pythreejs and activate

```
pip install pythreejs

jupyter nbextension enable --py widgetsnbextension
jupyter nbextension enable --py pythreejs
```

② Here we run directly in vscode. Pay attention: change the local repo
```
sys.path.append('../build/mechsystem')
```

<video width="600" controls>
  <source src="_static/ms_nodis_largek.mp4" type="video/mp4">
</video>


In the example code, the spring constant is too large such that we can't see an obviously 'spring' action.

If we change the spring constant with '1000' and '500', we can see:

<video width="600" controls>
  <source src="_static/ms_nodis_smallk.mp4" type="video/mp4">
</video>


### Q2: Add Lagragian constraints

1. In order to add distance constraints to the MassSpring system, which means:

$$
| x_{i} - x_{j} | = L_{0}
$$

A Lagrange multiplier 𝜆 enforces this constraint, leading to the Lagrange function：

$$
L(x,𝜆) = -U(x) + 𝜆g(x)
$$


The system changes from ODE to DAE:

$$
\begin{pmatrix}
M & G^T \\
G & 0
\end{pmatrix}
\begin{pmatrix}
\ddot{x} \\
\lambda
\end{pmatrix}
=
\begin{pmatrix}
F_\text{ext} \\
0
\end{pmatrix},
$$


In the implementation, each constraint contributes:
  - a constraint equation $g(x) = L-L_{0} = 0$
  - a constraint force $F_{i} ​= −λn, F_{j}​ = +λn$, n is unit direction vector



2. Modify **mass_spring.hpp** and classify new function **DistanceConstraint**

```
class DistanceConstraint
{
public:
    Connector c1;
    Connector c2;
    double rest_length;
};

# Add new container
std::vector<DistanceConstraint> m_constraints;

# Add function
void addDistanceConstraint(const DistanceConstraint &dc)
{
    m_constraints.push_back(dc);
}
```
Mianly focus on `MSS_Function`:

```
for (size_t i = 0; i < n_constraints; i++)
    {
        auto& dc = mss.constraints()[i];
        double lambda = x(D * n_masses + i); 

        Vec<D> p1 = (dc.c1.type == Connector::FIX) ? mss.fixes()[dc.c1.nr].pos : xmat.row(dc.c1.nr);
        Vec<D> p2 = (dc.c2.type == Connector::FIX) ? mss.fixes()[dc.c2.nr].pos : xmat.row(dc.c2.nr);
        
        # Vector p2 -> p1
        Vec<D> d = p1 - p2; 
        double L = norm(d);
        if (L < 1e-12) continue;
        
        # Gradient of C with respect to p1
        Vec<D> dir = d / L; 

        # Constraint force
        if (dc.c1.type == Connector::MASS) fmat.row(dc.c1.nr) -= lambda * dir;
        if (dc.c2.type == Connector::MASS) fmat.row(dc.c2.nr) += lambda * dir;

        # Constraint Equation: g(x) = L - L0 = 0
        f(D * n_masses + i) = (L - dc.rest_length);
    }
```

### Q3: Implement the exact derivative
Since we are solving DAE, F(x) is non-linear. We need Jacobian derivative:

$$
J = \frac{∂F}{∂x} 	
$$

The Jacobian has the block form

$$
J =
\begin{pmatrix}
\frac{\partial F}{\partial x} & \frac{\partial F}{\partial \lambda} \\
\frac{\partial g}{\partial x} & 0
\end{pmatrix},
$$

where

-  F(x) : physical forces (gravity + springs + constraint forces)  
-  g(x)=0 : constraint equations  
-  G = \nabla g(x) : gradient of the constraints  
- $\lambda$: Lagrange multipliers (constraint forces)

1. First, do the derivative of the **spring forces**. Taking the derivative yields the geometric stiffness matrix:

$$
K = k\, nn^T + \frac{k(L-L_0)}{L}(I - nn^T).
$$

```
double ninj = n(i)*n(j);
double Kij = k * ninj + force_over_L * ((i==j?1.0:0.0) - ninj);

# Assemble into Jacobian df
df(c1.nr*D + i, c1.nr*D + j) -= Kij;
df(c2.nr*D + i, c2.nr*D + j) -= Kij;
df(c1.nr*D + i, c2.nr*D + j) += Kij;
df(c2.nr*D + i, c1.nr*D + j) += Kij;
```
2. Then, do the derivative of the constraint forces:

```
# Geometric stiffness due to constraint tension (lambda)

double lambda_over_L = lambda / L;


# Hessian of the constraint: (lambda/L) * (I - n*n^T)

if (dc.c1.type == Connector::MASS) df(dc.c1.nr*D + i, dc.c1.nr*D + j) -= Hij;
if (dc.c2.type == Connector::MASS) df(dc.c2.nr*D + i, dc.c2.nr*D + j) -= Hij;
if (dc.c1.type == Connector::MASS && dc.c2.type == Connector::MASS) 

df(dc.c1.nr*D + i, dc.c2.nr*D + j) += Hij;
df(dc.c2.nr*D + i, dc.c1.nr*D + j) += Hij;
```

3. Finally, coupling between x and $\lambda$

```
df(dc.c1.nr*D + i, idx_lambda) -= n(i);  // dF/dlambda
df(idx_lambda, dc.c1.nr*D + i) += n(i);  // dC/dx
```


### Result: 3D Version_Mass-Spring system with Lagrangian Constraint

Add constrains directly in `mass_spring.ipnb` to see the 3D version

1. Delete the spring and change to distance constraint:

```
mss.add (Spring(1, 1000, (f1, mA)))
mss.add (Spring(1, 500, (mA, mB)))
```

2. Using `DistanceConstraint` and assume the fixed distance L = 1:

```
dc1 = DistanceConstraint(mA, mB, 1.0)
mss.addDistanceConstraint(dc1)

dc2 = DistanceConstraint(f1, mA, 1.0)
mss.addDistanceConstraint(dc2)
```
Search for every constraints and find the position:

```
for c in mss.constraints:
        if c.c1.type == 1:                      # FIX
            pA = mss.fixes[c.c1.nr].pos
        else:                                   # MASS
            pA = mss.masses[c.c1.nr].pos

        # get position of second constraint point
        if c.c2.type == 1:
            pB = mss.fixes[c.c2.nr].pos
        else:
            pB = mss.masses[c.c2.nr].pos
        springpos.append ([ pA, pB ])
```

3. When the Lagrange multiplier becomes very large, the Jacobian matrix used in the Newton solver becomes numerically ill-conditioned.

Error **Newton did not converge** exist.

We increase the tolerance and max. iteration to *1e-9* and *20*. And slow down the video *sleep(0.3)*

Finally we get:

<video width="600" controls>
  <source src="_static/ms_la_con.mp4" type="video/mp4">
</video>


Obviously, the constraint connection lines maintain a fixed length, regardless of how the mass blocks move. The segments do not stretch or contract noticeably like a soft spring, making the connections appear as metal rods.




