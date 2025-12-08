### Q4: Experiments - Crane Structure and Vibration

In this section, we build a complex mechanical structure (a crane) by using the Mass-Spring System combined with Lagrangian Constraints.

1. **Model Construction Strategy**

   We implemented a `build_crane` function in Python. The structure consists of two main types of connections:
   * Rigid Beams: Modeled using `DistanceConstraint`. We want the frame of the crane to be rigid and not oscillate like jelly.
   * Cables (Elastic): Modeled using `Spring`. We want to observe vibration when lifting a mass.

   We defined a helper function to easily switch between these two types:

```python
def add_beam(mss, node1, node2, length, stiffness, use_constraint=True):
    if use_constraint:
        # Use Lagrangian constraint for rigid parts (Steel)
        dc = DistanceConstraint(node1, node2, length)
        mss.addDistanceConstraint(dc)
    else:
        # Use Spring for elastic parts (Cables/Rubber)
        mss.add(Spring(length, stiffness, (node1, node2)))
```

2. **Building the Geometry**

The crane consists of a vertical tower (10 floors) and a horizontal arm (6 units). To ensure structural stability we added crossed diagonal supports to every face.

① The tower loop:
```python
for i in range(floors):
    # ... (node creation logic) ...
    # connecting floors
    add_beam(mss, c_curr, c_next, L, 0, use_constraint=True) #vertical union
    add_beam(mss, c_next, c_next_neighbor, L, 0, use_constraint=True)  #horizontal union
    
    diag_len = math.sqrt(2) * L 
    add_beam(mss, c_curr, c_next_neighbor, diag_len, stiffness_beam, use_constraint=False)  #diagonal (X shape) for stability
```
The stiffness for the structural parts is 0 because their rigidity will be actually being defined by the Lagrange Constraints, while for the elastic springs we fixed a value `stiffness_beam= 5000`. We coded the arm loop analougsly.

② The Load and Cables:

We attach a heavy mass ($m=5.0$) to the tip of the arm using elastic springs (stiffness_cable = 3000.0). This way we introduce the vibration.

```python
# Cables attached to the tip of the arm
mss.add(Spring(dist_h, stiffness_cable, (tip_node_1, load_mass)))
mss.add(Spring(dist_h, stiffness_cable, (tip_node_2, load_mass)))
```

3. **Simulation and Visualization**
   
To visualize the movement smoothly, we modified the simulation loop in the Jupyter Notebook. Since the code with the Lagrange Constraints is computationally heavy, we decouple the physics steps from the rendering steps using `steps_per_frame` so the computer do not struggle that much to calculate the physics and draw the crane quickly at the same time.

```python
steps_per_frame = 2  #render 1 frame every 2 physics calculations

for i in range(1000):
    mss.simulate(0.01, 2) 

    if i % steps_per_frame == 0:
        # we just actualise the graphic positions when we plot (1 over 2)
        sleep(0.01)
```

4. **Result**

The simulation shows the crane structure remaining rigid due to the DistanceConstraints (orange lines), while the load and arm bounce realistically due to the Springs (cyan lines).

https://github.com/user-attachments/assets/34ebc42f-3653-46d9-b584-9e0fe426b58c
