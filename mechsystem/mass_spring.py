# Define the Mass class (simplified for this example)
class Mass:
    def __init__(self, mass, pos):
        self.mass = mass  # Mass value
        self.pos = list(pos)  # Position as a list (mutable, so you can modify it)
        self.vel = [0.0, 0.0, 0.0]  # Velocity (initialized to zero)
        self.acc = [0.0, 0.0, 0.0]  # Acceleration (initialized to zero)

    def __repr__(self):
        return f"Mass(mass={self.mass}, pos={self.pos}, vel={self.vel})"

# Define the Fix class (for fixed points)
class Fix:
    def __init__(self, pos):
        self.pos = pos  # Fixed position

    def __repr__(self):
        return f"Fix(pos={self.pos})"

# Define the Spring class (simplified for this example)
class Spring:
    def __init__(self, length, stiffness, connectors):
        self.length = length  # Length of the spring
        self.stiffness = stiffness  # Stiffness of the spring
        self.connectors = connectors  # A list of two connectors (Mass or Fix)

    def __repr__(self):
        return f"Spring(length={self.length}, stiffness={self.stiffness}, connectors={self.connectors})"

# Define the DistanceConstraint class
class DistanceConstraint:
    def __init__(self, c1, c2, rest_length):
        self.c1 = c1  # Connector 1 (Mass or Fix)
        self.c2 = c2  # Connector 2 (Mass or Fix)
        self.rest_length = rest_length  # Rest length of the constraint

    def __repr__(self):
        return f"DistanceConstraint(c1={self.c1}, c2={self.c2}, rest_length={self.rest_length})"

# Define the MassSpringSystem3d class
class MassSpringSystem3d:
    def __init__(self):
        self.masses = []  # List of masses
        self.springs = []  # List of springs
        self.fixes = []  # List of fixed points
        self.constraints = []  # List of distance constraints
        self.gravity = (0, 0, -9.81)  # Default gravity in the z direction

    # Add a mass, fix, or spring to the system and return the object
    def add(self, obj):
        print(f"Adding object of type: {type(obj)}")  # Debugging line
        
        if isinstance(obj, Mass):
            self.masses.append(obj)
            print(f"Added Mass: {obj}")  # Debugging line
        elif isinstance(obj, Fix):
            self.fixes.append(obj)
            print(f"Added Fix: {obj}")  # Debugging line
        elif isinstance(obj, Spring):
            self.springs.append(obj)
            print(f"Added Spring: {obj}")  # Debugging line
        else:
            print(f"Invalid object type: {obj}")  # Debugging line
            return None
        
        return obj  # Return the object after adding it

    # Add a distance constraint
    def addDistanceConstraint(self, constraint):
        self.constraints.append(constraint)
        print(f"Added DistanceConstraint: {constraint}")  # Debugging line

    # Get the state of the system (for example, positions of the masses)
    def getState(self):
        return [(m.pos) for m in self.masses]

    # Simulate the system (This is a placeholder for simulation logic)
    def simulate(self, time_step, num_steps):
        # This function should implement the system's physics (e.g., using ODE solvers).
        # For now, we are just printing a message as a placeholder.
        print(f"Simulating for {num_steps} steps with time step {time_step}...")
        
        # Add simulation logic to update positions, velocities, etc.
        for m in self.masses:
            m.pos[0] += m.vel[0] * time_step
            m.pos[1] += m.vel[1] * time_step
            m.pos[2] += m.vel[2] * time_step


