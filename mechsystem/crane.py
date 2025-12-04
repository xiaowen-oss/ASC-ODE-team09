import math
import sys
# Asegúrate de que esta ruta apunta a donde tienes el archivo .so compilado
sys.path.append('../build/mechsystem') 

from mass_spring import *

def build_crane():
    mss = MassSpringSystem3d()
    mss.gravity = (0.0, 0.0, -9.81) # forces, assume only gravity

    # Parameters
    mass = 1.0    
    stiffness = 2000.0 # We simulate a strong string (beam)
    L = 1.0 #Easier to visualize
    floors = 10 # How tall we want the crane to be

    # 1. Create the "tower" of the crane
    nodes = []

    for i in range(floors + 1):
        z = i * L
        floor_nodes = []
        
        coords = [    # We put one corner in the origin and define the others in the x>0, y>0 plane
            (0, 0, z),
            (L, 0, z), 
            (L, L, z), 
            (0, L, z)   
        ]

        for pos in coords:
            if i == 0:
                c = mss.add(Fix(pos)) #the floor 0 does not move
            else:
                c = mss.add(Mass(mass, pos))

            floor_nodes.append(c) # we save the data corners

        nodes.append(floor_nodes)

    # We now connect every floor we have created
    for i in range(floors):
        current_floor = nodes[i]
        next_floor = nodes[i+1]

        for k in range(4):
            c_curr = current_floor[k]
            c_next = next_floor[k]
            
            c_curr_neighbor = current_floor[(k + 1) % 4] 
            c_next_neighbor = next_floor[(k + 1) % 4] #index for the next corner in the current flor
            # we have add %4 because the loop sums but for example the alst one is 3-0 not 3-4, so we want
            # the values in module 4

            mss.add(Spring(L, stiffness, (c_curr, c_next))) # vertical union (like 0 with 4)

            mss.add(Spring(L, stiffness, (c_next, c_next_neighbor)))  #horizontal union (like 0-1)

            diag_len = math.sqrt(2) * L  # we creatt the X in each face
            mss.add(Spring(diag_len, stiffness, (c_curr, c_next_neighbor)))
            mss.add(Spring(diag_len, stiffness, (c_curr_neighbor, c_next)))

    # 2. We constract the "arm" of the crane
    arm_length = 8  
    
    prev_face = [ #we "glue" the arm in the right face of last cube of the tower (floor m and floor m-1)
        nodes[floors-1][1],
        nodes[floors-1][2],
        nodes[floors][2],  
        nodes[floors][1] 
    ]
    
    for j in range(1, arm_length + 1):
        current_x = (L) + (j * L) 
        z_bot = (floors - 1) * L
        z_top = floors * L
        
        arm_coords = [
            (current_x, 0, z_bot),
            (current_x, L, z_bot),
            (current_x, L, z_top),
            (current_x, 0, z_top) 
        ]
        
        current_face = []
        for pos in arm_coords:
            c = mss.add(Mass(mass, pos)) 
            current_face.append(c)
            
        for k in range(4):
            c_old = prev_face[k]
            c_new = current_face[k]
            c_new_neighbor = current_face[(k+1)%4]
            c_old_neighbor = prev_face[(k+1)%4]  
            
            mss.add(Spring(L, stiffness, (c_old, c_new))) # in the x direction 
            
            mss.add(Spring(L, stiffness, (c_new, c_new_neighbor))) # the others (z and y)
            
            diag = math.sqrt(2) * L #fro the X structure
            mss.add(Spring(diag, stiffness, (c_old, c_new_neighbor)))
            mss.add(Spring(diag, stiffness, (c_old_neighbor, c_new)))
            
        prev_face = current_face

    # 3. Vibration
    
    # Lo ponemos en la punta (X=current_x), en el medio (Y=L/2) y colgando un poco (Z-1)
    load_pos = (current_x, L/2, z_bot - 1.0) #point where we are going to put the mass
    
    load_mass = mss.add(Mass(20.0, load_pos)) #we create the mass charge
    
    # La colgamos con cables desde la punta del brazo
    # Conectamos a las esquinas de abajo del final del brazo (prev_face[0] y prev_face[1])
    # Calculamos la distancia (Pitágoras)
    cable_len = math.sqrt( (L/2)**2 + 1.0**2 )
    
    mss.add(Spring(cable_len, stiffness/2, (prev_face[0], load_mass))) #more elastic so they rebote more
    mss.add(Spring(cable_len, stiffness/2, (prev_face[1], load_mass)))

    return mss

if __name__ == "__main__":
    my_crane = build_crane()  
    my_crane.simulate(0.1, 1000)
    
    print("Final position of the higher mass:") #to check that is has moved and the code works
    print(my_crane.masses[-1].pos)