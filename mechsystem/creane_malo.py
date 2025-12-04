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
    stiffness = 50000.0 # We simulate a strong string (beam)
    L = 1.0 #Easier to visualize
    floors = 10 # How tall we want the crane to be


    # 1. Create the "tower" of the crane
    nodes = []

    for i in range(floors + 1):
        z = i * L
        floor_nodes = [] # we save the corners of each floor
        
        # We put one corner in the origin and define the others in the x>0, y>0 plane
        c0 = mss.addMass([0, 0, z], mass)
        c1 = mss.addMass([L, 0, z], mass)
        c2 = mss.addMass([L, L, z], mass)
        c3 = mss.addMass([0, L, z], mass) #idx

        floor_nodes = [c0, c1, c2, c3]
        nodes.append(floor_nodes)

    for idx in nodes[0]:
        mss.fixMass(idx) #the floor does not move

    # We know connect every floor we have created
    for i in range(floors):
        current_floor = nodes[i]
        next_floor = nodes[i+1]

        for k in range(4):
            idx_curr = current_floor[k]
            idx_next = next_floor[k]
            
            idx_curr_neighbor = current_floor[(k + 1) % 4] 
            idx_next_neighbor = next_floor[(k + 1) % 4] #index for the next corner in the current flor
            # we have add %4 because the loop sums but for example the alst one is 3-0 not 3-4, so we want
            # the values in module 4

            mss.addSpring(idx_curr, idx_next, stiffness, L) # vertical union (like 0 with 4)

            mss.addSpring(idx_next, idx_next_neighbor, stiffness, L) #horizontal union (like 0-1)

            dist_diag = math.sqrt(2) * L # we creatt the X in each face
            mss.addSpring(idx_curr, idx_next_neighbor, stiffness, dist_diag)
            mss.addSpring(idx_curr_neighbor, idx_next, stiffness, dist_diag)

    return mss

if __name__ == "__main__":
    my_crane = build_crane()
    print("Sistema creado. Iniciando simulación...")
    
    # IMPORTANTE: Aquí llamamos a la función que abre la ventana.
    # Si test_mass_spring.py usaba otra cosa, avísame.
    mass_spring.simulate(my_crane)



    # --- 4. AÑADIR VIBRACIÓN (CARGA) ---
    # El ejercicio pide "simulate vibration".
    # Colgamos un peso extra en la punta del brazo para que tire hacia abajo.
    
    print("Añadiendo carga en la punta...")
    # Posición: Al final del brazo y un poco más abajo
    load_pos = (current_x, L/2, z_bot - 1.0) 
    
    # Masa de 10.0 (10 veces más pesada que el resto)
    tip_mass = mss.add(Mass(10.0, load_pos)) 
    
    # Lo enganchamos a la punta del brazo con 2 cables
    cable_len = math.sqrt( (L/2)**2 + 1.0**2 ) # Pitágoras simple
    mss.add(Spring(cable_len, stiffness, (prev_face[0], tip_mass)))
    mss.add(Spring(cable_len, stiffness, (prev_face[3], tip_mass)))
    
    # Damos un empujón inicial a la carga para que oscile
    # Nota: Como no tenemos setVelocity fácil, la gravedad hará el trabajo al soltarla.