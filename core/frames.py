import numpy as np
# import random
import math
import tqdm
from ase.data import atomic_numbers, atomic_names, atomic_masses, covalent_radii
from ase.build import molecule
# from xtb.ase.calculator import XTB
from ase import Atoms
# from ase.optimize import BFGS,MDMin,FIRE,GPMin

import core.tools as tools


def check_constraint(d_0, cluster, pos, n_atoms, gamma):
    min = 1.0 - gamma
    max = 1.0 + gamma
    idx=-1
    min_dist = 9999999.99999
    for j in range(len(cluster)):
        dij = np.linalg.norm(cluster[j]-pos)
        if dij < min_dist:
            min_dist = dij
            idx = j

    if idx < 0 or idx >= n_atoms:
        exit("INVALID INDEX!!!")
        
    if min_dist >= 2.0 * min * d_0 and min_dist <= 2.0 * max * d_0:
        return True
    return False

def genSamples(n_structures, shape, atom, n_atoms: int, 
               factor: float,
               gamma: float, 
               folder):
    print(f"\t\t\tBOX: {shape}")
    clusters_list=[]
    d_0=1.0*covalent_radii[atomic_numbers[atom]]

    # print("Raio: ", d_0)
    # exit()
    # This radius comes from http://dx.doi.org/10.1016/S0166-1280(01)00730-8
    sphere_radius =  2.0* d_0 * (factor + math.pow((3.0*n_atoms)/(4.0*np.pi*math.sqrt(2.0)), (1.0/3.0)))

    for n in tqdm.tqdm(range(n_structures)):
        if shape=="CUBE":
            tc = "RC"
            # box = 2.0 * d_0 * math.pow(n_atoms, (1.0 / 3.0))
            box = sphere_radius*math.pow(4.0*np.pi/3.0, (1.0/3.0))
            # print("CUBE: ", box)
            cluster=[]
            for i in range(n_atoms):
                ok = False
                while not ok:
                    pos=box*(np.random.rand(3)-np.ones(3)*0.5)
                    if i==0:
                        ok = True
                    else:
                        ok = check_constraint(d_0, cluster, pos, n_atoms, gamma)
                cluster.append(pos)

        if shape=="SPHERE":
            tc = "RS"
            # This radius comes from http://dx.doi.org/10.1016/S0166-1280(01)00730-8
            box = sphere_radius
            # print("SPHERE: ", box)
            cluster=[]
            for i in range(n_atoms):
                ok = False
                while not ok:
                    pos=box*(np.random.rand(3)-np.ones(3)*0.5)
                    while np.linalg.norm(np.zeros(3)-pos)>box:
                        pos=box*(np.random.rand(3)-np.ones(3)*0.5)
                    if i==0:
                        ok = True
                    else:
                        ok = check_constraint(d_0, cluster, pos, n_atoms, gamma)
                cluster.append(pos)

        # structure = Atoms(atom+str(n_atoms),
        #         positions=cluster,
        #         calculator=XTB(method="GFN1-xTB",max_iterations=1000,accuracy=1000))
        structure = Atoms(atom+str(n_atoms), positions=cluster)

        # /initial

        fname = f"{tc}_{atom}{n_atoms}_F{n}"
        list_atoms = structure.get_chemical_symbols()
        list_coords = structure.get_positions()
        # energyI = structure.get_potential_energy()
        tools.generateXYZ(list_atoms, list_coords, 0.0, fname, folder)

        # dyn = FIRE(structure, logfile="./log.txt")
        # # dyn = GPMin(structure, logfile="./log.txt")
        # # dyn = MDMin(structure, logfile="./log.txt")
        # # dyn = BFGS(structure, logfile="./log.txt")
        # dyn.run(fmax=0.1, steps=100)

        # fname = f"{tc}_{atom}{n_atoms}_{n}"
        # list_atoms = structure.get_chemical_symbols()
        # list_coords = structure.get_positions()
        # energyF = structure.get_potential_energy()
        # tools.generateXYZ(list_atoms, list_coords, energyF, fname, folder)
  


# import random
# import math

# def generate_atoms(n, d, v):
#     if n <= 0 or d <= 0 or v <= 0:
#         raise ValueError("Number of atoms, distance, and volume must be positive")

#     # Calculate the side length of the cube from the volume
#     side_length = v ** (1/3)

#     # Initialize the list of atoms; the first atom can be at any random position within the cube
#     atoms = [(random.uniform(0, side_length), random.uniform(0, side_length), random.uniform(0, side_length))]

#     while len(atoms) < n:
#         # Propose a new atom position within the cube
#         new_atom = (random.uniform(0, side_length), random.uniform(0, side_length), random.uniform(0, side_length))
        
#         # Check distance from all existing atoms
#         if all(distance(new_atom, existing_atom) >= d for existing_atom in atoms):
#             atoms.append(new_atom)

#     return atoms

# def distance(atom1, atom2):
#     return math.sqrt(sum((a - b) ** 2 for a, b in zip(atom1, atom2)))

# # Example usage
# n = 10  # Number of atoms
# d = 1.0 # Minimum distance between atoms
# v = 1000.0 # Volume of the cube

# atoms = generate_atoms(n, d, v)
# print(atoms)
