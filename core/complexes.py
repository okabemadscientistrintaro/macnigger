import ase
import numpy as np
from scipy.spatial import ConvexHull
import random
import math

# import tqdm
# import glob
# import math
# import ast
# import os

def fibonacci_sphere(sites:int):
    points = []
    phi = np.pi * (3. - np.sqrt(5.))  # golden angle in radians

    if sites == 1:
         return np.array([0.0,0.0,1.0])
    
    for i in range(sites):
        y = 1 - (i / float(sites - 1)) * 2  # y goes from 1 to -1
        radius = np.sqrt(1 - y * y)  # radius at y

        theta = phi * i  # golden angle increment

        x = np.cos(theta) * radius
        z = np.sin(theta) * radius

        points.append((x, y, z))

    return np.array(points)


# def move_towards_mass_center(center, core_coords, point, distance_to_surface):
#     dist = np.linalg.norm(core_coords-point, axis=1)
#     pcore = core_coords[np.argmin(dist)] # Core atom with the lowest distance to the siet
#     dp = center - pcore
#     dsite = np.dot(dp,point)    
#     core_dist_center = abs(dsite) / np.linalg.norm(dp) 
#     point_dist_center = np.linalg.norm(point - center)
#     direction_to_mass_center = (point - center) / point_dist_center
#     move_distance =  point_dist_center - core_dist_center - distance_to_surface
#     new_position = point - move_distance * direction_to_mass_center
#     return new_position

def move_towards_mass_center(center, core_coords, point, distance_to_surface):
    dist = np.linalg.norm(core_coords-point, axis=1)
    pcore = core_coords[np.argmin(dist)] # Core atom with the lowest distance to the siet
    project = (np.dot(pcore, point) / np.dot(point, point)) * point # projection of the first atom to the site vector
    norm_project = np.linalg.norm(project)
    dir_project = project / norm_project
    dsite = norm_project + distance_to_surface
    new_position = dir_project*dsite
    return new_position


def adjust_sites(core_coords, sites, ligdist, deformation):
    new_sites = []
    center = np.mean(core_coords, axis=0)
    dists = np.linalg.norm(core_coords - center, axis=1)
    radius = np.max(dists)
    # print("Radius: ", radius)
    nsites = []
    for site in sites:
        nsites.append((radius+ligdist)*site)
    nsites = np.array(nsites)    

    if deformation:
        for site in nsites:
            ns = move_towards_mass_center(center, core_coords, site, ligdist)
            new_sites.append(ns)
        return np.array(new_sites), radius
    else:
        return np.array(nsites), radius

def adjust_sites_oriented(core_coords, sites, ligdist, deformation):
    core_sel = [True]*len(core_coords)
    core_ids = []
    center = np.mean(core_coords, axis=0)
    dists = np.linalg.norm(core_coords - center, axis=1)
    radius = np.max(dists)

    nsites = []
    for i, site in enumerate(sites):
        site = site*radius
        dists = np.linalg.norm(core_coords - site, axis=1)
        sel_atom = np.argmin(dists)
        while core_sel[sel_atom] is False:
            dists[sel_atom] = 9999.9
            sel_atom = np.argmin(dists)
        core_ids.append(sel_atom)
        core_sel[sel_atom] = False
        nsites.append((radius+ligdist)*core_coords[sel_atom]/np.linalg.norm(core_coords[sel_atom]))
        # print(dists)
    nsites = np.array(nsites)    
    # print(nsites)

    if deformation:
        new_sites = []
        for i, site in enumerate(nsites):
            # direction_to_mass_center = site / np.linalg.norm(site)
            norm = np.linalg.norm(core_coords[core_ids[i]])
            new_position = (norm+ligdist)*(core_coords[core_ids[i]]/norm) 
            new_sites.append(new_position)
        return np.array(new_sites), radius

    else:
        return np.array(nsites), radius

def check_distribution(sites):
    center = np.mean(sites, axis=0)
    # print("Center: ", center)
    # print("Max_to_center: ", np.max(np.linalg.norm(sites - center, axis=1)))
    pdist = []
    for site in sites:
        dist = np.linalg.norm(sites - site, axis=1)
        dist[np.argmin(dist)] = np.max(dist)
        pdist.append(np.min(dist))
    print("Mean: ", np.mean(pdist))
    print("Max: ", np.max(pdist))
    print("Min: ", np.min(pdist))
    print("STD: ", np.std(pdist))



def rotate_atoms(atoms):
    theta = np.random.uniform(0, 2*np.pi)  # Random rotation angle
    phi = np.random.uniform(0, 2*np.pi)    # Random angle for axis direction
    z = np.random.uniform(0, 1)

    r = np.sqrt(1-z)
    x, y = r * np.cos(phi), r * np.sin(phi)

    # Rotation axis (unit vector)
    axis = np.array([x, y, z])
    axis = axis / np.linalg.norm(axis)

    # Axis-angle to rotation matrix
    cos_theta = np.cos(theta)
    sin_theta = np.sin(theta)
    one_minus_cos = 1 - cos_theta
    ux, uy, uz = axis

    # Constructing the rotation matrix
    rotation_matrix = np.array([
        [cos_theta + ux**2 * one_minus_cos, ux*uy*one_minus_cos - uz*sin_theta, ux*uz*one_minus_cos + uy*sin_theta],
        [uy*ux*one_minus_cos + uz*sin_theta, cos_theta + uy**2 * one_minus_cos, uy*uz*one_minus_cos - ux*sin_theta],
        [uz*ux*one_minus_cos - uy*sin_theta, uz*uy*one_minus_cos + ux*sin_theta, cos_theta + uz**2 * one_minus_cos]
    ])

    atoms = atoms.dot(rotation_matrix)
    return atoms


# def optimize_sites(sites):
#     num_p = len(sites)
#     for t in range(100):
#         step = np.zeros((num_p,3))
#         dmin = float('inf')
#         # ddd = []
#         for i in range(0,num_p-1):
#             for j in range(i+1,num_p):
#                 diff = sites[i,:] - sites[j,:]
#                 norm = np.linalg.norm(diff)
#                 # ddd.append(norm)
#                 dist = np.exp(-norm*4)
#                 step[i,:] += dist*diff / norm
#                 step[j,:] -= dist*diff / norm
#                 if norm<dmin: dmin=norm

#         sites += step
#         # soma = 0
#         for i in range(num_p):
#             norm = np.linalg.norm(sites[i,:])
#             sites[i,:] /= norm
#         #     soma += norm
#         # print("T: "+ str(t)+" -> "+str(soma/num_p)+" Step: "+str(dmin))
#         # print(f"T: {t}, {np.min(ddd)}, {np.mean(ddd)}, {np.std(ddd)}")

#     # return particles
#     return sites

def optimize_sites(sites):
    num_p = len(sites)
    for t in range(10):
        step = np.zeros((num_p,3))
        dmin = float('inf')

        # pdist = []
        # for idx, site in enumerate(sites):
        #     dist = np.linalg.norm(sites - site, axis=0)
        #     dist = np.delete(dist, idx)
    
        for i in range(0,num_p-1):
            for j in range(i+1,num_p):
                diff = sites[i,:] - sites[j,:]
                norm = np.linalg.norm(diff)
                # ddd.append(norm)
                dist = np.exp(-norm*4)
                step[i,:] += dist*diff / norm
                step[j,:] -= dist*diff / norm
                if norm<dmin: dmin=norm

        sites += step/num_p
        for i in range(num_p):
            norm = np.linalg.norm(sites[i,:])
            sites[i,:] /= norm
        #     soma += norm
        # print("T: "+ str(t)+" -> "+str(soma/num_p)+" Step: "+str(dmin))
        # print(f"T: {t}, {np.min(ddd)}, {np.mean(ddd)}, {np.std(ddd)}")

    # return particles
    return sites



def center_mol(coords):
    centroid = np.mean(coords, axis=0)
    translated_particles = coords - centroid
    return translated_particles




def rotation_matrix_from_vectors(vec1, vec2):
	""" Find the rotation matrix that aligns vec1 to vec2
	:param vec1: A 3d "source" vector
	:param vec2: A 3d "destination" vector
	:return mat: A transform matrix (3x3) which when applied to vec1, aligns it with vec2.
	"""
	a, b = (vec1 / np.linalg.norm(vec1)).reshape(3), (vec2 / np.linalg.norm(vec2)).reshape(3)
	#a, b = vec1, vec2
	v = np.cross(a, b)
	c = np.dot(a, b)
	s = np.linalg.norm(v)
	kmat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
	rotation_matrix = np.eye(3) + kmat + np.dot(kmat,kmat) * ((1 - c) / (s ** 2))
	return rotation_matrix

def randomdir():

    theta = np.random.uniform(0, 2 * np.pi)  # Azimuthal angle
    phi = np.random.uniform(0, np.pi)  # Polar angle
    
    x = np.sin(phi) * np.cos(theta)
    y = np.sin(phi) * np.sin(theta)
    z = np.cos(phi)
    return np.array([x,y,z])

	# theta=math.asin(random.uniform(-1.0, 1.0))
	# phi=random.uniform(0.0,2.0*math.pi)
	# x=math.cos(theta)*math.cos(phi)
	# y=math.sin(theta)
	# z=math.cos(theta)*math.sin(phi)

def positining_ligand(coords, direction, dist, core_coords, Dir):
	#coords - coordenares for the ligant
	#direction - directions defined by the functionalisation site
	#dist - bonding distance
	#RC - cluster radius
	#Dir - element of LIGANDS_ORIENTATION
    avcoords=np.mean(coords, axis=0)
    coords = coords - avcoords

    if len(Dir)==1:
        if Dir[0]==-1:
            orient = randomdir()
        elif Dir[0]==0:
            orient = (coords[1,:] - coords[0,:])
            orient = orient / np.linalg.norm(orient) 
        else:
            orient = (coords[0,:] - coords[1,:])
            orient = orient / np.linalg.norm(orient) 
    else:
        orient = -np.cross(np.array(coords[Dir[1]])-np.array(coords[Dir[0]]),np.array(coords[Dir[2]])-np.array(coords[Dir[0]]))
    
    matrix=rotation_matrix_from_vectors(direction,orient)
    aux=np.dot(coords,matrix)

    center_mol = np.mean(aux, axis=0)
    dist_atoms = np.linalg.norm(aux - center_mol, axis=1)
    radius = np.max(dist_atoms)

    aux= np.array([(p+direction*(dist+radius*2)) for p in aux])

    dist_atoms = np.linalg.norm(aux - direction*dist, axis=1)
    pmin_site = np.argmin(dist_atoms)

    # need to check sobreposition to all atoms in the core

    aux= np.array([(p+(-direction)*(dist_atoms[pmin_site])) for p in aux])
    return aux

# def positining_ligand(coords, direction, dist, RC, Dir):
# 	#coords - coordenares for the ligant
# 	#direction - directions defined by the functionalisation site
# 	#dist - bonding distance
# 	#RC - cluster radius
# 	#Dir - element of LIGANDS_ORIENTATION
#     avcoords=np.sum(coords, axis=0)
#     centered=[]
#     RL=0.000
#     for i in range(len(coords)):
#         v=[coords[i,j]-avcoords[j] for j in range(3)]
#         print(v)
#         centered.append(v)
#         for k in range(len(coords)):
#             v=[coords[i,j]-coords[k,j] for j in range(3)]
#             if RL<np.linalg.norm(v) : RL=np.linalg.norm(v)/2.0
#     if len(Dir)==1:
#         if Dir[0]==-1:
#             orient = randomdir()
#         elif Dir[0]==0:
#             orient = (coords[1,:] - coords[0,:])
#             orient = orient / np.linalg.norm(orient) 
#         else:
#             orient = (coords[0,:] - coords[1,:])
#             orient = orient / np.linalg.norm(orient) 
#     else:
#         orient = -np.cross(np.array(coords[Dir[1]])-np.array(coords[Dir[0]]),np.array(coords[Dir[2]])-np.array(coords[Dir[0]]))
    
#     matrix=rotation_matrix_from_vectors(direction,orient)
#     aux=np.dot(centered,matrix)
#     norm=np.linalg.norm(direction)
#     to_append=[(p+direction*(1.0/norm)*(dist+RL)) for p in aux]
#     return np.array(to_append)
