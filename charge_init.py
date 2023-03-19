# -*- coding: utf-8 -*-
"""
    Initiated the locations of a charge distribution
    
    """
    
import numpy as np
import matplotlib.pyplot as plt
from charge_class import point_charge
import random

        
    # Makes a 3D domain for charges to be constrained to
def define_grid(xrange = [-10, 10], yrange = [-10, 10], zrange = [-10, 10], res = (1, 1, 1)):
    dx = res[0]
    dy = res[1]
    dz = res[2]
                
    xn = int(np.floor((xrange[1] - xrange[0]) / dx))
    yn = int(np.floor((yrange[1] - yrange[0]) / dy))
    zn = int(np.floor((zrange[1] - zrange[0]) / dz))
    # print(xn, yn, zn)
    coord_grid = np.zeros((xn, yn, zn))
    
   
    # for num in x:
    #     print(num)
    # xarray = []
    # yarray = []
    # zarray = []
    

    for x in range(0, xn):
        coord_grid[x, :, :] = x * dx + xrange[0]
        # xarray.append(x)
        for y in range(0, yn):
            coord_grid[:, y, :] = y * dy + yrange[0]
            # yarray.append(y)
            for z in range(0, zn):
                coord_grid[z, :, :] = z * dz + zrange[0]
                # zarray.append(z)

    x = range(xrange[0], xrange[1] + dx, dx)
    y = range(yrange[0], yrange[1] + dy, dy)
    z = range(zrange[0], zrange[1] + dz, dz)
    
    return coord_grid, x, y, z

    # Initialize a point charge at every location of an inputted list of tuples
def initialize_charges(distrib, ics = None):
    # Initial axis limits set to 1st particle's position
    
    minx = distrib[0][0]
    maxx = distrib[0][0]
    miny = distrib[0][1]
    maxy = distrib[0][1]
    minz = distrib[0][2]
    maxz = distrib[0][2]
    
    
    me = 9.109 * 10 ** -31 # kg  electron mass
    qe = -1.602 * 10 ** -19 # C  electron charge
    
    mass = me
    q = qe
    
    q_id = 1
    q_dict = {}
    for coord in distrib:
        charge_inst = point_charge(q_id, mass, q,  coord, init_vel = [0,0,0])
        q_dict.update({q_id: charge_inst})
        q_id += 1
        
        # Update domain limits
        if coord[0] < minx:
            minx = coord[0]
        elif coord[0] > maxx:
            maxx = coord[0]
            
        if coord[1]< miny:
            miny = coord[1]
        elif coord[1] > maxy:
            maxy = coord[1]
            
        if coord[2]< minz:
            minz = coord[2]
        elif coord[2] > maxz:
            maxz = coord[2]
        
    return q_dict, [[minx, maxx], [miny, maxy], [minz, maxz]]

#########################################################
# Generate specific shapes interpolated to the grid

# Interpolates a solid sphere into coordinates within gridded domain
def sphere(radius, center, x, y, z, res, debug = 0, num = 0):
    interp_xcent = np.round(center[0] / res[0]) * res[0]
    interp_ycent = np.round(center[1] / res[1]) * res[1]
    interp_zcent = np.round(center[2] / res[2]) * res[2]
    
    ncenter = (interp_xcent, interp_ycent, interp_zcent)
    
    avgres = np.sqrt(res[0] ** 2 + res[1] ** 2 + res[2] ** 2)
    
    sphere_coords = []
    
    for xval in x:
        for yval in y:
            for zval in z:
                r = np.sqrt((xval - interp_xcent) ** 2 +
                            (yval - interp_ycent) ** 2 +
                            (zval - interp_zcent) ** 2 )
                # if cell is within 1 average resolution of radius
                if abs(r) < radius + avgres:
                    sphere_coords.append([xval, yval, zval])
                    
    # Remove all but a an inputted number of random particles
    if num != 0:
        randomkeys = random.sample(range(1, len(sphere_coords) + 1), num)
        temp_coords = [sphere_coords[ind] for ind in randomkeys]
        sphere_coords = temp_coords
        
        
    if debug == 1:
        print(f"sphere has points = {len(sphere_coords)}")
        plt.close("all")
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_zlabel("z")
        ax.set_xlim([min(x), max(x)])
        ax.set_ylim([min(y), max(y)]) 
        ax.set_zlim([min(z), max(z)]) 
        for tup in sphere_coords:
            if tup[0] == 3 and tup[2] == 0:
                ax.scatter(tup[0], tup[1], tup[2], )
            
    return sphere_coords, ncenter

# Interpolates a hollow sphere into coordinates within gridded domain
def sphere_shell(radius, center, x, y, z, res, debug = 0, num = 0):
    interp_xcent = np.round(center[0] / res[0]) * res[0]
    interp_ycent = np.round(center[1] / res[1]) * res[1]
    interp_zcent = np.round(center[2] / res[2]) * res[2]
    
    ncenter = (interp_xcent, interp_ycent, interp_zcent)
    
    avgres = np.sqrt(res[0] ** 2 + res[1] ** 2 + res[2] ** 2)
    
    sphere_coords = []
    
    for xval in x:
        for yval in y:
            for zval in z:
                r = np.sqrt((xval - interp_xcent) ** 2 +
                            (yval - interp_ycent) ** 2 +
                            (zval - interp_zcent) ** 2 )
                # if cell is within 1 average resolution of radius
                # if abs(radius - r) < 2*avgres:
                if radius - avgres / 2 < r < radius + avgres / 2:
                    sphere_coords.append([xval, yval, zval])
                    
    # Remove all but a an inputted number of random particles
    if num != 0:
        randomkeys = random.sample(range(1, len(sphere_coords) + 1), num)
        temp_coords = [sphere_coords[ind] for ind in randomkeys]
        sphere_coords = temp_coords
        
    if debug == 1:
        print(f"sphere has points = {len(sphere_coords)}")
        plt.close("all")
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_zlabel("z")
        ax.set_xlim([min(x), max(x)])
        ax.set_ylim([min(y), max(y)]) 
        ax.set_zlim([min(z), max(z)]) 
        for tup in sphere_coords:
            if tup[2] == 0:
            # if tup[0] == 3 and tup[2] == 0:
                ax.scatter(tup[0], tup[1], tup[2], )
        
    return sphere_coords, ncenter

# Initial distribution of 2 charges
def two_charges(res, p1 = [-10, 0, 0], p2 = [10, 0, 0]):
    coords = []
    
    interp_x = np.round(p1[0] / res[0]) * res[0]
    interp_y = np.round(p1[1] / res[1]) * res[1]
    interp_z = np.round(p1[2] / res[2]) * res[2]
    
    coords.append([interp_x, interp_y, interp_z])
    
    interp_x = np.round(p2[0] / res[0]) * res[0]
    interp_y = np.round(p2[1] / res[1]) * res[1]
    interp_z = np.round(p2[2] / res[2]) * res[2]
    
    coords.append([interp_x, interp_y, interp_z])

    return coords

def randomize_properties(q_dict, version = 1, optarg = 0):
    num = len(q_dict)
    # print(num)
    # Will flip charge of half particles at random
    if version == 1:
        randomkeys = random.sample(range(1, num), int(np.ceil(num / 2)))
        for key in randomkeys:
            q_dict[key].charge *= -1
            if q_dict[key].charge < 0:
                q_dict[key].color = "blue"
            elif q_dict[key].charge > 0:
                q_dict[key].color = "red"
            
    # Randomizes mass
    if version == 2:
        for key in q_dict:
            q_dict[key].mass = 10 ** random.randint(-5, 15)
    
    # Remove charge from random sample of half particles
    if version == 3:
        randomkeys = random.sample(range(1, num), int(np.ceil(num / 2)))
        for key in randomkeys:
            q_dict[key].charge = 0
            q_dict[key].color = "grey"
            
    # Remove charge from all particles
    if version == 4:
        for key in q_dict:
            q_dict[key].charge = 0
            q_dict[key].color = "grey"
            
    # Make all masses equal to inputted value or 1000 kg by default
    if version == 5:
        for key in q_dict:
            if optarg == 0:
                q_dict[key].mass = 1000
            else:
                q_dict[key].mass = optarg
                
    # Make all masses zero
    if version == 6:
        for key in q_dict:
            if optarg == 0:
                q_dict[key].mass = 0
                q_dict[key].color = "yellow"
                
    # Set random speeds to half of the particles
    if version == 7:
        randomkeys = random.sample(range(1, num), int(np.ceil(num / 2)))
        for key in randomkeys:
            q_dict[key].vel[0] = random.randint(-10,10)
            q_dict[key].vel[1] = random.randint(-10,10)
            q_dict[key].vel[2] = random.randint(-10,10)


    return
