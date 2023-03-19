# -*- coding: utf-8 -*-
"""
 Main function to generate a charge distribution and have it propagate in time
"""

import charge_init as ci
import charge_class as qc
import plot_nbody as pltn
import numpy as np


if __name__ == "__main__":

    #############################################################
    # 2 particles
    ############################################################
    coords = ci.two_charges(res= (1, 1, 1), p1 = [-10, 0, 0], p2 = [10, 0, 0])
    # ci.randomize_properties(q_dict, version = 1) # flip charge of half particles
    q_dict = ci.initialize_charges(coords)
    q_dict, domain_limits = ci.initialize_charges(coords)
    ci.randomize_properties(q_dict, version = 4) # Remove charge from all particles
    ci.randomize_properties(q_dict, version = 5, optarg = 1e13) # Increases mass of all particles to 1000 kg
    timestep = .1
    qsystem = qc.system(q_dict, timestep = timestep, nsteps = int(np.ceil(20/timestep)), domain_limits = domain_limits)
    qsystem.progress_time()
    pltn.animate_pos(qsystem, last_step = -1,  axis_lims =[[-70, 70], [-70,70], [-70,70]])
    
    ############################################################
    # Shell of mixed charges with random velocities
    ############################################################
    # coord_grid, x, y, z = ci.define_grid(xrange = [-10, 10], yrange = [-10, 10], zrange = [-10, 10], res = (1, 1, 1))

    # sphere_coords, ncenter = ci.sphere_shell(2, [0,0,0], x, y, z, res = (1, 1, 1))
    # q_dict, domain_limits = ci.initialize_charges(sphere_coords)
    # ci.randomize_properties(q_dict, version = 7) # Randomizes valocity of half particles to between -10 and 10 m/s
    # ci.randomize_properties(q_dict, version = 1) # half at random flips cahrge
    # ci.randomize_properties(q_dict, version = 3) # hald at random become neutral
    # timestep = .1
    # qsystem = qc.system(q_dict, timestep = timestep, nsteps = int(np.ceil(20/timestep)), domain_limits = domain_limits)
    # qsystem.progress_time()
    # pltn.animate_pos(qsystem, last_step = -1,  axis_lims =[[-70, 70], [-70,70], [-70,70]])
    
    
    ############################################################
    # Shell of positive charges with random velocities in constant E-field
    ############################################################
    # coord_grid, x, y, z = ci.define_grid(xrange = [-10, 10], yrange = [-10, 10], zrange = [-10, 10], res = (1, 1, 1))

    # sphere_coords, ncenter = ci.sphere_shell(2, [0,0,0], x, y, z, res = (1, 1, 1))
    # q_dict, domain_limits = ci.initialize_charges(sphere_coords)
    # ci.randomize_properties(q_dict, version = 7) # Randomizes valocity of half particles to between -10 and 10 m/s
    # timestep = .05
    # qsystem = qc.system(q_dict, timestep = timestep, nsteps = int(np.ceil(6/timestep)), domain_limits = domain_limits, E_amb = [-1e-10,0,0])
    
    # qsystem.progress_time()
    
    # pltn.animate_pos(qsystem, last_step = -1,  axis_lims =[[-70, 70], [-70,70], [-70,70]])
    
    # ############################################################
    # # Shell of uncharged masses charges with random velocities and mass
    # ############################################################
    # coord_grid, x, y, z = ci.define_grid(xrange = [-100, 100], yrange = [-100, 100], zrange = [-100, 100], res = (1, 1, 1))
    # sphere_coords, ncenter = ci.sphere_shell(50, [0,0,0], x, y, z, res = (1, 1, 1), num = 50)
    # q_dict, domain_limits = ci.initialize_charges(sphere_coords)
    # ci.randomize_properties(q_dict, version = 7) # Randomizes valocity of half particles to between -10 and 10 m/s
    # ci.randomize_properties(q_dict, version = 4) # Remove charge from all particles
    # # ci.randomize_properties(q_dict, version = 5, optarg = 1e13) # Increases mass of all particles to 1000 kg
    # ci.randomize_properties(q_dict, version = 2) # randomizes mass to be between 10e-5 kg and 10e15 kg
    # timestep = .1
    # qsystem = qc.system(q_dict, timestep = timestep, nsteps = int(np.ceil(20/timestep)), domain_limits = domain_limits)
    # qsystem.progress_time()
    # pltn.animate_pos(qsystem, last_step = -1,  axis_lims =[[-70, 70], [-70,70], [-70,70]])

    # ############################################################
    # # Shell of uncharged masses charges with random velocities and fixed mass
    # ############################################################
    # coord_grid, x, y, z = ci.define_grid(xrange = [-100, 100], yrange = [-100, 100], zrange = [-100, 100], res = (1, 1, 1))
    # sphere_coords, ncenter = ci.sphere_shell(50, [0,0,0], x, y, z, res = (1, 1, 1), num = 50)
    # q_dict, domain_limits = ci.initialize_charges(sphere_coords)
    # ci.randomize_properties(q_dict, version = 7) # Randomizes valocity of half particles to between -10 and 10 m/s
    # ci.randomize_properties(q_dict, version = 4) # Remove charge from all particles
    # ci.randomize_properties(q_dict, version = 5, optarg = 1e12) # Increases mass of all particles to 1000 kg
    # # ci.randomize_properties(q_dict, version = 2) # randomizes mass to be between 10e-5 kg and 10e15 kg
    # timestep = .1
    # qsystem = qc.system(q_dict, timestep = timestep, nsteps = int(np.ceil(20/timestep)), domain_limits = domain_limits)
    # qsystem.progress_time()
    # pltn.animate_pos(qsystem, last_step = -1,  axis_lims =[[-70, 70], [-70,70], [-70,70]])



    # timestep = .1
    # qsystem = qc.system(q_dict, timestep = timestep, nsteps = int(np.ceil(20/timestep)), domain_limits = domain_limits)
    # qsystem.progress_time()
    # pltn.animate_pos(qsystem, last_step = -1,  axis_lims =[[-70, 70], [-70,70], [-70,70]])