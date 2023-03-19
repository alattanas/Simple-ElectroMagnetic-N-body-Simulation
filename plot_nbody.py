# -*- coding: utf-8 -*-
"""
    Plot functions for charges in N-body system 
"""

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from celluloid import Camera

# fig = plt.figure()
# camera = Camera(fig)
# plt.plot(...)
# plt.fancy_stuff()
# camera.snap()
# animation = camera.animate()
# animation.save('animation.mp4')
plt.close("all")

# Animates entries of the timedev dictionary of a system instance
def animate_pos(qsystem, last_step = -1, axis_lims = []):
    print("Plotting Routine Started")
    if last_step == -1:
        last_step = len(qsystem.timedev)
        
    dt = qsystem.timestep

    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    
    # ax.set_title("N-body with dt = {} at time = {} s".format(dt, 0))
    ax.set_xlabel("X position [m]")
    ax.set_ylabel("Y position [m]")
    ax.set_zlabel("Z position [m]")
    # print(qsystem.domain_limits)
    if not axis_lims:
        axis_lims = axis_limits(qsystem.domain_limits)
    # print(axis_lims)
    ax.set_xlim(axis_lims[0])
    ax.set_ylim(axis_lims[1])
    ax.set_zlim(axis_lims[2])

    camera = Camera(fig)
    for tkey in qsystem.timedev:
        # print(f"animating timestep {tkey} when time = {tkey * dt}")
        
        for qkey in qsystem.timedev[tkey]:
            particle = qsystem.timedev[tkey][qkey]
            ax.scatter(particle.pos[0], particle.pos[1], particle.pos[2], color = particle.color)
            # plt.show()
        title_string = str("N-body with dt = {} at time = {: .3f} s".format(dt, round(tkey * dt, 2)))
        # This is the only way to update title with calluloid animations
        ax.text2D(.5, 1, title_string, fontsize = 14, ha='center', transform=ax.transAxes)
        print(f"tkey = {tkey}, time = {tkey * dt}")
        camera.snap()

    animation = camera.animate(interval = 50)
    animation.save('animation.mp4')
    animation.save('animation.gif')
    return camera

    # Calculates plot axis limits to give some relative buffer away from eaxh min and max value plotted
def axis_limits(minmax_coords, buffer = .1, mindist = 10):
    xdist = max(minmax_coords[0][1] - minmax_coords[0][0], mindist)
    ydist = max(minmax_coords[1][1] - minmax_coords[1][0], mindist)
    zdist = max(minmax_coords[2][1] - minmax_coords[2][0], mindist)
    
    xmin = minmax_coords[0][0] - buffer * xdist
    xmax = minmax_coords[0][1] + buffer * xdist
    ymin = minmax_coords[1][0] - buffer * ydist
    ymax = minmax_coords[1][1] + buffer * ydist
    zmin = minmax_coords[2][0] - buffer * zdist
    zmax = minmax_coords[2][1] + buffer * zdist
    
    return [[xmin, xmax], [ymin, ymax], [zmin, zmax]]
    