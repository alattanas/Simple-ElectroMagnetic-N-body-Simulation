# -*- coding: utf-8 -*-
"""
Created on Wed Jun 29 00:01:33 2022

@author: Alex
"""

import numpy as np

def calc_lorentz(pos1, charge1, mass1, pos2, charge2, mass2):
    
    # Currently ignoring magnetic fields so F = qE
    G = 6.674 * 10 ** (-11)
    
    eps0 = 8.854 * 10 ** (-12)
    k = 1/(4 * np.pi * eps0)
    
    radius = np.sqrt((pos1[1] - pos2[1]) ** 2 +
                     (pos1[2] - pos2[2]) ** 2 +
                     (pos1[3] - pos2[3]) ** 2)
    
    r = [pos1[0] - pos2[0], pos1[1] - pos2[1], pos1[2] - pos2[2]]
    rmag = np.sqrt(r[0] ** 2 + r[1] ** 2 + r[2] ** 2)
    
    # E2 = [k * charge2 * rcomp / (rmag ** 3) for rcomp in r]
    Fe = [k * charge1 * charge2 * rcomp / (rmag ** 3) for rcomp in r]

    # Add gravity because why not
    Fg = [-1 * G * mass1 * mass2 * rcomp / (rmag ** 3) for rcomp in r]
    
    Feg = [Fe[ind] + Fg[ind] for ind in range(3)]
    
    return Feg
    