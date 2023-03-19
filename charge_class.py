# -*- coding: utf-8 -*-
"""
    Defines the class of charged particles
"""
import numpy as np
from copy import copy, deepcopy

class system:
    def __init__(self, q_dict, timestep = 1, nsteps = 10, domain_limits = [[-20, 20], [-20, 20], [-20, 20]], scale_limits=1, E_amb = [0,0,0]):
        self.timestep = timestep
        self.nsteps = nsteps
        self.q_dict = q_dict
        self.timedev = {}
        self.timealias = {}
        self.domain_limits = domain_limits
        self.scale_limits = scale_limits
        
        if E_amb:
            self.E_amb = E_amb
        else:
            self.E_amb = [0,0,0]

        
    def progress_time(self):
        # minx = self.q_dict
        tkey = 0
        self.timedev.update({tkey: self.q_dict})
        self.timealias.update({tkey: tkey * self.timestep})
        tkey += 1
        print(f"nsteps = {self.nsteps}")
        for tind in range(self.nsteps):
            self.timealias.update({tkey: tkey * self.timestep})
            # print(f"tkey = {tkey}")
            # print(f" timedev keys: {self.timedev.keys()}")
            tempqdict = deepcopy(self.timedev[tkey - 1])
            # get forces for each charge particle adn find next pos and vel
            self._develop_forces(tempqdict)
            
            # update positions and velocities after all forces determined
            self._resolve_movements(tempqdict)
            
            self.timedev.update({tkey: tempqdict})
            
            tkey += 1
            
            # Check for updates to domain limits
            if self.scale_limits == 1:
                for qkey in tempqdict:
                    coord = tempqdict[qkey].pos
                    if coord[0] < self.domain_limits[0][0]:
                        self.domain_limits[0][0] = coord[0]
                    elif coord[0] > self.domain_limits[0][1]:
                        self.domain_limits[0][1] = coord[0]
                        
                    if coord[1]< self.domain_limits[1][0]:
                        self.domain_limits[1][0] = coord[1]
                    elif coord[1] > self.domain_limits[1][1]:
                        self.domain_limits[1][1] = coord[1]
                        
                    if coord[2]< self.domain_limits[2][0]:
                        self.domain_limits[2][0] = coord[2]
                    elif coord[2] > self.domain_limits[2][1]:
                        self.domain_limits[2][1] = coord[2]
            
            
    def _resolve_movements(self, timeq_dict):
        for qkey in timeq_dict:
            timeq_dict[qkey]._update_loc()
        return
        
    # Calculates all forces acting on each particle
    def _develop_forces(self, tempqdict):
        for qkey in tempqdict:
            q1 = tempqdict[qkey]
            Ftot = [0, 0, 0]
            Ftot[0] += self.E_amb[0] * q1.charge
            Ftot[1] += self.E_amb[1] * q1.charge
            Ftot[2] += self.E_amb[2] * q1.charge
            # Accumulating forces from every other particle and summing them
            for sourcekey in tempqdict:
                
                if sourcekey != qkey:
                    q2 = tempqdict[sourcekey]
                    f_on_q1 = self._calc_lorentz(q1, q2)
                    Ftot[0] += f_on_q1[0]
                    Ftot[1] += f_on_q1[1]
                    Ftot[2] += f_on_q1[2]
                    
                else:
                    continue
            q1.Ftot = Ftot
            # print(f" particle {q1.q_id} has tot forces {Ftot}")
            # Calculates next location and velocities of particle for end of step
            q1._force_to_loc(self.timestep)
            
        return
            
                    
                    

    # Calculates electrostatic and gravitational forces acting on q1
    def _calc_lorentz(self, q1, q2):
        # Currently ignoring magnetic fields so loretnz force is F = qE
        G = 6.674 * 10 ** (-11)
        
        eps0 = 8.854 * 10 ** (-12)
        k = 1/(4 * np.pi * eps0)
        
        pos1 = q1.pos
        charge1 = q1.charge
        mass1 = q1.mass
        
        pos2 = q2.pos
        charge2 = q2.charge
        mass2 = q2.mass
        
        
        r = [pos1[0] - pos2[0], pos1[1] - pos2[1], pos1[2] - pos2[2]]
        rmag = np.sqrt(r[0] ** 2 + r[1] ** 2 + r[2] ** 2)
        
        # E2 = [k * charge2 * rcomp / (rmag ** 3) for rcomp in r]
        Fe = [k * charge1 * charge2 * rcomp / (rmag ** 3) for rcomp in r]
        
        # Add gravity because why not
        Fg = [-1 * G * mass1 * mass2 * rcomp / (rmag ** 3) for rcomp in r]
        for comp in range(3):
            if np.isnan(Fg[comp]):
                Fg[comp] = 0
            if np.isnan(Fe[comp]):
                Fe[comp] = 0
                
        Feg = [Fe[ind] + Fg[ind] for ind in range(3)]
        
        return Feg
                    

class point_charge:
    
    def __init__(self, q_id, mass, charge,  init_loc, init_vel = (0,0,0)):
        self.pos = list(init_loc)
        self.vel = list(init_vel)
        self.mass = mass
        self.charge = charge
        self.q_id = q_id
        if charge == 0:
            self.color = "grey"
        elif charge > 0:
            self.color = "red"
        else:
            self.color = "blue"
    # Numerical equations of motion using total forces on particles to find next location
    def _force_to_loc(self, timestep):
        c = 2.99e8 #m/s speed of light
        rel_thresh = .7 * c
        maxv = .95 * c
        
        self.vel_next = [0,0,0]
        self.pos_next = [0,0,0]
        
        # Checking for relativistic effects assuming relavaent at .7 C
        ax = self.Ftot[0] / self.mass
        ay = self.Ftot[1] / self.mass
        az = self.Ftot[2] / self.mass
        
        vx_next = min(self.vel[0] + ax * timestep, maxv)
        vy_next = min(self.vel[1] + ay * timestep, maxv)
        vz_next = min(self.vel[2] + az * timestep, maxv)
        
        x_next = self.pos[0] + self.vel[0] * timestep + ax * (timestep ** 2) / 2
        y_next = self.pos[1] + self.vel[1] * timestep + ay * (timestep ** 2) / 2
        z_next = self.pos[2] + self.vel[2] * timestep + az * (timestep ** 2) / 2
        
        self.vel_next = [vx_next, vy_next, vz_next]
        self.pos_next = [x_next, y_next, z_next]
    
    # Finishes updating new positions and velocities at end of step
    def _update_loc(self):
        self.vel = self.vel_next
        self.pos = self.pos_next
        