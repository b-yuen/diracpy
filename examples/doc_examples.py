#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Examples used in documentation
"""

# Created on Tue Jun 11 16:28:56 2024
# @author: benjaminyuen

import os
import sys
sys.path.insert(0, os.path.abspath('/Users/benjaminyuen/Documents/Documents - Benjamin’s MacBook Pro/Ben/physics/Birmingham/QuantumCode/diracpy_project/src'))

import diracpy as dp
import numpy as np
import matplotlib.pyplot as plt


#%% unitary evolution

#define system
atom = dp.two_level_subspace(index=0)
cav = dp.fock_subspace(index=1)
Delta, g = 0, np.pi
H_0 = Delta * atom.sigma_z
V = g * (cav.a * atom.sigma_plus + cav.adag * atom.sigma_minus)
H = H_0 + V
psi0 = dp.ket(['e',0])
system = dp.qsys(H, initialstates=[psi0], n_int=2)

# construct unitary solver
times = np.linspace(0,2,100)
rho0 = psi0 * psi0.conj()
usolver = dp.unitaryevolution(rho0, times, system)
# solve
usolver.solve()

# plot soln
fig, ax = plt.subplots(1,1)
labels = [state for state in system.basis]
for i in range(system.dim):
    ax.plot(times, 
            np.real(usolver.soln[:,i,i]),
            label=labels[i])
plt.legend(loc=1)
plt.show()













