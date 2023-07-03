#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 26 14:37:59 2022

@author: benjaminyuen
"""

# # Add diracpy parent directory to system path for testing
# import sys
# # sys.path.insert(0, '../src/')
# sys.path.insert(0, '/Users/benjaminyuen/Documents/Documents - Benjamin’s MacBook Pro/Ben/physics/Birmingham/QuantumCode/diracpy_project/src')

import diracpy as dp
import numpy as np
import matplotlib.pyplot as plt

#subspaces
pump = dp.floquet_subspace(index = 0)
cav = dp.fock_subspace(index = 1)
atom = dp.three_level_subspace(index = 2)


#parameters 
delta_c = 20.
g = 1.76
delta_p = delta_c + 1 + 0.99 * g**2/10.5
rabi = 0.1 # 0.001 for adiabatic passage
kappa = 0.006
effrabi = 0.5 * rabi * g / ((delta_c + delta_p) / 2)

# pump pulse function
pulse_T = 0.5 * 2 * np.pi / effrabi
def pulse_f(t):
#        T = 2 * 2 * np.pi / self.effrabi
        if t<pulse_T/2:
#            val = 0.5 * np.sin(np.pi * t / self.T) ** 2
            val = np.sin(2 * np.pi * t / pulse_T) ** 2
        else:
            val = 0
        return val

#Quantum system and full hamiltonian operator
ham0 = delta_p * pump.n + delta_c * cav.n + atom.sigma_z
ham_v_pump = 0.5 * rabi * (pump.a * atom.sigma_plus1 + pump.adag * atom.sigma_minus1
                           + pump.a * atom.sigma_plus2 + pump.adag * atom.sigma_minus2)
ham_v_cav = g * (cav.a * atom.sigma_plus1 + cav.adag * atom.sigma_minus1
                 + cav.a * atom.sigma_plus2 + cav.adag * atom.sigma_minus2)
#ham = ham0 + ham_v_pump + ham_v_cav
initialbasis = [(0,0,'g1'),(-1,0,'e'),(-1,1,'g2')]
# initialbasis = [dp.ket(state) for state in initialbasis] # initial basis can be tuples of kets
# static qunatum system:
#testsys = dp.qsys(ham, initialstates = initialbasis, n_int = 0, jump_ops = [np.sqrt(kappa) * cav.a])
# Dynamic pumping
testsys = dp.qsys_t(ham0 + ham_v_cav, ham_v_pump, pulse_f,
                    initialstates = initialbasis, n_int = 0, jump_ops = [np.sqrt(kappa) * cav.a])

##Quantum system and effective hamiltonian operator
#ham0 = delta_p * pump.n + delta_c * cav.n + atom.sigma_z
#ham_v_eff = effrabi * (cav.adag * pump.a * atom.sigma_minus2 * atom.sigma_plus1
#                              + pump.adag * cav.a * atom.sigma_minus1 * atom.sigma_plus2)
#initialbasis = [(0,0,'g1'),(-1,0,'e'),(-1,1,'g2')]
#testsys = dp.qsys_t(ham0, ham_v_eff, pulse_f,
#                    initialstates = initialbasis, n_int = 0, jump_ops = [np.sqrt(kappa) * cav.a])

# initial conditions
#psi0 = (dp.ket((0,0,'g1')) + dp.ket((-1,1,'g2'))) * (1 / np.sqrt(2))
psi0 = dp.ket((0,0,'g1'))
rho0_op = psi0 * psi0.conj()
rho0 = testsys.matrix(rho0_op)
#testsys.matprint(rho0)

#Damping terms
#kappa = 0.005
#testsys.jump_operator(kappa * cav.a)

# Time evoluton
#effrabi = g * rabi / (delta_c + delta_p)
cycles = 0.5
Tmax = cycles * 2 * np.pi / effrabi
times = np.linspace(0, Tmax, 500)
#soln = dp.unitaryevolution(rho0, times, testsys)
soln = dp.lindbladint(rho0, times, testsys)
soln.solve()

#plot solution
fig, ax = plt.subplots(testsys.dim+1, 1, figsize = (4,8))
for i in range(testsys.dim):
    ax[i+1].plot(times, np.real(soln.soln[:, i, i]))
#    ax[i+1].set_ylim(-0.05, 1.05)
ax[0].plot(times, [pulse_f(t) for t in times])
    































