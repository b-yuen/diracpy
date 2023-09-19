#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 24 14:57:04 2023

@author: benjaminyuen

Performance test of diracpy operators
"""

import diracpy as dp
import numpy as np
import timeit


# #%% Test 1, JC interaction terms, two methods.
# # operators
# dipole = dp.two_level_subspace(index = 0)
# field = dp.fock_subspace(index = 1)
# g = np.pi

# #initial state (ket) and final state (bra)
# psi1 = dp.ket(['e',0])
# psi2_star = dp.bra(['g',1])

# composite_op = g * field.adag * dipole.sigma_minus + g * field.a * dipole.sigma_plus
# oplist1 = [g, field.adag, dipole.sigma_minus]
# oplist2 = [g, field.a, dipole.sigma_plus]

# def method1_code():
#     me = psi2_star * composite_op * psi1

# def method2_code():
#     phi1 = psi1
#     for op in oplist1:
#         phi1 = op * phi1
#     phi2 = psi1
#     for op in oplist2:
#         phi2 = op * phi2
#     psi = phi1 + phi2
#     me = psi2_star * psi

# # Number of times to repeat each test
# num_repeats = 1000

# # Measure execution time for Test 1
# method1_time = timeit.timeit(method1_code, number=num_repeats)
# print('Method 1 average execution time:', method1_time / num_repeats)

# # Measure execution time for Test 2
# method2_time = timeit.timeit(method2_code, number=num_repeats)
# print('Method 2 average execution time:', method2_time / num_repeats)

#%% Test 2: same interaction as in JC model but for 10 different modes

dipole = dp.two_level_subspace(index = 0)
field = []
gs = []
nmodes=120
for i in range(nmodes):
    field.append(dp.fock_subspace(index = i+1))
    gs.append((i+1) * np.pi + 1.j)

psi1 = dp.ket(['e'] + [0 for i in range(nmodes)])
psi2_star = dp.bra(['g',1] + [0 for i in range(nmodes-1)])
    
def method1_code():
    composite_op = gs[0] * (field[0].adag * dipole.sigma_minus + field[0].a * dipole.sigma_plus)
    for i in range(1,nmodes):
        composite_op = composite_op +  (gs[i] * 
                                        (field[i].adag * dipole.sigma_minus + 
                                         field[i].a * dipole.sigma_plus))
    me = psi2_star * composite_op * psi1
    
# def method2_code():
#     adag_list = [op.adag for op in field]
#     a_list = [op.a for op in field]
#     me = 0
#     for i in range(nmodes):
#         me += psi2_star * gs[i] * (adag_list[i] * dipole.sigma_minus + a_list[i] * dipole.sigma_plus)  * psi1
#         # me += psi2_star * gs[i] * a_list[i] * dipole.sigma_plus * psi1
        
# def method2_code():
#     me = 0
#     for i in range(nmodes):
#         me += psi2_star * gs[i] * (field[i].adag * dipole.sigma_minus + field[i].a * dipole.sigma_plus)  * psi1

def method2_code():
    adag_list = [op.adag for op in field]
    a_list = [op.a for op in field]
    psia = dipole.sigma_minus * psi1
    psib = dipole.sigma_plus * psi1
    me = 0
    for i in range(nmodes):
        # b1 = psi2_star * adag_list[i]
        # b2 = psi2_star * a_list[i]
        # me += gs[i] * (b1 * psia)
        # me += gs[i] * (b2 * psib)
        k1 = adag_list[i] * psia
        k2 = a_list[i] * psib
        me += gs[i] * (psi2_star * k1)
        me += gs[i] * (psi2_star * k2)
        # me += gs[i] * (psi2_star * adag_list[i] * psia)
        # me += gs[i] * (psi2_star * a_list[i] * psib)
        
# Number of times to repeat each test
num_repeats = 100

# # Measure execution time for Test 1
# method1_time = timeit.timeit(method1_code, number=num_repeats)
# print('Method 1 average execution time:', method1_time / num_repeats)

# Measure execution time for Test 2
method2_time = timeit.timeit(method2_code, number=num_repeats)
print('Method 2 average execution time:', method2_time / num_repeats)
    
        
    




    

