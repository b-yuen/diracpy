#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 26 14:37:59 2022

@author: benjaminyuen
"""

# Add diracpy parent directory to system path for testing
import sys
sys.path.insert(0, '../src/')

import diracpy as dp
import numpy as np

##-----------------------------------------------------
## FOR FOCK SUBSPACE TESTS BELOW FIRST INSTANTIATE THESE SUBSPACES
##-----------------------------------------------------
#pump = dp.floquet_subspace(index = 0)
#cav = dp.fock_subspace(index = 1)
#atom = dp.two_level_subspace(index = 2)

##-----------------------------------------------------
## test qop.__str__()
#psi0 = dp.ket([[0,1,'g'],[0,0,'e']],[0, np.sqrt(2)])
#psi0_str = str(psi0)
#print(psi0_str)

##-----------------------------------------------------
## test conjugate action of operator product
#psi0 = dp.ket([0,0,'g'])
#op = cav.adag * atom.sigma_minus
#psi1 = op.conj() * psi0
#print('psi0 = ', end='')
#psi0.print()
#print('op.conj() * psi0  = atom.sigma_plus * cav.a * psi0 = ', end='')
#psi1.print()
#print('psi0.conj * op = psi0.conj * cav.adag * atom.sigma_minus = ', end='')
#psi2 = psi0.conj() * op
#psi2.print()

##-----------------------------------------------------
## test combining operators from different subspaces when this gives an empty ket
#psi0 = dp.ket([0,1,'g'])
#psi1 = cav.a * atom.sigma_minus * psi0
#print('psi0 = ', end='')
#psi0.print()
#print('cav.a * atom.sigma_minus * psi0 = ', end='')
#psi1.print()

##-----------------------------------------------------
## test combining operators from different subspaces
#psi0 = dp.ket([0,1,'g'])
#psi1 = cav.a * atom.sigma_plus * psi0
#print('psi0 = ', end='')
#psi0.print()
#print('cav.a * atom.sigma_minus * psi0 = ', end='')
#psi1.print()

##-----------------------------------------------------
## test action on empty kets
##psi0 = dp.ket([0,0,'g'])
#psi0 = dp.ket()
##psi1 = atom.sigma_minus * atom.sigma_plus * psi0
#psi1 = atom.sigma_minus * psi0
##psi1 = cav.a * psi0
#print('psi0 = ', end='')
#psi0.print()
#print('psi1 = ', end='')
#psi1.print()

##-----------------------------------------------------
## test pump.n operator matrix elements
#n=-5
#psi0 = dp.ket([0,n,'g'])
#psi1 = dp.bra([0,n,'g'])
#print('psi0 = ', end='')
#psi0.print()
#print('psi1 = ', end='')
#psi1.print()
#print('psi1 * pump.n * psi0 = ', end='')
#print(psi1 * pump.n * psi0)

##-----------------------------------------------------
## test pump.adag (creation) operator on bras
#psi0 = dp.bra([0,1,'g'], 1+2j)
#psi1 =  psi0 * pump.adag
#print('psi0 = ', end='')
#psi0.print()
#print('psi0 * pump.adag = ', end='')
#psi1.print()

##-----------------------------------------------------
## test pump.adag (creation) operator on kets
#psi0 = dp.ket([1,0,'g'])
#psi1 = pump.adag * psi0
#print('psi0 = ', end='')
#psi0.print()
#print('pump.adag * psi0 = ', end='')
#psi1.print()

##-----------------------------------------------------
## test pump.a (annihilation) operator on bras
#psi0 = dp.bra([0,-1,'g'], 1+2j)
#psi1 =  psi0 * pump.a
#print('psi0 = ', end='')
#psi0.print()
#print('psi0 * pump.a = ', end='')
#psi1.print()

##-----------------------------------------------------
## test pump.a (annihilation) operator on kets
#psi0 = dp.ket([2,0-,'g'])
#psi1 = pump.a * psi0
#print('psi0 = ', end='')
#psi0.print()
#print('pump.a * psi0 = ', end='')
#psi1.print()

##-----------------------------------------------------
## test cav.n operator matrix elements
#n=10
#psi0 = dp.ket([n,0,'g'])
#psi1 = dp.bra([n,0,'g'])
#print('psi0 = ', end='')
#psi0.print()
#print('psi1 = ', end='')
#psi1.print()
#print('psi1 * cav.n * psi0 = ', end='')
#print(psi1 * cav.n * psi0)

##-----------------------------------------------------
## test cav.adag (creation) operator on bras
#psi0 = dp.bra([0,0,'g'], 1+2j)
#psi1 =  psi0 * cav.adag
#print('psi0 = ', end='')
#psi0.print()
#print('psi0 * cav.adag = ', end='')
#psi1.print()

##-----------------------------------------------------
## test cav.adag (creation) operator on kets
#psi0 = dp.ket([1,0,'g'])
#psi1 = cav.adag * psi0
#print('psi0 = ', end='')
#psi0.print()
#print('cav.adag * psi0 = ', end='')
#psi1.print()

##-----------------------------------------------------
## test cav.a (annihilation) operator on bras
#psi0 = dp.bra([0,0,'g'], 1+2j)
#psi1 =  psi0 * cav.a
#print('psi0 = ', end='')
#psi0.print()
#print('psi0 * cav.a = ', end='')
#psi1.print()

##-----------------------------------------------------
## test cav.a (annihilation) operator on kets
#psi0 = dp.ket([2,1,'g'])
#psi1 = cav.a * psi0
#print('psi0 = ', end='')
#psi0.print()
#print('cav.a * psi0 = ', end='')
#psi1.print()

##-----------------------------------------------------
## FOR THREE LEVEL SYSTEM TESTS BELOW FIRST INSTANTIATE THESE SUBSPACES
##-----------------------------------------------------
#pump = dp.floquet_subspace(index = 0)
#cav = dp.fock_subspace(index = 1)
#atom = dp.three_level_subspace(index = 2)

##-----------------------------------------------------
## test three level system sigma_z operator
#psi0 = dp.ket([0,0,'g1'],1)
#print('psi0 = ', end='')
#psi0.print()
#print('psi0.conj() * atom.sigma_z * psi0 = ', end='')
#print(psi0.conj() * atom.sigma_z * psi0)

##-----------------------------------------------------
## test three level system sigma_plus1 and sigma_minus1 operator matrix elements
#psi0 = dp.ket([0,0,'e'])
#psi1 = dp.bra([0,0,'g1'])
#print('psi0 = ', end='')
#psi0.print()
#print('psi1 = ', end='')
#psi1.print()
#print('psi1 * atom.sigma_minus1 * psi0 = ', end='')
#print(psi1 * atom.sigma_minus1 * psi0)
#print('psi1 * atom.sigma_plus1 * psi0 = ', end='')
#print(psi1 * atom.sigma_plus1 * psi0)

##-----------------------------------------------------
## test three level system sigma_plus1 and sigma_minus1 operators
#psi0 = dp.ket([0,0,'e'])
#psi1 = atom.sigma_plus1 * psi0
#psi1.print()
#psi1 = atom.sigma_minus1 * psi0
#psi1.print()


#-----------------------------------------------------
# FOR ALL TESTS BELOW FIRST INSTANTIATE THESE SUBSPACES
#-----------------------------------------------------
pump = dp.floquet_subspace(index = 0)
cav = dp.fock_subspace(index = 1)
atom = dp.two_level_subspace(index = 2)

#-----------------------------------------------------
# test scalar multiplication of operators
psi0 = dp.ket([0,2,0])
new_op = np.sqrt(2) * cav.a
psi1 = new_op * psi0
print("psi0 = ", end='')
psi0.print()
print("psi1 = np.sqrt(2) * cav.a * psi0 = ", end='')
psi1.print()


##-----------------------------------------------------
## test of qop * ket non-Hermitian operator
#psi0 = dp.bra([0,1,'g'],1+1j)
#print('psi0')
#psi0.print()
#psi1 = psi0 * atom.sigma_plus
#psi2 = psi0 * atom.sigma_minus
#print('psi0 * atom.sigma_plus')
#psi1.print()
#print('psi0 * atom.sigma_minus')
#psi2.print()

##-----------------------------------------------------
## test of qop * ket non-Hermitian operator
#psi0 = dp.ket([0,1,'g'],1)
#print('psi0')
#psi0.print()
#psi1 = atom.sigma_plus * psi0
#psi2 = atom.sigma_minus * psi0
#print('atom.sigma_plus * psi0')
#psi1.print()
#print('atom.sigma_minus * psi0')
#psi2.print()

#-----------------------------------------------------
## test of bra * operator for Hermitian operator
#psi0 = dp.bra([0,1,'g'],1)
#psi0.print()
#psi1 = psi0 * atom.sigma_z
#print(psi1.type)
#psi1.print()

##-----------------------------------------------------
## test of bra ket overlap when bra and ket are empty (shoudl be zero)
## not that the empty state is different from vacuum
## e.g. annihilation of vacumm state gives an empty state with scalar coefficient zero
## and is clearly different from the vacuum state.
#psi0 = dp.ket()
##psi1 = psi0.conj()
#psi1 = dp.bra()
#overlap = psi1 * psi0
#print(overlap)

##-----------------------------------------------------
## test of bra ket overlap
#psi0 = dp.ket( [0,1,'g'] , 1)
##psi1 = psi0.conj()
#psi1 = dp.bra( [[0,0,'g'], [0,1,'g'] ], [0.2,0.5])
#psi2 = dp.bra( [0,1,'g'], 1)
#print(psi1.type)
#overlap = (psi1 + 3 * psi2) * psi0
#print(overlap)


##-----------------------------------------------------
## test operator addition rule
#psi0 = dp.ket( [0,1,'g'] , 1)
#psi0.print()
#
#psi1 = atom.sigma_z * psi0
#psi1.print()
#
#op_test = atom.sigma_z + atom.sigma_z
#
#psi1 = op_test * psi0
#psi1.print()

##-----------------------------------------------------
## test operator acting on superposition states
#psi0 = dp.ket([ [0,1,'g'] , [0,0,'e'] ] , [np.sqrt(1/2),np.sqrt(1/2)])
#psi0.print()
#
#psi1 = atom.sigma_z * psi0
#psi1.print()

##-----------------------------------------------------
## test addition of and empty bra to a bra
#psi0 = dp.bra([0,1,'g'],1)
#psi1 = dp.bra()
#psi2 = psi0 + 2 * psi1
#print('psi0 = ', end='')
#psi0.print()
#print('psi1 = ', end='')
#psi1.print()
#print('psi0 + psi1 = ', end='')
#psi2.print()

##-----------------------------------------------------
## test addition of and empty ket to a ket
#psi0 = dp.ket([0,1,'g'],1)
#psi1 = dp.ket()
#psi2 = psi0 + 2 * psi1
#print('psi0 = ', end='')
#psi0.print()
#print('psi1 = ', end='')
#psi1.print()
#print('psi0 + psi1 = ', end='')
#psi2.print()

##-----------------------------------------------------
## test addition bras
#psi0 = dp.bra([0,1,'g'],1)
#psi1 = dp.bra([0,0,'e'],1)
#psi2 = psi0 + 2 * psi1
#print('psi0 = ', end='')
#psi0.print()
#print('psi1 = ', end='')
#psi1.print()
#print('psi0 + psi1 = ', end='')
#psi2.print()

##-----------------------------------------------------
## test addition kets
#psi0 = dp.ket([0,1,'g'],1)
#psi1 = dp.ket([0,0,'e'],1)
#psi2 = psi0 + 2 * psi1
#print('psi0 = ', end='')
#psi0.print()
#print('psi1 = ', end='')
#psi1.print()
#print('psi0 + psi1 = ', end='')
#psi2.print()

##-----------------------------------------------------
## test scalar multiplication of kets
#psi0 = dp.ket([0,0,'g'],1)
#psi1 = 2 * psi0
#print('psi0 = ', end='')
#psi0.print()
#print('psi1 = ', end='')
#psi1.print()
#print('0 * psi1 = ', end='')
#(0 * psi1).print()

##-----------------------------------------------------
## test scalar multiplication of bras
#psi0 = dp.bra([0,0,'g'],1)
#psi1 = 2 * psi0
#print('psi0 = ', end='')
#psi0.print()
#print('psi1 = ', end='')
#psi1.print()































