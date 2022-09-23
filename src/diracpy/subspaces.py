#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 15 15:44:22 2022

@author: benjaminyuen
"""

from diracpy.states_operators import ket
#from diracpy.states import bra
from diracpy.states_operators import qop
from numpy import sqrt

class subspace:
    def __init__(self, **kwargs):
        self.index = kwargs.get('index', 0)
        self.basis = []
        
    def modify_bstate(self, bstate, new_label):
        # produces a copy of a basis state tuple with the element at index self.index changed to new_label
        # modifies bstate made by concatenation
        return bstate[:self.index] + (new_label,) + bstate[self.index+1:]
        
        
class two_level_subspace(subspace):
    # need to build index into ket referencing so that it can handle multi-indexed basis states
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.subspacebasis = [('e',),('g',)]
        self.sigma_z_dict = {'e':0.5, 'g':-0.5}
        self.sigma_z = qop(self.sigma_z_action)
        self.sigma_plus = qop(self.sigma_plus_action, self.sigma_minus_action)
        self.sigma_minus = qop(self.sigma_minus_action, self.sigma_plus_action)

    def sigma_z_action(self, ket_in):
        # assumes that the input ket_in is for a single basis state but with arbitrary cnum.
        try:
            basisstate_in = list(ket_in.vec.keys())[0]
            c_in = ket_in.vec[basisstate_in]
            c_out = self.sigma_z_dict[basisstate_in[self.index]] * c_in
            ket_out = ket(basisstate_in, c_out)
        except IndexError:
            ket_out = ket()
        return ket_out
    
    def sigma_minus_action(self, ket_in):
        try:
            basisstate_in = list(ket_in.vec.keys())[0]
            c_in = ket_in.vec[basisstate_in]
            ket_out = ket()
            if basisstate_in[self.index] == 'e':
                basisstate_out = self.modify_bstate(basisstate_in, 'g')
                ket_out.update( basisstate_out, c_in )
        except IndexError:
            ket_out = ket()
        return ket_out
    
    def sigma_plus_action(self, ket_in):
        try:
            basisstate_in = list(ket_in.vec.keys())[0]
            c_in = ket_in.vec[basisstate_in]
            ket_out = ket()
            if basisstate_in[self.index] == 'g':
                basisstate_out = self.modify_bstate(basisstate_in, 'e')
                ket_out.update( basisstate_out, c_in )
        except IndexError:
            ket_out = ket()
        return ket_out
        
class three_level_subspace(subspace):
    # need to build index into ket referencing so that it can handle multi-indexed basis states
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.subspacebasis = [('e',),('g1',),('g2',)]
        self.sigma_z_dict = {'g1':-0.5, 'g2':0.5, 'e':0}
        self.sigma_z = qop(self.sigma_z_action)
        self.sigma_minus1 = qop(self.sigma_minus1_action, self.sigma_plus1_action)
        self.sigma_plus1 = qop(self.sigma_plus1_action, self.sigma_minus1_action)
        self.sigma_minus2 = qop(self.sigma_minus2_action, self.sigma_plus2_action)
        self.sigma_plus2 = qop(self.sigma_plus2_action, self.sigma_minus2_action)
    
    def sigma_minus1_action(self, ket_in):
        try:
            basisstate_in = list(ket_in.vec.keys())[0]
            c_in = ket_in.vec[basisstate_in]
            ket_out = ket()
            if basisstate_in[self.index] == 'e':
                basisstate_out = self.modify_bstate(basisstate_in, 'g1')
                ket_out.update( basisstate_out, c_in )
        except IndexError:
            ket_out = ket()
        return ket_out
    
    def sigma_plus1_action(self, ket_in):
        try:
            basisstate_in = list(ket_in.vec.keys())[0]
            c_in = ket_in.vec[basisstate_in]
            ket_out = ket()
            if basisstate_in[self.index] == 'g1':
                basisstate_out = self.modify_bstate(basisstate_in, 'e')
                ket_out.update( basisstate_out, c_in )
        except IndexError:
            ket_out = ket()
        return ket_out
    
    def sigma_minus2_action(self, ket_in):
        try:
            basisstate_in = list(ket_in.vec.keys())[0]
            c_in = ket_in.vec[basisstate_in]
            ket_out = ket()
            if basisstate_in[self.index] == 'e':
                basisstate_out = self.modify_bstate(basisstate_in, 'g2')
                ket_out.update( basisstate_out, c_in )
        except IndexError:
            ket_out = ket()
        return ket_out
    
    def sigma_plus2_action(self, ket_in):
        try:
            basisstate_in = list(ket_in.vec.keys())[0]
            c_in = ket_in.vec[basisstate_in]
            ket_out = ket()
            if basisstate_in[self.index] == 'g2':
                basisstate_out = self.modify_bstate(basisstate_in, 'e')
                ket_out.update( basisstate_out, c_in )
        except IndexError:
            ket_out = ket()
        return ket_out
    
    def sigma_z_action(self, ket_in, print_output = False):
        try:
            basisstate_in = list(ket_in.vec.keys())[0]
            key = basisstate_in[self.index]
            ket_out = self.sigma_z_dict[key] * ket_in
        except IndexError:
            ket_out = ket()
        return ket_out
        
class fock_subspace(subspace):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
#        self.index = 0
        if 'max_n' in kwargs:
            self.basis = [(n,) for n in range(kwargs['max_n']+1)]
        self.a = qop(self.a_action, self.adag_action)
        self.adag = qop(self.adag_action, self.a_action)
        self.n = qop(self.n_action)
        
    def a_action(self, ket_in):
        try:
            basisstate_in = list(ket_in.vec.keys())[0]
            n_in = basisstate_in[self.index]
            c_in = ket_in.vec[basisstate_in]
            ket_out = ket()
            if n_in > 0:
                c_out = sqrt(n_in) * c_in
                basisstate_out = self.modify_bstate(basisstate_in, n_in - 1)
                ket_out.update( basisstate_out, c_out)
        except IndexError:
            ket_out = ket()
        return ket_out
    
    def adag_action(self, ket_in):
        try:
            basisstate_in = list(ket_in.vec.keys())[0]
            n_in = basisstate_in[self.index]
            c_in = ket_in.vec[basisstate_in]
            ket_out = ket()
            c_out = sqrt(n_in + 1) * c_in
            basisstate_out = self.modify_bstate(basisstate_in, n_in + 1)
            ket_out.update( basisstate_out, c_out)
        except IndexError:
            ket_out = ket()
        return ket_out
    
    def n_action(self, ket_in):
        try:
            basisstate_in = list(ket_in.vec.keys())[0]
            n_in = basisstate_in[self.index]
            ket_out = n_in * ket_in
        except IndexError:
            ket_out = ket()
        return ket_out
    
class floquet_subspace(subspace):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
#        self.index = 0
        try:
            self.basis = [(n,) for n in range(kwargs['n_range'])]
        except KeyError:
            pass
        self.a = qop(self.a_action, self.adag_action)
        self.adag = qop(self.adag_action, self.a_action)
        self.n = qop(self.n_action)
        
    def a_action(self, ket_in):
        try:
            basisstate_in = list(ket_in.vec.keys())[0]
            n_in = basisstate_in[self.index]
            c_in = ket_in.vec[basisstate_in]
            ket_out = ket()
            c_out = c_in
            basisstate_out = self.modify_bstate(basisstate_in, n_in - 1)
            ket_out.update( basisstate_out, c_out)
        except IndexError:
            ket_out = ket()
        return ket_out
    
    def adag_action(self, ket_in):
        try:
            basisstate_in = list(ket_in.vec.keys())[0]
            n_in = basisstate_in[self.index]
            c_in = ket_in.vec[basisstate_in]
            ket_out = ket()
            c_out = c_in
            basisstate_out = self.modify_bstate(basisstate_in, n_in + 1)
            ket_out.update( basisstate_out, c_out)
        except IndexError:
            ket_out = ket()
        return ket_out
    
    def n_action(self, ket_in):
        try:
            basisstate_in = list(ket_in.vec.keys())[0]
            n_in = basisstate_in[self.index]
            ket_out = n_in * ket_in
        except IndexError:
            ket_out = ket()
        return ket_out