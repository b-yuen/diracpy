#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 15 21:31:18 2022

@author: benjaminyuen
"""

from diracpy.states_operators import ket
from diracpy.states_operators import bra
from diracpy.states_operators import qop
import numpy as np

class qsys:
    def __init__(self, hamiltonian_operator, **kwargs):
        # hamiltonian operator should be of type dp.qop
        
        # key word arguments are: 
        #
        # 'initialstates', a list of initial states from which to form a basis
        # these states are either input as an interable or a dp.ket object
        #
        # 'n_int' is an integer for the order of interactions at which additional
        # basis states should be inluded in basis. I.e. those connected to
        # state in the 'intialstates' list by interactions up to 
        # the 'hamiltonian_operator' to the power 'n_int'.
        #
        # 'jump_ops' are dp.qop objects or a list of such objects which decribe
        # quantum jump operators for the system if it is an open quantum system
        # quantum jump operators are also known as the Lindblad lowering operators
        # decay coefficients should be built into these operators, e.g, for
        # a system which decays at rate kappa via annihilation operator fock_subspace.a
        # the jump operator should be np.sqrt(kappa) * fock_subspace.a
        
        if not isinstance(hamiltonian_operator, qop):
            raise TypeError("qsys first positional argument should be of type diracpy.states_operators.qop")
        self.ham_op = hamiltonian_operator
        try:
            self.jump_ops = list(kwargs.get('jump_ops', []))
        except TypeError:
            self.jump_ops = [ kwargs['jump_ops'] ]
        if 'initialstates' in kwargs:
#            print('building hamiltonian...')
            n_int = kwargs.get('n_int', 0)
            self.build(kwargs['initialstates'], n_int)
        self.make_lindblads()
        
    def build(self, initialstates, n_int):
        # build basis of coherent coupled states to order n_int
        self.n_int = n_int
        self.basis = [ket(state) for state in initialstates]
        for _ in range(n_int):
            new_states = []
            for state in self.basis:
                state_out = self.ham_op * state
                state_out_components = [ket(component_state) for component_state in state_out.vec.keys()]
                [new_states.append(new_state) for new_state in state_out_components if new_state not in self.basis]
            self.basis += new_states
        # add states these decay to
        basis_incomplete = True
        while basis_incomplete:
            new_states = []
            for state in self.basis:
                for jump_op in self.jump_ops:
                    state_out = jump_op * state
                    state_out_components = [ket(component_state) for component_state in state_out.vec.keys()]
                    [new_states.append(new_state) for new_state in state_out_components if new_state not in self.basis]
            basis_incomplete = bool( len(new_states)>0 )
            unique_new_states = []
            for state in new_states:
                if state not in unique_new_states:
                    unique_new_states.append(state)
            self.basis += unique_new_states
                
        
#        while len(decayedbasis) > 0:
#
#                newbasisstates = []
#                for state in decayedbasis:
#                    newbasisstates += self.decaycoupledstates(state)
#
#                decayedbasis = []
#                [decayedbasis.append(state) for state in newbasisstates if state not in decayedbasis];
#                basis += decayedbasis
        
        self.adjoint_basis = [bra(state) for state in self.basis]
        self.dim = len(self.basis)
        self.hmatrix = self.matrix(self.ham_op)
            
    def matrix(self, operator):
        if not isinstance(operator, qop):
            raise TypeError("first positional argument should be of type diracpy.states_operators.qop")
        matrix_out = np.zeros([self.dim, self.dim], complex)
        for i, basis_bra in enumerate(self.adjoint_basis):
            for j, basis_ket in enumerate(self.basis):
                matrix_out[i,j] = basis_bra * operator * basis_ket
        return matrix_out
    
    def print_basis(self):
        [state.print() for state in self.basis]
        
    def print_ham(self):
        self.matprint(self.hmatrix)
                
    def matprint(self, mat, fmt="g"):
        col_maxes = [max([len(("{:"+fmt+"}").format(x)) for x in col]) for col in mat.T]
        for x in mat:
            for i, y in enumerate(x):
                print(("{:"+str(col_maxes[i])+fmt+"}").format(y), end="  ")
            print("")
                
    def ham(self, t):
        return self.hmatrix
    
    def add_jump_op(self, jump_qop):
        if not isinstance(jump_qop, qop):
            raise TypeError("first positional argument should be of type diracpy.states_operators.qop")
        # jump_qop should be the lowering operator associated with the quantum jump
        self.jump_ops.append(jump_qop)
        
    def make_lindblads(self):
        self.lindbladgamma = {}
        self.lindbladraising = {}
        self.lindbladlowering = {}
        for index, jump_qop in enumerate(self.jump_ops):
            self.lindbladgamma[index] = 1
            self.lindbladlowering[index] = self.matrix(jump_qop)
            self.lindbladraising[index] = self.matrix(jump_qop.conj())
            
class qsys_t(qsys):
    def __init__(self, static_operator, dynamic_operators, dynamic_coefficients, **kwargs):
        # 'static operator' should be a dp.qop object describing the static part of the Hamiltonian
        #
        # 'dynamic_operators' should be a dp.qop object of list of objects describing the
        # dynamic part of the Hamiltonian.
        #
        # 'dynamic_coefficients' is a scalar value function or list of functions
        # with single scalar input (time parameter). These should correspond to the
        # time dependent coefficients of the dynamic operators respectively.
        # 
        # for key word arguments see comments for parent class qsys.
        self.static_operator = static_operator
        self.dynamic_operators = dynamic_operators
        self.dynamic_coefficients = dynamic_coefficients
        self.static_hamiltonian = self.static_operator
        try:
            for op in self.dynamic_operators:
                self.static_hamiltonian = self.static_hamiltonian + op
        except TypeError:
            self.static_hamiltonian = self.static_hamiltonian + self.dynamic_operators
        super().__init__(self.static_hamiltonian, **kwargs)
        self.ham_matrices()
        
    def ham_matrices(self):
        self.static_op_matrix = self.matrix(self.static_operator)
        try:
            self.dynamic_op_matrices = [self.matrix(op) for op in self.dynamic_operators]
        except TypeError:
            self.dynamic_op_matrices = [self.matrix( self.dynamic_operators )]
            self.dynamic_coefficients = [self.dynamic_coefficients]
        # Check 'dynamic_coefficients' are correctly specified
        try:
            _testval = 0
            for func in self.dynamic_coefficients:
                _testval += func(0)
        except TypeError:
            raise TypeError("'dynamic_coefficients' should be a scalar function or list of functions of 1 positional argument")
        if len(self.dynamic_coefficients) == len(self.dynamic_op_matrices):
            self.num_dynamic_components = len(self.dynamic_coefficients)
        else:
            raise IndexError("'dynamic_operators' and 'dynamic_coefficients' inputs should be same length")
            
        
    def ham(self, t):
        current_hmatrix = self.static_op_matrix.copy()
        for i in range(self.num_dynamic_components):
            current_hmatrix += self.dynamic_coefficients[i](t) * self.dynamic_op_matrices[i]
        return current_hmatrix
            
        
        
        
        
        
        
        

        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
    