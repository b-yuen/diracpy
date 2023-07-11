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
import time

class qsys:
    def __init__(self, hamiltonian_operator=None, H0terms=None, Vterms=None,
                 jump_ops=None, n_int=0, initialstates=None, **kwargs):
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
        
        
        self.n_int = n_int
        self._input_hamiltonian(hamiltonian_operator, H0terms, Vterms)
        self.jump_ops = self._input_jump_ops(jump_ops)
        self.build_sys(initialstates)
        self.make_lindblads()
        
        # if not isinstance(hamiltonian_operator, qop):
        #     raise TypeError("qsys first positional argument should be of type diracpy.states_operators.qop")
        # self.ham_op = hamiltonian_operator
        # try:
        #     self.jump_ops = list(kwargs.get('jump_ops', []))
        # except TypeError:
        #     self.jump_ops = [ kwargs['jump_ops'] ]
        # if 'initialstates' in kwargs:
        #     print('building hamiltonian...')
        #     self.build(kwargs['initialstates'])
        #     print('built')
        # self.make_lindblads()
        
    # returns empty list if no jump_op_in given, or a list of jump operators otherwise
    def _input_jump_ops(self, jump_op_in):
        if jump_op_in == None:
            return []
        else:
            try:
                return list(jump_op_in)
            except:
                return [jump_op_in]
        
    # Check at least on of input arguments hamiltonian_operator, H0terms or Vterms are populated.
    # If populated check hamiltonian is a qop or H0terms and/or V0terms are list or array of qops.
    def _input_hamiltonian(self, hamiltonian_operator, H0terms, Vterms):
        if isinstance(hamiltonian_operator, qop):
            self.ham_op = hamiltonian_operator
            self.ham_array_flag = False
        else:
            self.H0terms = []
            self.Vterms = []
            if H0terms == None:
                pass
            elif isinstance(H0terms[0], qop):
                self.H0terms = H0terms
                self.ham_array_flag = True
            else:
                raise TypeError("H0terms should be a list or array of type diracpy.states_operators.qop")
            if Vterms == None:
                pass
            elif isinstance(Vterms[0], qop):
                self.Vterms = Vterms
                self.ham_array_flag = True
            else:
                raise TypeError("Vterms should be a list or array of type diracpy.states_operators.qop")
        try:
            self.ham_array_flag
        except:
            NameError("No hamiltonian has been specified. At least one of hamiltonian_operator, H0terms or Vterms must be given.")
            
    def build_sys(self, initial_states):
        # print(initial_states)
        if initial_states == None:
            raise NameError("System not built, no initial basis states given.")
        else:
            self._build(initial_states)
        
    def _build(self, initialstates):
        # build basis of coherent coupled states to order n_int
        # self.n_int = n_int
        # self.basis = [ket(state) for state in initialstates]
        # for _ in range(n_int):
        #     new_states = []
        #     for state in self.basis:
        #         state_out = self.ham_op * state
        #         state_out_components = [ket(component_state) for component_state in state_out.vec.keys()]
        #         [new_states.append(new_state) for new_state in state_out_components if new_state not in self.basis]
        #     self.basis += new_states
        # print("building...")
        print("building sys basis...")
        t1 = time.time()
        
        self.basis = self._build_coherent_sys(initialstates)
        
        # add states these decay to
        if len(self.jump_ops) == 0:
            basis_incomplete = False
        else:
            basis_incomplete = True
        while basis_incomplete:
            print('...appending decay subspaces...')
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
        
        t2 = time.time()
        print("...system basis built in {} seconds".format(t2-t1))

        self.adjoint_basis = [bra(state) for state in self.basis]
        self.dim = len(self.basis)
        print("defining hmatrix...")
        t3 = time.time()
        self.hmatrix = self.make_hmatrix()
        t4 = time.time()
        print("...hmatrix evaluated in {} seconds".format(t4-t3))
        
    # builds the basis states coherently coupled to initialstates via hamiltonian  
    def _build_coherent_sys(self, initialstates):
        # build sys when hamiltonian given via arrays H0terms and or Vterms
        if self.ham_array_flag == True:
            builder = self._build_from_ham_array
        # build sys when hamiltonian given via single hamiltonian_operator
        else:
            builder = self._build_from_single_ham
        return builder(initialstates)
        
    # This method more efficiently builds basis for extensive operators when given as H0terms and Vterms
    # Builds basis states coherently coupled to initial states via operators in Vterms to order n_int
    def _build_from_ham_array(self, initialstates):
        basis = [ket(state) for state in initialstates]
        for _ in range(self.n_int):
            for Vop in self.Vterms:
                new_states = []
                for state in basis:
                    state_out = Vop * state
                    state_out_components = [ket(component_state) for component_state in state_out.vec.keys()]
                    [new_states.append(new_state) for new_state in state_out_components if new_state not in basis]
                basis += new_states
        return basis
            
    # Build basis coherently coupled to initial states via hamiltonian_operator to order n_int.
    def _build_from_single_ham(self, initialstates):
        basis = [ket(state) for state in initialstates]
        for _ in range(self.n_int):
            new_states = []
            for state in basis:
                state_out = self.ham_op * state
                state_out_components = [ket(component_state) for component_state in state_out.vec.keys()]
                [new_states.append(new_state) for new_state in state_out_components if new_state not in basis]
            basis += new_states
        return basis
    
    def make_hmatrix(self):
        if self.ham_array_flag == True:
            matrix_evaluator = self._hmatrix_from_ham_array
        else:
            matrix_evaluator = self._hmatrix_from_ham_op
        return matrix_evaluator()
            
    def _hmatrix_from_ham_array(self):
        hmatrix = np.zeros([self.dim, self.dim], complex)
        for term_op in self.H0terms:
            term_matrix = self.matrix(term_op)
            hmatrix += term_matrix
        for term_op in self.Vterms:
            term_matrix = self.matrix(term_op)
            hmatrix += term_matrix
        return hmatrix
    
    def _hmatrix_from_ham_op(self):
        hmatrix = self.matrix(self.ham_op)
        return hmatrix
            
    def matrix(self, operator):
        if not isinstance(operator, qop):
            raise TypeError("first positional argument should be of type diracpy.states_operators.qop")
        matrix_out = np.zeros([self.dim, self.dim], complex)
        for i, basis_bra in enumerate(self.adjoint_basis):
            _brapsi = basis_bra * operator
            for j, basis_ket in enumerate(self.basis):
                # matrix_out[i,j] = basis_bra * operator * basis_ket
                matrix_out[i,j] = _brapsi * basis_ket
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
            
        
        
        
        
        
        
        

        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
    