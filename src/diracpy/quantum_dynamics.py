# -*- coding: utf-8 -*-
# """
# Spyder Editor

# This is a temporary script file.
# """

# This script defines classes which can be used to
# numerically solve quantum dynamics

import numpy as np
import scipy
from scipy.integrate import odeint
from random import random
from multiprocessing import Pool
import time

# The vonneumannint class numerically integrates the vonneumann equation
# using odeint. A Hamiltonian matrix and the initial value of the density
# matrix must be given to instantiate this class, together with a list of
# times at which to calculate the time-dependent density matrix. Thus, the 
# vonneumannint class has no knowledge of the specific system to be solved.
# It is only assumed that this system can be represented by on a finite
# dimensional Hilbert space such that the Hamiltonian and density matrix
# can be represented by a finite dimensional matrix.
# Input variables z0 (initial density matrix)  and hmatrix (the Hamitlonian)
# are specified as n by n complex numpy arrays, and t is a np.array.

# class lindbladint:
    
#     def __init__(self, z0, t, ham_obj):
#         self.z0 = z0
#         self.t = t
#         self.ham = ham_obj
#         self.y0matrix = self.c2rmatrix(self.z0)
#         self.rhmatrix = self.c2rmatrix(-1.j * self.ham.hmatrix)
#         self.rpumpmatrices = tuple([self.c2rmatrix(-1.j * self.ham.pumpmatrices[i]) for i in range(self.ham.numatoms)])
        
#         self.vdim = np.size(self.rhmatrix)
#         self.mshape = np.shape(self.rhmatrix)
# #         self.dim = np.shape(self.ham.hmatrix)[0]
#         self.dim = self.ham.dim
        
#         self.y0 = self.y0matrix.reshape(self.vdim)

class vonneumannint:
    
    def __init__(self, z0, t, ham_obj):
        self.z0 = z0
        self.t = t
        self.ham = ham_obj
        self.y0matrix = self.c2rmatrix(self.z0)
        self.rhmatrix = self.c2rmatrix(-1.j * self.ham.hmatrix)
        if 'pumpmatrices' in self.ham.__dict__:
            self.rpumpmatrices = tuple([self.c2rmatrix(-1.j * self.ham.pumpmatrices[i]) for i in range(self.ham.numatoms)])
            for k in range(self.ham.numatoms):
                self.rhmatrix += self.ham.rabis[k] * self.rpumpmatrices[k]
        
        self.vdim = np.size(self.rhmatrix)
        self.mshape = np.shape(self.rhmatrix)
#         self.dim = np.shape(hmatrix)[0]
        self.dim = self.ham.dim
        
        self.y0 = self.y0matrix.reshape(self.vdim)
        
        self.realsoln = odeint(self.derivs, self.y0, self.t, args = (self.rhmatrix,))
#         self.realsoln = odeint(self.derivs, self.y0, self.t, 
#                                args = (self.rhmatrix, self.rpumpmatrices, 
#                                        self.ham.rabis, self.ham.numatoms))
        self.soln = self.rvec2cmatrixsoln()
    
    def c2rmatrix(self, matrix):
        dim = np.shape(matrix)[0]
#         Should be a square matrix. May want to put in some error handling
        realm = np.zeros([2*dim,2*dim])
        for i in range(dim):
            for j in range(dim):
                realm[2*i,2*j] = np.real(matrix[i,j])
                realm[2*i,2*j+1] = -np.imag(matrix[i,j])
                realm[2*i+1,2*j] = np.imag(matrix[i,j])
                realm[2*i+1,2*j+1] = np.real(matrix[i,j])
        return realm
    
#     def derivs(self, y, t, rhmatrix, rpumpmatrices, rabis, numatoms):
    def derivs(self, y, t, rhmatrix):
        ym = y.reshape(self.mshape)
#         local_rhmatrix = rhmatrix
#         for k in range(numatoms):
#             local_rhmatrix += rabis[k] * rpumpmatrices[k]
        derivs = ( rhmatrix @ ym - ym @ rhmatrix)
        return derivs.reshape(self.vdim)
    
    def rvec2cmatrixsoln(self):
        nrows, ncols = np.shape(self.realsoln)
        csoln = np.zeros([nrows, self.dim, self.dim], complex)
        realsoln = self.realsoln.reshape(nrows, 2 * self.dim, 2 * self.dim)
        for i in range(nrows):
            for j in range(self.dim):
                for k in range(self.dim):
                    csoln[i, j, k] = realsoln[i, 2*j, 2*k] + 1.j * realsoln[i, 2*j, 2*k+1]            
        return csoln

# The lindbladint class numerically integrates the Lindblad master equation.
# It takes a very similar for to the vonneumannint class, but also includes
# an arbitrary number of terms of Lindblad form to accomodate for dissipation
# in the system due to the interaction with an enviroment.
# Input variables are the initial density matrix, a list of times, and a
# Hamiltonian object. The initial density matrix should be an n by n
# numpy array, and the times a 1-d numpy array. The Hamiltonian object must
# have the following four attributes. 1. ham_obj.hmatrix must be an n by n
# complex numpy array representing the system Hamiltonian. 
# 2. ham_obj.lindbladgamma must be a dictionary of coefficients (i.e. damping
# rates) for the linblad terms. 3. ham_obj.lindbladraising must be a dictionary
# of raising operators for the linblad terms, represent by n by n matrices 
# using numpy arrays. Indexing for the ham_obj.lindbladraising dictionary 
# should correspond to the indexing of ham_obj.lindbladgamma. 
# 4. ham_obj.lindbladlowering must be a dictionary of lowering operators for
# the lindblad terms represented by n by n matrices using numpy arrays. Again,
# indexing for this dictionary should correspond with indexing for the
# ham_obj.lindbladgamma and ham_obj.lindbladraising dictionaries.
        
class lindbladint:
    
    def __init__(self, z0, t, ham_obj):
        self.z0 = z0
        self.t = t
        self.ham = ham_obj
        self.y0matrix = self.c2rmatrix(self.z0)
#        self.rihmatrix = self.c2rmatrix(-1.j * self.ham.ham(0))
#         self.rhmatrix = self.c2rmatrix(-1.j * self.ham.hmatrix)
#         self.rpumpmatrices = tuple([self.c2rmatrix(-1.j * self.ham.pumpmatrices[i]) for i in range(self.ham.numatoms)])
        
#         self.vdim = np.size(self.rhmatrix)
#         self.mshape = np.shape(self.rhmatrix)
        self.dim = self.ham.dim
        self.mshape = (2 * self.dim, 2 * self.dim)
        self.vdim = (2 * self.dim) ** 2
        
        self.y0 = self.y0matrix.reshape(self.vdim)
        
#         self.generate_rabis()
        
#     def rabi_t(self):
#         self.rabis = []
#         for k in range(self.ham.numatoms):
#             rabidictionary = {}
#             for t in self.t:
#                 rabidictionary[t] = self.ham.rabis[k]
#             self.rabis += [rabidictionary]

#    def generate_rabis(self):
#        self.rabis = []
#        for k in range(self.ham.numatoms):
#            rabidictionary = {0 : self.ham.rabis[k]}
#            self.rabis += [rabidictionary]
#
#    def lookup_rabi(self, t, rabi_dictionary):
#        change_times = list(rabi_dictionary.keys())
#    #     change_times.sort() # uncomment if times are not sorted
#        t_key = change_times[0]
#        for tau in change_times:
#            if t < tau:
#                break
#            t_key = tau
#        return rabi_dictionary[t_key]
        
    def solve(self):
#         self.rabis = []
#         for k in range(self.ham.numatoms):
#             rabidictionary = {}
#             for t in range(self.t):
#                 rabidictionary[t] = self.ham.rabis[k]
#             self.rabis += [rabidictionary] 

#         self.rabi_t()
        
        self.realsoln = odeint(self.derivs, self.y0, self.t)
#         self.realsoln = odeint(self.derivs, self.y0, self.t, 
#                                args = (self.rhmatrix, self.rpumpmatrices, 
#                                        self.rabis, self.ham.numatoms,))
#         self.realsoln = odeint(self.derivs, self.y0, self.t, args = (self.rhmatrix,))
        self.soln = self.rvec2cmatrixsoln()
        
    
    def c2rmatrix(self, matrix):
        dim = np.shape(matrix)[0]
#         Should be a square matrix. May want to put in some error handling
        realm = np.zeros([2*dim,2*dim])
        for i in range(dim):
            for j in range(dim):
                realm[2*i,2*j] = np.real(matrix[i,j])
                realm[2*i,2*j+1] = -np.imag(matrix[i,j])
                realm[2*i+1,2*j] = np.imag(matrix[i,j])
                realm[2*i+1,2*j+1] = np.real(matrix[i,j])
        return realm
    
    def com(self, matrix1, matrix2):
        return matrix1 @ matrix2 - matrix2 @ matrix1
    
    def acom(self, matrix1, matrix2):
        return matrix1 @ matrix2 + matrix2 @ matrix1
    
    def lindbladsops(self, ym, raisingm, loweringm, gamma):
        return gamma * ( loweringm @ ym @ raisingm - 0.5 * self.acom( raisingm @ loweringm, ym))
    
    def derivs(self, y, t):
        ym = y.reshape(self.mshape)
        sysderivs = self.com(self.c2rmatrix(-1.j * self.ham.ham(t)), ym)
#        For time independent Hamiltonians it is faster to use the following
#        sysderivs = self.com(self.rihmatrix, ym)
        
        lindbladderivs = np.zeros(self.mshape)
        for i in self.ham.lindbladgamma:
            gamma = self.ham.lindbladgamma[i]
            rrm = self.c2rmatrix(self.ham.lindbladraising[i])
            rlm = self.c2rmatrix(self.ham.lindbladlowering[i])
            lindbladderivs += self.lindbladsops(ym, rrm, rlm, gamma)
            
        derivs = sysderivs + lindbladderivs
#         derivs = sysderivs
        return derivs.reshape(self.vdim)
    
    def rvec2cmatrixsoln(self):
        nrows, ncols = np.shape(self.realsoln)
        csoln = np.zeros([nrows, self.dim, self.dim], complex)
        realsoln = self.realsoln.reshape(nrows, 2 * self.dim, 2 * self.dim)
        for i in range(nrows):
            for j in range(self.dim):
                for k in range(self.dim):
                    csoln[i, j, k] = realsoln[i, 2*j, 2*k] + 1.j * realsoln[i, 2*j, 2*k+1]            
        return csoln
    
       
class liouville:

    # z0 is the initial state of the system described by its density matrix.
    # z0 is represented by a 2d complex numpy array with the dimensions of the Hamiltonian
    # times is linearly spaced numpy array of times for which the system is solved for
    # This class is used to solve the dynamics for an open quantum system using the Liouvillian.
    # Note; to quickly access the steady state times can be a array of shape (2,) with a large final time 
    # as the evolution is 

    
    def __init__(self, z0, times, ham_obj):
        self.t = times
        self.z = z0 
        self.ham = ham_obj
        self.dim = self.ham.dim
        self.identity = np.identity(self.dim)
        self.flsize = self.dim**2
        self.vz = z0.T.reshape(z0.size)
        self.liouvillian = self.generate_liouvillian()
        
    def lindblad_sop(self):
        '''converts dissipative part of master equation into superoperator'''
        lindbladsop = np.zeros((self.ham.hmatrix.size, self.ham.hmatrix.size), complex)
        for i in self.ham.lindbladgamma:
            gamma = self.ham.lindbladgamma[i]
            lbr = self.ham.lindbladraising[i]
            lbl = self.ham.lindbladlowering[i]
            lindbladsop += gamma*(np.kron(lbl.conj(), lbl) - 0.5*np.kron(self.identity, lbr@lbl) - 0.5*np.kron((lbr@lbl).T, self.identity))
        return lindbladsop
    
    def system_sop(self):
        '''converts unitary part of master equation into superoperator'''
        return -1j*(np.kron(self.identity, self.ham.hmatrix) - np.kron(self.ham.hmatrix.T, self.identity))
        
    def generate_liouvillian(self):
        '''returns the liouvillian superoperator'''
        return self.system_sop() + self.lindblad_sop()
    
    def eigensolve(self):
        '''returns liouvillian eigenvalues and eigevectors'''
        self.evals, self.levecs, self.revecs = scipy.linalg.eig(self.liouvillian, right = True, left = True)
              
    def expL_op(self, t):
        L = np.zeros([self.dim**2, self.dim**2], complex)
        for i in range(self.dim**2):
            L[i,i] = np.exp(self.evals[i]*t)
        L = self.revecs @ L @ np.linalg.inv(self.revecs)
        return L
    
    def solve(self):
        '''calculates time evolution of z0 using liouvillian'''
        self.eigensolve()
        dt = self.t[1] - self.t[0]
        num_t = np.size(self.t)
        self.soln = np.zeros([num_t, self.dim, self.dim], complex)
        if np.linalg.det(self.revecs) == 0:
            print('The Liouvillian does not have linearly independent right eigenvectors: L is not diagonalisable')
        else: 
            for i, t in enumerate(self.t):
                self.soln[i] = self.vz.reshape((self.dim, self.dim)).T
                self.vz = self.expL_op(dt)@self.vz     
    
    
# test comment    
class unitaryevolution:
    # z0 is the initial state of the system described by its density matrix.
    # z0 is represented by a 2d complex numpy array with the dimensions of the Hamiltonian
    # times is linearly spaced numpy array of times for which the system is solved for
    # This class is used to solve the dynamics for a static Hamiltonian. It uses the 
    # the Hamiltonian function (of time) in the ham_obj, evaluated at time zero.
    def __init__(self, z0, times, ham_obj_or_matrix):
        
        self.t = times
        self.z = z0
        # self.ham = ham_obj
        # self.dim = self.ham.dim
        # self.hmatrix = self.ham.ham(0)
        self._get_ham(ham_obj_or_matrix)
        
    def _get_ham(self, ham_in):
        if type(ham_in) == np.ndarray:
            _get = self._get_hmatrix
        else:
            _get = self._get_ham_obj
        return _get(ham_in)
    
    def _get_hmatrix(self, ham_in):
        self.hmatrix = ham_in
        self.dim = len(ham_in)
        
    def _get_ham_obj(self, ham_in):
        self.dim = ham_in.dim
        self.hmatrix = ham_in.ham(0)
        
    def eigensolve(self):
        self.evals, self.evecs = np.linalg.eigh(self.hmatrix)
        
    def u_op(self, t):
        u = np.zeros([self.dim, self.dim], complex)
        for i in range(self.dim):
            u[i,i] = np.exp(-1.j * self.evals[i] * t)
        u = self.evecs @ u @ self.hc(self.evecs)
        return u
        
    def hc(self, np2darray):
        return np.transpose(np.conj(np2darray))
    
    def solve(self):
        self.eigensolve()
        dt = self.t[1] - self.t[0]
        num_t = np.size(self.t)
        u_dt = self.u_op(dt)
        ud_dt = self.hc(u_dt)
        self.soln = np.zeros([num_t, self.dim, self.dim], complex)
        for i, t in enumerate(self.t):
            self.soln[i] = self.z
            self.z = u_dt @ self.z @ ud_dt
            
    
class schrodint:
    # This class solves the Schrodinger equation given an initial wavevector
    # psi0, a list of times t, and a hamiltonian object.
    # The Hamiltonian object must
    # have the following four attributes. 1. ham_obj.ham.ham(0) must be an n by n
    # complex numpy array representing the system Hamiltonian. 
    # 2. ham_obj.lindbladgamma must be a dictionary of coefficients (i.e. damping
    # rates) for the linblad terms. 3. ham_obj.lindbladraising must be a dictionary
    # of raising operators for the linblad terms, represent by n by n matrices 
    # using numpy arrays. Indexing for the ham_obj.lindbladraising dictionary 
    # should correspond to the indexing of ham_obj.lindbladgamma. 
    # 4. ham_obj.lindbladlowering must be a dictionary of lowering operators for
    # the lindblad terms represented by n by n matrices using numpy arrays. Again,
    # indexing for this dictionary should correspond with indexing for the
    # ham_obj.lindbladgamma and ham_obj.lindbladraising dictionaries.
    def __init__(self, psi0, t, ham_obj):
        self.psi0 = psi0
        self.t = t
        self.ham = ham_obj
        self.dim = self.ham.dim
        self.y0 = self.c2r(self.psi0)
        self.rhmatrix = self.c2r(-1.j * self.ham.ham(0))
        self.ntimes = len(self.t)
    
    def solve(self):    
        self.realsoln = odeint(self.derivs, self.y0, self.t, args = (self.rhmatrix,))
        self.soln = self.reformatsolution(self.realsoln)
    
    def c2rmatrix(self, matrix):
#        dim = np.shape(matrix)[0]
        dim = self.dim
#         Should be a square matrix. May want to put in some error handling
        realm = np.zeros([2*dim,2*dim])
        for i in range(dim):
            for j in range(dim):
                realm[2*i,2*j] = np.real(matrix[i,j])
                realm[2*i,2*j+1] = -np.imag(matrix[i,j])
                # realm[2*i,2*j+1] = np.imag(matrix[i,j]) # brake to test unit_tests
                realm[2*i+1,2*j] = np.imag(matrix[i,j])
                realm[2*i+1,2*j+1] = np.real(matrix[i,j])
        return realm
    
    def c2rvector(self, vector):
#        dim = np.size(vector)
        dim = self.dim
#         Should be a square matrix. May want to put in some error handling
        realv = np.zeros([2*dim])
        for j in range(dim):
            realv[2*j] = np.real(vector[j])
            realv[2*j+1] = np.imag(vector[j])
        return realv
    
    def c2r(self, array):
        shape = np.shape(array)
        if len(shape) == 2:
            rarray = self.c2rmatrix(array)
        elif len(shape) == 1:
            rarray = self.c2rvector(array)
        return rarray
    
    def derivs(self, y, t, rhmatrix):
#    def derivs(self, y, t):
        derivs = rhmatrix @ y
#        derivs = self.c2r(-1.j * self.ham.ham(t)) @ y
        return derivs
        
    def reformatsolution(self, real_soln):
        n_times = len(real_soln)
        complex_soln = np.zeros([n_times, self.dim], complex)
        for i, c_vector in enumerate(complex_soln):
            for j in range(self.dim):
                c_vector[j] = real_soln[i,2*j] + 1.j * real_soln[i,2*j+1]
        return complex_soln
    
class quantumjumps(schrodint):
    # This class can be used to solve the quantm dynamics in an open quantum 
    # system using the Monte Carlo wavefunction method, a.k.a. quantum jumps.
    # It is initialised with an initial wavevector (psi0), a list of times t 
    # at which the solution is given, and a hamiltonian
    # object. This class is a child class of schrodint, using its methods to 
    # solve the time evolution for a given wavevector. 
    #
    # The Hamiltonian object must
    # have the following four attributes. 1. ham_obj.ham(0) must be an n by n
    # complex numpy array representing the system Hamiltonian. 
    # 2. ham_obj.lindbladgamma must be a dictionary of coefficients (i.e. damping
    # rates) for the linblad terms. 3. ham_obj.lindbladraising must be a dictionary
    # of raising operators for the linblad terms, represent by n by n matrices 
    # using numpy arrays. Indexing for the ham_obj.lindbladraising dictionary 
    # should correspond to the indexing of ham_obj.lindbladgamma. 
    # 4. ham_obj.lindbladlowering must be a dictionary of lowering operators for
    # the lindblad terms represented by n by n matrices using numpy arrays. Again,
    # indexing for this dictionary should correspond with indexing for the
    # ham_obj.lindbladgamma and ham_obj.lindbladraising dictionaries.
    def __init__(self, psi0, t, ham_obj, **kwargs):
        super().__init__(psi0, t, ham_obj)
        self.generate_nonhermitian_ham()
        self.t_max = self.t[-1]
        # Initiallise a dictionary of basis state wavevectors
        self.bstate_evolution = np.zeros([self.dim, self.ntimes, self.dim], complex)
        self.sampleratio = kwargs.pop('sampleratio', 1)
        
    def gen_bstate_evolution(self, **kwargs):
        # Before calculating the quantum trajectory from a given initial wavevector,
        # we calculate the non-unitary evolution of the basis states. Basis state
        # are written here like wavevectors with only one non-zero coefficient
        # which is set to 1.
        # map over bstate_i_evolution with multiprocessing capability
#        n_proc = kwargs.pop('n_proc', 1)
##        self.bstate_evolution = {}
#        with Pool(processes = n_proc) as p:
#            p.map(self.bstate_i_evolution, range(self.dim))
        # Old version with no multiprocessing
        # Initiallise a dictionary of basis state wavevectors
        self.bstate_evolution = {}
        for i in range(self.dim):
            # Make a wavevector psi0 which represents the i^th basis state.
            psi0 = np.zeros([self.dim], complex)
            psi0[i] = 1
            # convert from complex vector of dimension n to real one of dimension 2n.
            real_psi0 = self.c2r(psi0)
            # Calculate the non-unitary evolution of the state which starts in the 
            # i^th basis state.
            self.bstate_evolution[i] = self.deterministic_evolve(real_psi0, self.t[0])
            
    def bstate_i_evolution(self, i):
        # Before calculating the quantum trajectory from a given initial wavevector,
        # we calculate the non-unitary evolution of the basis states. Basis state
        # are written here like wavevectors with only one non-zero coefficient
        # which is set to 1.
        # Make a wavevector psi0 which represents the i^th basis state.
        psi0 = np.zeros([self.dim], complex)
        psi0[i] = 1
        # convert from complex vector of dimension n to real one of dimension 2n.
        real_psi0 = self.c2r(psi0)
        # Calculate the non-unitary evolution of the state which starts in the 
        # i^th basis state.
#        self.bstate_evolution[i] = self.deterministic_evolve(real_psi0, self.t[0])
        
        # Make new list of times with self.ntimes * self.sampleratio points
        # This is becuase odeint needs a smaller time step to find accurate solutions,
        # But such fine time steps mean that the results cannot be conveniently serialised
        # during multiprocessing (pool). However, the quantum jump trajectories do not require
        # such a small times step, so the output of odeint is down sampled to the original
        # list of times self.t
        times = np.linspace(0, self.t_max, (self.ntimes-1) * self.sampleratio + 1)
        realsoln = odeint(self.derivs, real_psi0, times, args = (self.rhmatrix,))
        soln = self.reformatsolution(realsoln)
        soln = soln[0::self.sampleratio]
#        print("Length of times is {}, and length of bstate_i_solution is {}" (np.size(times), np.size(soln)))
        
#        return self.deterministic_evolve(real_psi0, self.t[0])
        return soln
    
    def bstate_i_long_evolution(self, i):
        t1 = time.time()
        tstep = self.t[1] - self.t[0]
        times = np.linspace(0, tstep, self.sampleratio + 1)
        psi0 = np.zeros([self.dim], complex)
        psi0[i] = 1
        # convert from complex vector of dimension n to real one of dimension 2n.
        real_psi0 = self.c2r(psi0)
        downsampled_soln = np.zeros([self.ntimes, self.dim], complex)
        downsampled_soln[0] = psi0
        for j, t in enumerate(self.t[:-1]):
            realsoln = odeint(self.derivs, real_psi0, times, args = (self.rhmatrix,))
            soln = self.reformatsolution(realsoln)
            real_psi0 = realsoln[-1]
            downsampled_soln[j+1] = soln[-1]
        t2 = time.time()
        print("bstate {0} evolution calculated in {1} seconds".format(i, t2-t1))
        return downsampled_soln
            
            
            
            
    def quantum_trajectory(self, bstate_evolution):
        # This method builds a wavefunction trajectory with quantum jumps.
        # tpsi_amps are the list of coefficients to buiid the wavefunction from
        # using the time dependent state vectors of self.bstate_evolution.
        tpsi_amps = self.psi0
        # t_q_index is the index of the last quantum jump and is initially set to zero.
        t_q_index = 0
        # eta is a random variable [0,1] used to decide whether a jump has happened.
        eta = random()
        # initiallise array that describes the quantum trajectory.
        psi_array = np.zeros([self.ntimes, self.dim], complex)
        psi_array[0] = self.psi0
        
        # For loop runs over each time step, evolving the wavefunction forward,
        # and choosing whether a quantum jump is made.
        for km1, t_k in enumerate(self.t[1:]):
            # k is the current undex of psi_array for this timestep
            k = km1 + 1
            # tau_index is the number of time steps since last quantum jump
            tau_index = k - t_q_index
            # calculate the non-normalised wavevector tpsi_k for this (the k^th) step.
            tpsi_k = self.calc_tpsi_k(tpsi_amps, tau_index, bstate_evolution)
            # calculate the normalisation constants for tpsi_k
            tpsi_norm_sq = self.sqnorm(tpsi_k)
            tpsi_norm = np.sqrt(tpsi_norm_sq)
            # Choose whether a quantum jump has taken place
            if eta < tpsi_norm_sq:
                # No jump
                # normalise tpsi_k and record normalised wavefunction in psi_array
                psi_array[k] = tpsi_k / tpsi_norm
            else:
                # Quantum jump
                # Evaluate the quantum jump using self.quantum_jump method
                # This method returns the normalised wavefunction after the jump,
                # which is also equal to the new value of tpsi_amps.
                tpsi_amps = self.quantum_jump(tpsi_k)
                psi_array[k] = tpsi_amps
                # reset counter which marks the index of the last quantum jump
                t_q_index = k
        # returns quantum trajectory
        return psi_array
    
    def calc_rho(self, n):
        # estimate rho from n trajectories 
        bstate_evolution = self.bstate_evolution.copy()
        mean_rho = np.zeros([self.ntimes, self.dim, self.dim], complex)
        mean_rho_sq = np.zeros([self.ntimes, self.dim, self.dim], complex)
        for i in range(n):
            trajectory = self.quantum_trajectory(bstate_evolution)
            for k, psi_k in enumerate(trajectory):
                rho_k = np.outer(psi_k, np.conj(psi_k))
                mean_rho[k] += rho_k
                mean_rho_sq[k] += rho_k**2
        mean_rho = mean_rho/n
        mean_rho_sq = mean_rho_sq/n
        
        return mean_rho, mean_rho_sq
            
    def calc_tpsi_k(self, tpsi_amps, tau_index, bstate_evolution):
        # Calculates the non-unitary evolution of wavevector initially in state
        # psi0 = tpsi_amps, propagating forward by tau_index steps in time.
        # initiallise output state array
        tpsi_k = np.zeros([self.dim], complex)
        # Build tpsi_k from the evolution of the basis states after tau_index steps.
        for i, amp in enumerate(tpsi_amps):
            tpsi_k += amp * bstate_evolution[i][tau_index]
        return tpsi_k
    
    def quantum_jump(self, tpsi_in):
        # this method performs a quantum jump of tpsi_in, and returns the state
        # after the jump
        # Jump can be to any state produced by acting on tpsi_in with any of the 
        # lindbladlowering operators.
        # First calculate the probability of jumping to each of these states. These
        # probabilites sum to 1 as it is already decided a jump must occur when this
        # method is run.
        probs = self.jump_probs(tpsi_in)
        # Choose which lindbladlowering operator to apply based on a random choice,
        # weighted by the probs just defined above.
        lindblad_keys = list(self.ham.lindbladgamma.keys())
        lindblad_key = self.random_choice_weighted(lindblad_keys, probs)
        # perform quantum jump by acting with this lindbladlowering operator
        tpsi_out = self.ham.lindbladlowering[lindblad_key] @ tpsi_in
        # normalise tpsi_out.
        tpsi_out = tpsi_out / np.sqrt(self.sqnorm(tpsi_out))
        
        return tpsi_out
    
    def jump_probs(self, tpsi_k):
        # This method calculates the jump probabilites to each state a quantum
        # jump can be made to.
        probabilities = []
        for i in self.ham.lindbladgamma:
            # Evaluate the probability to make the jump using the i^th lindbladlowering operator.
            gamma_i = self.ham.lindbladgamma[i]
            a_i = self.ham.lindbladlowering[i]
            a_tpsi_k = a_i @ tpsi_k
            probabilities.append(gamma_i * self.sqnorm(a_tpsi_k))
        # normalise the probailites to 1.
        total = sum(probabilities)
        if total != 0:
            probabilities = [p / total for p in probabilities]
       
        return probabilities
    
    def random_choice_weighted(self, choices, weights):
        # returns random choice from choices, where each choice has probability given by weights
        # Assumes that sum of weights is 1.
        cumulative_prob = weights
        i = 1
        while i < len(cumulative_prob):
            cumulative_prob[i] += cumulative_prob[i-1]
            i += 1
        x = random()
        i = 0
        cp = cumulative_prob[i]
        while x > cp:
            i += 1
            cp = cumulative_prob[i]
        
        return choices[i]
            
    def deterministic_evolve(self, real_psi_i, t_i):
        # Evaluates the non-unitary evolution from initial state real_psi_i 
        # (expressed with real numbers with double the dimension to the complex wavector).
        # Assumes the initial state is at time t_i in the list of times self.t /
        i = np.where(np.isclose(self.t, t_i))[0][0]
        remainingtimes = self.t[i:]
        realsoln = odeint(self.derivs, real_psi_i, remainingtimes, args = (self.rhmatrix,))
        soln = self.reformatsolution(realsoln)
        return soln
        
        
    def generate_nonhermitian_ham(self):
        # Adds the non-Hermitian terms to the Hamiltonian matrix to solve the
        # non-unitary time evoltuion
        nonh_term = np.zeros([self.dim, self.dim], complex)
        
        for i in self.ham.lindbladgamma:
            gamma_i = self.ham.lindbladgamma[i]
            adag_i = self.ham.lindbladraising[i]
            a_i = self.ham.lindbladlowering[i]
            newterm = -0.5 * gamma_i * adag_i @ a_i
            nonh_term += newterm
            
        self.rhmatrix += self.c2r(nonh_term)
        
    def sqnorm(self, psi):
        # Returns the square norm of psi
        n2 = 0
        for c in psi:
            n2 += abs(c)**2
        return n2
            
'''
Solves a systme with a non-Hermitian (also works for Hermitian) Hamiltonian
by applying a similarity transformation to the state vector psi into a basis
where the Hamiltonian is diagonal.

input psi0 is the state vector at t=0.

Note this may not work in all cases since non-Hermitian matrices are not
always diagonalizeable. This may be the case when eigenvalues of the 
non-Hermitian matrix are degenerate - further work on this is needed and no
checks on the Hamiltonian matrix are made.
'''
class non_hermitian_unitaryevolution:
    def __init__(self, psi0, times, ham_obj_or_matrix):
        
        self.t = times
        self.psi0 = psi0
        # self.ham = ham_obj
        # self.dim = self.ham.dim
        # self.hmatrix = self.ham.ham(0)
        self._get_ham(ham_obj_or_matrix)
        
    def _get_ham(self, ham_in):
        if type(ham_in) == np.ndarray:
            _get = self._get_hmatrix
        else:
            _get = self._get_ham_obj
        return _get(ham_in)
    
    def _get_hmatrix(self, ham_in):
        self.hmatrix = ham_in
        self.dim = len(ham_in)
        
    def _get_ham_obj(self, ham_in):
        self.dim = ham_in.dim
        self.hmatrix = ham_in.ham(0)
        
    def eigensolve(self):
        self.evals, self.evecs = np.linalg.eig(self.hmatrix)
        # hhbar = self.hmatrix @ np.conj(self.hmatrix)
        # self.evals, self.evecs = np.linalg.eig(hhbar)
        # self.evals = np.sqrt(self.evals)
        self.s = np.transpose(self.evecs)
        self.si = np.linalg.inv(self.s)
        
    def u_op(self, t):
        u = np.zeros([self.dim, self.dim], complex)
        for i in range(self.dim):
            u[i,i] = np.exp(-1.j * self.evals[i] * t)
        u = self.si @ u @ self.s
        return u
    
    def u_op_matrixelement(self, i, j, t):
        diag_elements = np.exp(-1.j * self.evals * t)
        u = np.diag(diag_elements)
        u_ij = self.si[i] @ u @ self.s[:,j]
        return u_ij
        
    def hc(self, np2darray):
        return np.transpose(np.conj(np2darray))
    
    def solve(self):
        self.eigensolve()
        dt = self.t[1] - self.t[0]
        num_t = np.size(self.t)
        u_dt = self.u_op(dt)
        self.soln = np.zeros([num_t, self.dim], complex)
        psi = self.psi0
        for i, t in enumerate(self.t):
            self.soln[i] = psi
            psi = u_dt @ psi
        
            
        
            
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    