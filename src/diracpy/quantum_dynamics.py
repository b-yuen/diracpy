#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""Solves quantum dyanmics.

Solves quantum dyanmics of :mod:`quantum systems <diracpy.quantum_systems>`
objects. Multiple dyanmical models are given by the different
classes in this module. Each class shares a commont interface comprising of
initial state, a list of times to solve the dynamics for, and the 
:class:`diracpy.quantum_systems.qsys` object that describes the system.
The inital state is either the inital state vector or the density matrix, 
depending on the dynamical model. In both cases, these should be vectorised
in the basis constructed by the given :class:`diracpyt.quantum_systems.qsys` 
object. For systems with a time dependent hamiltonian, a 
:class:`diracpy.quantum_systems.qsys_t` object is given instead.
In this case only dynamical models that use scipy.integrate.odeint rather
matrix diagonalisatin should be used -- that are the :class:`vonneumannint`, 
:class:`lindbladint`,  :class:`schrodint` or :class:`quantumjumps` models.

Classes
-------
    vonneumannint
    lindbladint
    schrodint
    quantumjumps
    unitaryevolution
    non_hermitian_unitaryevolution
    liouville
"""


# Things to do:
# Make lindblad a inherit from vonneuman int.

# This script defines classes which can be used to
# numerically solve quantum dynamics

import numpy as np
import scipy
from scipy.integrate import odeint
#from random import random
import random
from multiprocessing import Pool
import time
import diracpy.quantum_systems

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
    """
    Solve von Neumann equation.
    
    Solves the dynamics for a given :class:`diracpy.quantum_systems.qsys` object
    by numerically integrating the von Neumann equation.
    
    Attributes
    ----------
    t : numpy.ndarray
        1d array of times to solve the system for.
    qsys : :class:`diracpy.quantum_systems.qsys` or :class:`diracpy.quantum_systems.qsys_t`
        The qsys object that describes the quantum system to be solved
    rho0 : numpy.ndarray
        2d array representing the initial density matrix of the system in the
        basis given by qsys. If the initial density operator is given as
        a :class:`diracpy.states_operators.qop` object then the density
        matrix is generated using the `qsys.vectorize` method.
    soln : numpy.ndarray
        3d array containing the density matrices of the solved system for
        each time of the attribute `t`.
    """
    
    def __init__(self, rho0, t, qsys, solve=True):
        """
        Solve von Neumann equation.
        
        Initialises the object and solves the system for the intial density
        operator rho0, at each time in t. The system is solve by expanding
        the von Neumann equation into an equivalent set of coupled ODEs which
        are then solved using `scipy.integrate.odeint <https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.odeint.html>`_.

        Parameters
        ----------
        rho0 : :class:`diracpy.states_operators.qop` or numpy.ndarray
            Initial density operator or density matrix of the system.
        t : numpy.ndarray
            1D array of times (`int` or `float`) to solve the system for.
        qsys : :class:`diracpy.quantum_systems.qsys`
            Quantum system to solve, which contains all information on the
            Hamiltonian and the basis used for the calculation.
        solve : bool, optional
            System solved on initialisation if True, and not solved otherwise.
            The default is True.

        Returns
        -------
        None.
        """
        # setup simulation
        self.t = t
        self.qsys = qsys
        self._generate_rhmatrix()
        # define dimensions for equivalent real variable system
        self._dim = self.qsys.dim
        self._mshape = (2 * self._dim, 2 * self._dim)
        self._vdim = (2 * self._dim) ** 2
        # define initial state
        self.set_initial_state(rho0)
        self._y0matrix = self._c2rmatrix(self.rho0)
        self._y0 = self._y0matrix.reshape(self._vdim)
        # solve
        if solve:
            self.solve()
            
    def set_initial_state(self, rho0):
        """
        Define initial state.
        
        Define the initial condition of the von Neumann equation given by
        the inital density matrix of density operator.

        Parameters
        ----------
        rho0 : :class:`diracpy.states_operators.qop` or numpy.ndarray
            Initial density operator or density matrix of the system.

        Returns
        -------
        None.
        """
        try:
            rho = self.qsys.vectorize(rho0)
        except AttributeError:
            # handles depricated qsys functionality where qsys was an object
            # describing the Hamiltonian but had no method vectorize.
            rho = rho0
        self.rho0 = rho
        
    def solve(self):
        """
        Solve system.
        
        Solves the system for array of times, `t`, using
        `scipy.integrate.odeint <https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.odeint.html>`_.
        The solution is an array of density matrices. These are stored in
        the `soln` attribute which is a 3d numpy.ndarray.
        
        Returns
        -------
        None.
        """
        realsoln = odeint(self._derivs, self._y0, self.t)
        self.soln = self._rvec2cmatrixsoln(realsoln)
        
    def _generate_rhmatrix(self):
        # self._rhmatrix = self._c2rmatrix(-1.j * self.qsys.hmatrix)
        self._rhmatrix = self._c2rmatrix(-1.j * self.qsys.ham(0))
        
        if 'pumpmatrices' in self.qsys.__dict__:
            self._rhmatrix = self._append_pump(self._rhmatrix)
                
    def _append_pump(self, rhmatrix):
        # depricated - used in predecessor or diracpy - kept for backwards
        # compatibility
        rhmatrix_out = rhmatrix
        rpumpmatrices = tuple([self._c2rmatrix(-1.j * self.qsys.pumpmatrices[i]) for i in range(self.qsys.numatoms)])
        for k in range(self.qsys.numatoms):
            rhmatrix_out += self.qsys.rabis[k] * rpumpmatrices[k]
        return rhmatrix_out
    
    def _c2rmatrix(self, matrix):
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
    
    def _com(self, matrix1, matrix2):
        return matrix1 @ matrix2 - matrix2 @ matrix1
    
    # def _derivs(self, y, t):
    #     ym = y.reshape(self._mshape)
    #     derivs = ( self._rhmatrix @ ym - ym @ self._rhmatrix)
    #     return derivs.reshape(self._vdim)
    
    def _derivs(self, y, t):
        ym = y.reshape(self._mshape)
        if type(self.qsys) == diracpy.quantum_systems.qsys:
            derivs = self._com(self._rhmatrix, ym)
        else:
            derivs = self._com(self._c2rmatrix(-1.j * self.qsys.ham(t)), ym)
#        For time independent Hamiltonians it is faster to use the following
#        derivs = self._com(self.rihmatrix, ym)
        return derivs.reshape(self._vdim)
    
    def _rvec2cmatrixsoln(self, realsoln):
        nrows, ncols = np.shape(realsoln)
        csoln = np.zeros([nrows, self._dim, self._dim], complex)
        realsoln = realsoln.reshape(nrows, 2 * self._dim, 2 * self._dim)
        for i in range(nrows):
            for j in range(self._dim):
                for k in range(self._dim):
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
        
class lindbladint(vonneumannint):
    """
    Solve Lindblad master equation.
    
    Solves the dynamics for a given :class:`diracpy.quantum_systems.qsys` or
    the :class:`diracpy.quantum_systems.qsys_t` object
    by numerically integrating the Lindblad master equation using 
    `scipy.integrate.odeint <https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.odeint.html>`_.
    This class
    inherits from :class:`vonneumannint` due to the similarity between these
    equations of motion. The differences are: (1) that :class:`lindbladint`
    accomodates open quantum systems by including the optional jump_ops
    when initialising the quantum system 
    (see documentation :mod:`diracpy.quantum_systems`), and (2) that 
    :class:`time dependent systems <diracpy.quantum_systems.qsys_t>`
    can be solved with class. Note that the system is not solved on
    initililsation here --  the solve method must subsequently be used.
    
    Attributes
    ----------
    t : numpy.ndarray
        1d array of times to solve the system for.
    qsys : :class:`diracpy.quantum_systems.qsys` or :class:`diracpy.quantum_systems.qsys_t`
        The qsys object that describes the quantum system to be solved
    rho : numpy.ndarray
        2d array representing the initial density matrix of the system in the
        basis given by qsys. If the initial density operator is given as
        a :class:`diracpy.states_operators.qop` object then the density
        matrix is generated using the `qsys.vectorize` method.
    soln : numpy.ndarray
        3d array containing the density matrices of the solved system for
        each time of the attribute `t`.
    
    """
    
    def __init__(self, rho0, t, qsys):
        """
        Initialise the Lindblad equation solver.
        
        Initialises the object with the intial density
        operator rho0, a list of times to solve for, and the
        quantum system to be solved. This sets up the set of coupled
        ordinary differential equations equivalent to the Lindblad
        master equation.
        
        Parameters
        ----------
        rho0 : :class:`diracpy.states_operators.qop` or numpy.ndarray
            Initial density operator or density matrix of the system.
        t : numpy.ndarray
            1D array of times (`int` or `float`) to solve the system for.
        qsys : :class:`diracpy.quantum_systems.qsys` or :class:`diracpy.quantum_systems.qsys_t`
            Quantum system to solve, which contains all information on the
            Hamiltonian and the basis used for the calculation.
        solve : bool, optional
            System solved on initialisation if True, and not solved otherwise.
            The default is True.

        Returns
        -------
        None.
        """
        super().__init__(rho0, t, qsys, solve=False)
    
    def _com(self, matrix1, matrix2):
        return matrix1 @ matrix2 - matrix2 @ matrix1
    
    def _acom(self, matrix1, matrix2):
        return matrix1 @ matrix2 + matrix2 @ matrix1
    
    def _lindbladsops(self, ym, raisingm, loweringm, gamma):
        return gamma * ( loweringm @ ym @ raisingm - 0.5 * 
                        self._acom( raisingm @ loweringm, ym))
    
    def _derivs(self, y, t):
        ym = y.reshape(self._mshape)
        if type(self.qsys) == diracpy.quantum_systems.qsys:
            sysderivs = self._com(self._rhmatrix, ym)
        else:
            sysderivs = self._com(self._c2rmatrix(-1.j * self.qsys.ham(t)), ym)
        # sysderivs = self._com(self._c2rmatrix(-1.j * self.qsys.ham(t)), ym)
#        For time independent Hamiltonians it is faster to use the following
#        sysderivs = self._com(self.rihmatrix, ym)
        
        lindbladderivs = np.zeros(self._mshape)
        for i in self.qsys.lindbladgamma:
            gamma = self.qsys.lindbladgamma[i]
            rrm = self._c2rmatrix(self.qsys.lindbladraising[i])
            rlm = self._c2rmatrix(self.qsys.lindbladlowering[i])
            lindbladderivs += self._lindbladsops(ym, rrm, rlm, gamma)
            
        derivs = sysderivs + lindbladderivs
        return derivs.reshape(self._vdim)

            
    
class schrodint:
    """
    Solve the Schrodinger equation.
    
    This class solves the Schrodinger equation given an initial wavevector
    psi0, a list of times t, and a :class:`diracpy.quantum_systems.qsys` 
    or :class:`diracpy.quantum_systems.qsys_t` object that contains the
    hamiltonian and system basis. The Schrodinger equatoin is solved by
    numerical intergrate using 
    `scipy.integrate.odeint <https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.odeint.html>`_.
    
    Attributes
    ----------
    t : numpy.ndarray
        1d array of times to solve the system for.
    qsys : :class:`diracpy.quantum_systems.qsys` or :class:`diracpy.quantum_systems.qsys_t`
        The qsys object that describes the quantum system to be solved
    psi0 : numpy.ndarray
        1d array representing the initial state of the system in the
        basis given by qsys. If the initial density operator is given as
        a :class:`diracpy.states_operators.ket` object then this vector
        is generated using the `qsys.vectorize` method.
    soln : numpy.ndarray
        2d array containing the state vectors of the solved system for
        each time of the attribute `t`.
        
    Examples
    --------
    Solve the Jaynes-Cummings model for an initially excited atom.
    
    First we need to setup the :class:'qsys <quantum_systems.qsys>'.
    
    >>> # define state space
    >>> atom = dp.two_level_subspace(index=0)
    >>> cav = dp.fock_subspace(index=1)
    >>> # define model parameters
    >>> Delta, g = 0, np.pi
    >>> # define Hamiltonian
    >>> H_0 = Delta * atom.sigma_z
    >>> V = g * (cav.a * atom.sigma_plus + cav.adag * atom.sigma_minus)
    >>> H = H_0 + V
    >>> # define intial state and qsys.
    >>> psi0 = dp.ket(['e',0])
    >>> system = dp.qsys(H, initialstates=[psi0], n_int=2)
    
    Construct :class:`schrodint` solver object for this system:
        
    >>> schsolver = dp.schrodint(psi0, times, system)
    
    Solve Schrodinger equation:
        
    >>> schsolver.solve()
    
    Return solution:
        
    >>> schsolver.soln
    
    """
    
    def __init__(self, psi0, t, qsys):
        """
        Initialise the Schrodinger equation solver.
        
        Initialises the object with the intial state vector
        psi0, a list of times to solve for, and the
        quantum system to be solved. This sets up the set of coupled
        ordinary differential equations equivalent to the Schrodinger
        equation.
        
        Parameters
        ----------
        psi0 : :class:`diracpy.states_operators.ket` or numpy.ndarray
            Initial state ket or vector of the system.
        t : numpy.ndarray
            1D array of times (`int` or `float`) to solve the system for.
        qsys : :class:`diracpy.quantum_systems.qsys` or :class:`diracpy.quantum_systems.qsys_t`
            Quantum system to solve, which contains all information on the
            Hamiltonian and the basis used for the calculation.
        dim : int
            The dimension of the system state space.

        Returns
        -------
        None.
        """
        
        # self.psi0 = psi0
        self.t = t
        self.qsys = qsys
        # system dimension
        self.dim = self.qsys.dim
        # Static Hamiltonian matrix in equivalent real system of equations
        self._rhmatrix = self._c2r(-1.j * self.qsys.ham(0))
        # number of timesteps
        self._ntimes = len(self.t)
        # define initial state
        self.set_initial_state(psi0)
        self._y0 = self._c2r(self.psi0)
    
    def set_initial_state(self, psi0):
        """
        Define initial state.
        
        Define the initial condition of the Schrodinger equation given by
        the inital state vector, `psi0`.

        Parameters
        ----------
        psi0 : :class:`diracpy.states_operators.ket` or numpy.ndarray
            Initial state vector of the system.

        Returns
        -------
        None.
        """
        self.psi0 = self.qsys.vectorize(psi0)
    
    def solve(self):
        """
        Solve Schrodinger equation.
        
        Solve the Schrodinger equation using `scipy.integrate.odeint`.
        Solution returned to `soln` attribute.

        Returns
        -------
        None.

        """
        # self.realsoln = odeint(self._derivs, self._y0, self.t, args = (self._rhmatrix,))
        self.realsoln = odeint(self._derivs, self._y0, self.t)
        self.soln = self._reformatsolution(self.realsoln)
    
    def _c2rmatrix(self, matrix):
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
    
    def _c2rvector(self, vector):
#        dim = np.size(vector)
        dim = self.dim
#         Should be a square matrix. May want to put in some error handling
        realv = np.zeros([2*dim])
        for j in range(dim):
            realv[2*j] = np.real(vector[j])
            realv[2*j+1] = np.imag(vector[j])
        return realv
    
    def _c2r(self, array):
        shape = np.shape(array)
        if len(shape) == 2:
            rarray = self._c2rmatrix(array)
        elif len(shape) == 1:
            rarray = self._c2rvector(array)
        return rarray
    
    def _derivs(self, y, t):
        derivs = self._rhmatrix @ y
        if type(self.qsys) == diracpy.quantum_systems.qsys:
            derivs = self._rhmatrix @ y
        else:
            rhm_t = self._c2r(-1.j * self.qsys.ham(t))
            derivs = rhm_t @ y
        return derivs
        
    def _reformatsolution(self, real_soln):
        # n_times = len(real_soln)
        # complex_soln = np.zeros([n_times, self.dim], complex)
        ntimes = len(real_soln)
        complex_soln = np.zeros([ntimes, self.dim], complex)
        for i, c_vector in enumerate(complex_soln):
            for j in range(self.dim):
                c_vector[j] = real_soln[i,2*j] + 1.j * real_soln[i,2*j+1]
        return complex_soln
    
    
    
class quantumjumps(schrodint):
    """
    Qunatum jump simulation.
    
    This class solves the quantum dynamics stochastically using the quantum
    jump (a.k.a. Mote Carlo Wavefunction) method. This simulation builds
    quantum trajectories with periods of deterministic evolution
    interspersed with quantum jumps. The density matrix of the system is
    estimated by the mean of a specified number of trajectories.
    
    The simulation is built for a initial state psi0, a list of times t, and
    a :class:`diracpy.quantum_systems.qsys` that contains the Hamiltonian
    and system basis.
    
    Once the object is constructed the deterministic evolution is next
    calculated using :func:`gen_bstate_evolution`, before finally running
    the stochastic simulation with :func:`calc_rho`.
    
    Attributes
    ----------
    t : numpy.ndarray
        1d array of times to solve the system for.
    qsys : :class:`diracpy.quantum_systems.qsys` or :class:`diracpy.quantum_systems.qsys_t`
        The qsys object that describes the quantum system to be solved
    psi0 : numpy.ndarray
        1d array representing the initial state of the system in the
        basis given by qsys. If the initial density operator is given as
        a :class:`diracpy.states_operators.ket` object then this vector
        is generated using the `qsys.vectorize` method.
    soln : numpy.ndarray
        2d array containing the state vectors of the solved system for
        each time of the attribute `t`.
    bstate_evolution : dict of 2d numpy.ndarray
        Dictionary containing the time evolution of each basisstate.
        Automatically populated using :func:`gen_bstate_evolution`, but
        can be manually populated using :func:`bstate_i_evolution` or
        :func:`bstate_i_long_evolution`.
        
    Examples
    --------
    Solve the open Jaynes-Cummings model for an initially excited atom.
    
    First we need to setup the :class:'qsys <quantum_systems.qsys>'.
    
    >>> # define state space
    >>> atom = dp.two_level_subspace(index=0)
    >>> cav = dp.fock_subspace(index=1)
    >>> # define model parameters
    >>> Delta, g, gamma = 0, np.pi, 0.5
    >>> # define Hamiltonian
    >>> H_0 = Delta * atom.sigma_z
    >>> V = g * (cav.a * atom.sigma_plus + cav.adag * atom.sigma_minus)
    >>> H = H_0 + V
    >>> # define the jump operator
    >>> jump = np.sqrt(gamma)*cav.a
    >>> # define intial state and qsys.
    >>> psi0 = dp.ket(['e',0])
    >>> system = dp.qsys(H, initialstates=[psi0], n_int=2, jump_ops=jump)
    
    Construct :class:`quantumjumps` solver object for this system:
        
    >>> jpsolver = dp.quantumjumps(psi0, times, system)
    
    Solve the deterministic evolution
        
    >>> jpsolver.gen_bstate_evolution()
    
    Run stochastic simulation with random quantum jumps from 1000 quantum
    trajectories.
        
    >>> mean_rho = jpsolver.calc_rho(1000)
    
    The mean value of rho (an unbaised estimator for rho) at each time is
    returned to `mean_rho`.
    
    """
    
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
    def __init__(self, psi0, t, qsys, test=False, **kwargs):
        """
        Initialise the quantum jump solver.
        
        Initialises the object with the intial state vector
        psi0, a list of times to solve for, and the
        quantum system to be solved. This sets up the set of coupled
        ordinary differential equations for the deterministic evolution
        and creates a random generator for the jump processes.
        
        Parameters
        ----------
        psi0 : :class:`diracpy.states_operators.ket` or numpy.ndarray
            Initial state ket or vector of the system.
        t : numpy.ndarray
            1D array of times (`int` or `float`) to solve the system for.
        qsys : :class:`diracpy.quantum_systems.qsys` or :class:`diracpy.quantum_systems.qsys_t`
            Quantum system to solve, which contains all information on the
            Hamiltonian and the basis used for the calculation.
        dim : int
            The dimension of the system state space.
        test : bool, optional.
            When true then the random generator for calculating jump
            processes is seeded with 0 such that the stochastic simulation
            gives repeatable results, enabling tests for this class.
            The default it false.

        Returns
        -------
        None.
        """
        super().__init__(psi0, t, qsys)
        self._generate_nonhermitian_ham()
        self.t_max = self.t[-1]
        # Initiallise a dictionary of basis state wavevectors
        self.bstate_evolution = np.zeros([self.dim, self._ntimes, self.dim], complex)
        self.sampleratio = kwargs.pop('sampleratio', 1)
        self._random = self._random_generator(test)
        # self._test_mode(test)
        
    def gen_bstate_evolution(self):
        """
        Calculate deterministic evolution for basis states.
        
        Calculates the determinstic evolution for each basis state, which is
        how the basis states evolve between jumps. Since all states are 
        a superposition of basis states, the determinstic evolution of any
        state can be constructed from the evolution of basis states.
        
        The evolution of each basis states is determined by a non-Hermitian
        version of the system Hamiltonian, and the evolution found using
        scipy.integrate.odeint.
        
        The result is returned to the attribute `bstate_evolution` which
        is an 3d numpy.ndarray of state vectors (axis 2) at each time 
        (axis 1) for each basis state (axis 0).

        Returns
        -------
        None.

        """
        # Before calculating the quantum trajectory from a given initial wavevector,
        # we calculate the non-unitary evolution of the basis states. 
        # Basis state are written here like wavevectors with only one 
        # non-zero coefficient which is set to 1.
        for i in range(self.dim):
            # Make a wavevector psi0 which represents the i^th basis state.
            psi0 = np.zeros([self.dim], complex)
            psi0[i] = 1
            # convert from complex vector of dimension n to real one of dimension 2n.
            real_psi0 = self._c2r(psi0)
            # Calculate the non-unitary evolution of the state which starts in the 
            # i^th basis state.
            self.bstate_evolution[i] = self._deterministic_evolve(real_psi0, self.t[0])
            
    def bstate_i_evolution(self, i):
        r"""
        Calculate evolution of :math:`i^\mathrm{th}` basis state.
        
        Lower level variation of :func:`gen_bstate_evolution` which allows
        one to calculate the evolution specifically for the 
        :math:`i^\mathrm{th}`. This is useful when used in conjunction with 
        multiprocessing using the
        `pathos.ParallelPool <https://pathos.readthedocs.io/en/latest/pathos.html#pathos.parallel.ParallelPool>`_.

        Parameters
        ----------
        i : int
            The index of the basis state for which to find the time evolution.

        Returns
        -------
        soln : numpy.ndarray
            2d array containing the basis state's time evolution. Each
            entry is the state vector for the evolving state at one instance
            in time.

        """
        # Before calculating the quantum trajectory from a given initial wavevector,
        # we calculate the non-unitary evolution of the basis states. Basis state
        # are written here like wavevectors with only one non-zero coefficient
        # which is set to 1.
        # Make a wavevector psi0 which represents the i^th basis state.
        psi0 = np.zeros([self.dim], complex)
        psi0[i] = 1
        # convert from complex vector of dimension n to real one of dimension 2n.
        real_psi0 = self._c2r(psi0)
        # Calculate the non-unitary evolution of the state which starts in the 
        # i^th basis state.
#        self.bstate_evolution[i] = self._deterministic_evolve(real_psi0, self.t[0])
        
        # Make new list of times with self._ntimes * self.sampleratio points
        # This is becuase odeint needs a smaller time step to find accurate solutions,
        # But such fine time steps mean that the results cannot be conveniently serialised
        # during multiprocessing (pool). However, the quantum jump trajectories do not require
        # such a small times step, so the output of odeint is down sampled to the original
        # list of times self.t
        times = np.linspace(0, self.t_max, (self._ntimes-1) * self.sampleratio + 1)
        # realsoln = odeint(self._derivs, real_psi0, times, args = (self._rhmatrix,))
        realsoln = odeint(self._derivs, real_psi0, times)
        soln = self._reformatsolution(realsoln)
        soln = soln[0::self.sampleratio]
#        print("Length of times is {}, and length of bstate_i_solution is {}" (np.size(times), np.size(soln)))
        
#        return self._deterministic_evolve(real_psi0, self.t[0])
        return soln
    
    def bstate_i_long_evolution(self, i, sampleratio, verbose=False):
        r"""
        Calculate evolution of $i^\mathrm{th}$ basis state.
        
        Calculates basis state evolution for state `i`, but with
        time steps that are `sampleratio` times smaller than those in
        the attribute `t`. This is because the stochastic simulation
        can sometimes require a much lower time resolution than `odeint`
        requires to accurately capture the deterministic time evolution.
        Such a situation occurs when far off-resonant states ar present,
        and duration of the the time evolution required is long by 
        comparison.

        Parameters
        ----------
        i : int
            The index of the basis state for which to find the time evolution.
        verbose : bool, optional
            Prints updates of calculation to screen so progress over 
            a large number of states can be monitored. The default is False.

        Returns
        -------
        soln : numpy.ndarray
            2d array containing the basis state's time evolution. Each
            entry is the state vector for the evolving state at one instance
            in time.

        """
        t1 = time.time()
        tstep = self.t[1] - self.t[0]
        times = np.linspace(0, tstep, sampleratio + 1)
        psi0 = np.zeros([self.dim], complex)
        psi0[i] = 1
        # convert from complex vector of dimension n to real one of dimension 2n.
        real_psi0 = self._c2r(psi0)
        downsampled_soln = np.zeros([self._ntimes, self.dim], complex)
        downsampled_soln[0] = psi0
        for j, t in enumerate(self.t[:-1]):
            # realsoln = odeint(self._derivs, real_psi0, times, args = (self._rhmatrix,))
            realsoln = odeint(self._derivs, real_psi0, times)
            soln = self._reformatsolution(realsoln)
            real_psi0 = realsoln[-1]
            downsampled_soln[j+1] = soln[-1]
        t2 = time.time()
        if verbose:
            print("bstate {0} evolution calculated in {1} seconds".format(i, t2-t1))
        return downsampled_soln
            
    # def _test_mode(self, test):
    #     self.test = test
    #     if test:
    #         self._seed_generator = _counter()
    #     else:
    #         pass
            
    # def _generate_eta(self):
    #     if self.test:
    #         eta = random([self._seed_generator()])
    #     else:
    #         eta = random()
    #     return eta

    def _random_generator(self, test):
        random_generator = random.Random()
        if test:
            random_generator.seed(0)
        else:
            pass
        return random_generator
            
            
    # def quantum_trajectory(self, bstate_evolution):
    def quantum_trajectory(self):
        """
        Calculate quantum trajectory.
        
        Calculates a single quantum trajectory from the intial state of
        the system. Each trajectory is given be periods of deterministic 
        evolution interspersed randomly with quantum jumps. The jump 
        probability distribution is calculated from the current state of 
        the simulated wave function trajectory at each time step.

        Returns
        -------
        psi_array : numpy.ndarray
            Quantum trajectory.

        """
        bstate_evolution = self.bstate_evolution
        # This method builds a wavefunction trajectory with quantum jumps.
        # tpsi_amps are the list of coefficients to buiid the wavefunction from
        # using the time dependent state vectors of self.bstate_evolution.
        tpsi_amps = self.psi0
        # t_q_index is the index of the last quantum jump and is initially set to zero.
        t_q_index = 0
        # eta is a random variable [0,1] used to decide whether a jump has happened.
        # eta = random()
        eta = self._random.random()
        # initiallise array that describes the quantum trajectory.
        psi_array = np.zeros([self._ntimes, self.dim], complex)
        psi_array[0] = self.psi0
        
        # For loop runs over each time step, evolving the wavefunction forward,
        # and choosing whether a quantum jump is made.
        for km1, t_k in enumerate(self.t[1:]):
            # k is the current undex of psi_array for this timestep
            k = km1 + 1
            # tau_index is the number of time steps since last quantum jump
            tau_index = k - t_q_index
            # calculate the non-normalised wavevector tpsi_k for this (the k^th) step.
            tpsi_k = self._calc_tpsi_k(tpsi_amps, tau_index, bstate_evolution)
            # calculate the normalisation constants for tpsi_k
            tpsi_norm_sq = self._sqnorm(tpsi_k)
            tpsi_norm = np.sqrt(tpsi_norm_sq)
            # Choose whether a quantum jump has taken place
            if eta < tpsi_norm_sq:
                # No jump
                # normalise tpsi_k and record normalised wavefunction in psi_array
                psi_array[k] = tpsi_k / tpsi_norm
            else:
                # Quantum jump
                # Evaluate the quantum jump using self._quantum_jump method
                # This method returns the normalised wavefunction after the jump,
                # which is also equal to the new value of tpsi_amps.
                tpsi_amps = self._quantum_jump(tpsi_k)
                psi_array[k] = tpsi_amps
                # reset counter which marks the index of the last quantum jump
                t_q_index = k
                # Update eta for next quantum jump
                # eta = random()
                eta = self._random.random()
        # returns quantum trajectory
        return psi_array
    
    def calc_rho(self, n):
        """
        Run quantum jump simulation.
        
        Run simulation, taking the average of n quantum trajectories.
        Each trajectory is given be periods of deterministic evolution
        interspersed randomly with quantum jumps. The jump probability
        distribution is calculated from the current state of the simulated
        wave function trajectory at each time step.

        Parameters
        ----------
        n : int
            Number of quantum trajectories to use in simulation.

        Returns
        -------
        mean_rho_sq : numpy.ndarray
            Mean value of density matrix for each time.

        """
        # estimate rho from n trajectories 
        # estimate of variance of rho not yet implemented
        bstate_evolution = self.bstate_evolution.copy()
        mean_rho = np.zeros([self._ntimes, self.dim, self.dim], complex)
        #mean_rho_sq = np.zeros([self._ntimes, self.dim, self.dim], float)
        for i in range(n):
            # trajectory = self.quantum_trajectory(bstate_evolution)
            trajectory = self.quantum_trajectory()
            for k, psi_k in enumerate(trajectory):
                rho_k = np.outer(psi_k, np.conj(psi_k))
                mean_rho[k] += rho_k
                #mean_rho_sq[k] += np.real(np.multiply(rho_k, rho_k.conj()))
        mean_rho = mean_rho/n
        #mean_rho_sq = mean_rho_sq/n
        
        return mean_rho#, mean_rho_sq
            
    def _calc_tpsi_k(self, tpsi_amps, tau_index, bstate_evolution):
        # Calculates the non-unitary evolution of wavevector initially in state
        # psi0 = tpsi_amps, propagating forward by tau_index steps in time.
        # initiallise output state array
        tpsi_k = np.zeros([self.dim], complex)
        # Build tpsi_k from the evolution of the basis states after tau_index steps.
        for i, amp in enumerate(tpsi_amps):
            tpsi_k += amp * bstate_evolution[i][tau_index]
        return tpsi_k
    
    def _quantum_jump(self, tpsi_in):
        # this method performs a quantum jump of tpsi_in, and returns the state
        # after the jump
        # Jump can be to any state produced by acting on tpsi_in with any of the 
        # lindbladlowering operators.
        # First calculate the probability of jumping to each of these states. These
        # probabilites sum to 1 as it is already decided a jump must occur when this
        # method is run.
        probs = self._jump_probs(tpsi_in)
        # Choose which lindbladlowering operator to apply based on a random choice,
        # weighted by the probs just defined above.
        lindblad_keys = list(self.qsys.lindbladgamma.keys())
        lindblad_key = self._random_choice_weighted(lindblad_keys, probs)
        # perform quantum jump by acting with this lindbladlowering operator
        tpsi_out = self.qsys.lindbladlowering[lindblad_key] @ tpsi_in
        # normalise tpsi_out.
        tpsi_out = tpsi_out / np.sqrt(self._sqnorm(tpsi_out))
        
        return tpsi_out
    
    def _jump_probs(self, tpsi_k):
        # This method calculates the jump probabilites to each state a quantum
        # jump can be made to.
        probabilities = []
        for i in self.qsys.lindbladgamma:
            # Evaluate the probability to make the jump using the i^th lindbladlowering operator.
            gamma_i = self.qsys.lindbladgamma[i]
            a_i = self.qsys.lindbladlowering[i]
            a_tpsi_k = a_i @ tpsi_k
            probabilities.append(gamma_i * self._sqnorm(a_tpsi_k))
        # normalise the probailites to 1.
        total = sum(probabilities)
        if total != 0:
            probabilities = [p / total for p in probabilities]
       
        return probabilities
    
    def _random_choice_weighted(self, choices, weights):
        # returns random choice from choices, where each choice has probability given by weights
        # Assumes that sum of weights is 1.
        cumulative_prob = weights
        i = 1
        while i < len(cumulative_prob):
            cumulative_prob[i] += cumulative_prob[i-1]
            i += 1
        # x = random()
        x = self._random.random()
        i = 0
        cp = cumulative_prob[i]
        while x > cp:
            i += 1
            cp = cumulative_prob[i]
        
        return choices[i]
            
    def _deterministic_evolve(self, real_psi_i, t_i):
        # Evaluates the non-unitary evolution from initial state real_psi_i 
        # (expressed with real numbers with double the dimension to the complex wavector).
        # Assumes the initial state is at time t_i in the list of times self.t /
        i = np.where(np.isclose(self.t, t_i))[0][0]
        remainingtimes = self.t[i:]
        # realsoln = odeint(self._derivs, real_psi_i, remainingtimes, args = (self._rhmatrix,))
        realsoln = odeint(self._derivs, real_psi_i, remainingtimes)
        soln = self._reformatsolution(realsoln)
        return soln
        
        
    def _generate_nonhermitian_ham(self):
        # Adds the non-Hermitian terms to the Hamiltonian matrix to solve the
        # non-unitary time evoltuion
        nonh_term = np.zeros([self.dim, self.dim], complex)
        
        for i in self.qsys.lindbladgamma:
            gamma_i = self.qsys.lindbladgamma[i]
            adag_i = self.qsys.lindbladraising[i]
            a_i = self.qsys.lindbladlowering[i]
            newterm = -0.5 * gamma_i * adag_i @ a_i
            nonh_term += newterm
            
        self._rhmatrix += self._c2r(nonh_term)
        
    def _sqnorm(self, psi):
        # Returns the square norm of psi
        n2 = 0
        for c in psi:
            n2 += abs(c)**2
        return n2   
    
    
# test comment    
class unitaryevolution:
    """
    Unitary evolution.
    
    Finds the time evolution of the density matrix/operator by a unitary
    transformation. This method can be applied when the Hamiltonian has no
    explicit time dependence. The unitary evolution operator is found by
    diagonalising the Hamiltonian.
    
    Attributes
    ----------
    t : numpy.ndarray
        1d array of times to solve the system for.
    qsys : :class:`diracpy.quantum_systems.qsys`
        The qsys object that describes the quantum system to be solved
    rho0 : numpy.ndarray
        2d array representing the initial density matrix of the system in the
        basis given by qsys. If the initial density operator is given as
        a :class:`diracpy.states_operators.qop` object then the density
        matrix is generated using the `qsys.vectorize` method.
    soln : numpy.ndarray
        3d array containing the density matrices of the solved system for
        each time of the attribute `t`.
    """
    
    # z0 is the initial state of the system described by its density matrix.
    # z0 is represented by a 2d complex numpy array with the dimensions of the Hamiltonian
    # times is linearly spaced numpy array of times for which the system is solved for
    # This class is used to solve the dynamics for a static Hamiltonian. It uses the 
    # the Hamiltonian function (of time) in the ham_obj, evaluated at time zero.
    def __init__(self, rho0, times, qsys):
        """
        Construct unitary solver.
        
        Initialises the object with times, initial state (density operator)
        and the :class:`diracpy.quantum_systems.qsys` that describes the
        system.
        
        Parameters
        ----------
        rho0 : :class:`diracpy.states_operators.qop` or numpy.ndarray
            Initial density operator or density matrix of the system.
        times : numpy.ndarray
            1D array of times (`int` or `float`) to solve the system for.
        qsys : :class:`diracpy.quantum_systems.qsys`
            Quantum system to solve, which contains all information on the
            Hamiltonian and the basis used for the calculation.

        Returns
        -------
        None.
        """
        self.t = times
        self.set_initial_state(rho0)
        # self.z = z0
        # self.ham = ham_obj
        # self.dim = self.ham.dim
        # self.hmatrix = self.ham.ham(0)
        self._get_ham(qsys)
        
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
        
    def set_initial_state(self, rho0):
        """
        Define initial state.
        
        Define the initial state using the density matrix. When a 
        :class:`diracpy.states_operators.qop` is give, it is vecrtotized 
        by this method.

        Parameters
        ----------
        rho0 : :class:`diracpy.states_operators.qop` or numpy.ndarray
            Initial density operator or density matrix of the system.

        Returns
        -------
        None.
        """
        if type(rho0) == np.ndarray:
            rho = rho0
        else:
            rho = self.qsys.vectorize(rho0)
        self.rho0 = rho
        
    def _eigensolve(self):
        self.evals, self.evecs = np.linalg.eigh(self.hmatrix)
        
    def _u_op(self, t):
        u = np.zeros([self.dim, self.dim], complex)
        for i in range(self.dim):
            u[i,i] = np.exp(-1.j * self.evals[i] * t)
        u = self.evecs @ u @ self._hc(self.evecs)
        return u
        
    def _hc(self, np2darray):
        return np.transpose(np.conj(np2darray))
    
    def solve(self):
        """
        Solve system.
        
        Solves the system for array of uniformly spaces times, `t`, by 
        applying a unitary transformation to the density matrix at each
        time step.
        
        The solution is stored in the `soln` attribute of the object.

        Returns
        -------
        None.

        """
        
        self._eigensolve()
        self.soln = self._evolve(self.rho0)
            
    def _evolve(self, rho0):
        dt = self.t[1] - self.t[0]
        num_t = np.size(self.t)
        u_dt = self._u_op(dt)
        ud_dt = self._hc(u_dt)
        z = self.rho0
        soln = np.zeros([num_t, self.dim, self.dim], complex)
        for i, t in enumerate(self.t):
            soln[i] = z
            z = u_dt @ z @ ud_dt
        return soln
        
            
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
        
            
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    