#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""Builds quantum systems.

This module contains classes to build a quantum system from a Hamiltonian,
an initial state or set of states, and an optional set of Lindblad operators.
This module provides an interface between the system definition in diracpy
notation using :class:`bra's <diracpy.states_operators.bra>`,
:class:`kets's <diracpy.states_operators.ket>` and
:class:`qop's <diracpy.states_operators.qop>`, and the solvers of the
quantum dyanmics in the :mod:`diracpy.quantum_dynamics`.

Classes
-------
    qsys
    qsys_t

Notes
-----
Diracpy allows one to define the system through its Hamiltonian abstractly,
in terms of operators, rather than needing to specify a basis and Hamiltonian
matrix. Similary, states of the system (such as the initial state), can be
entered immediately within the context of a state space; the complete basis
of the Hilbert space does not need to be explicitly defined, only the
index values for the state at hand. However, when it comes to solving the
system dynamics, almost all methods required the problem to be vectorised
within an explicitly defined basis.

This module handles the vectorisation problem automatically, given an
an initial state or set of states. All states that interact with these
initial states up to :math:`n^{\mathrm {th}}` order in the given
Hamiltonian are included (see for example chapter 3 of 
:cite:`cohen1998atom`
for a detailed discussion of truncation of the state space by interaction
Hamiltonian order). 
If Lindblad operators are given too, the basis is
extended accordingly.

Examples
--------
To build a :class:`quantum system <qsys>` we first define the Hamiltonian.
For the Jaynes-Cummings model for example, where

.. math::
    
    H = \omega a^{\dagger} a + \omega_0 \sigma_z +
    g(a \sigma_+ + a^{\dagger} \sigma_-)
    
we use

>>> atom = dp.two_level_subspace(index=0)
>>> cavity = dp.fock_subspace(index=1)
>>> omega, omega_0, g = 100, 100, np.pi
>>> H_0 = omega * cavity.n + omega_0 * atom.sigma_z
>>> V = g * (cavity.a * atom.sigma_plus + cavity.adag * atom.sigma_minus)
>>> H = H_0 + V

where we have chosen :math:`\omega=\omega_0=100` and
:math:`g=\pi`. We furthermore choose the intial state
:math:`\vert \psi_0 \rangle = \vert e,0 \rangle` such that the atom
is excited and the cavity is in vacuum.

>>> psi0 = dp.ket(['e',0])

We can now build the quantum system,

>>> system = dp.qsys(hamiltonian_operator=H, initialstates=[psi0], n_int=2)

We immediately have access to the Hamiltonian matrix

>>> system.print_ham()
     50+0j  3.14159+0j  
3.14159+0j       50+0j

and can list the basis that the Hamiltonian matrix is specified in

>>> [print(state) for state in system.basis];
ket['e', 0]
ket['g', 1]

"""

# Created on Sat Jan 15 21:31:18 2022
# @author: benjaminyuen

from diracpy.states_operators import ket
from diracpy.states_operators import bra
from diracpy.states_operators import qop
import numpy as np
import time

class qsys:
    r"""Create quantum system for a static Hamiltonian.
    
    This class constructs quantum system objects for time-independent
    Hamiltonians. The quantum system consists of the Hamiltonian,
    any Lindblad operators, a basis for the Hilbert space, and the
    Hamiltonian's matrix in this basis.
    
    Attributes
    ----------
    n_int : int
        The interaction order to which the basis is found
    ham_op : :class:`diracpy.states_operators.qop`, list or 1d np.array of qop
        The hamiltonian used to define system. This can be a single
        :class:`diracpy.states_operators.qop` object, or a list or 1d np.array
        of :class:`diracpy.states_operators.qop` objects
        corresponding to terms that comprise the Hamiltonian.
    jump_ops : list or 1d np.array of :class:`diracpy.states_operators.qop` objects
        Jump operators should specify the lindblad lowering term or terms
        used in the master equation only. The square root of the decay
        coefficient should be included in this term.
    basis : list of :classL`diracpy.states_operators.ket`.
        Basis used for vecrtorisation of the quantum systems Hilbert space.
        All subsequent vectors and matices will be in this basis with same
        component indexing as 'basis'.
    adjoint_basis : list of :classL`diracpy.states_operators.bra`.
        The basis adjoint to :attr:`basis`.
    dim : int
        Dimension of the Hilbert space of the 'qsys' instance.
    hmatrix : numpy.ndarray
        Numpy array of the Hamiltonian matrix in the :attr:`basis`.
    lindbladlowering : dict
        Dictionary of lindblad lowering operators of type 
        :class:`diracpy.states_operators.qop`. These are automatically
        generated when a list of :attr:`jump_ops` are specified on
        intsantiation of a qsys object. Alternatively, this
        dictionary can be manually populated after the object is initialised.
        Used in :class:`quantum_dynamics.lindblad`, 
        :class:`quantum_dynamics.quantum_jumps`, and
        :class:`quantum_dynamics.liouvillian`.
    lindbladraising : dict
        Dictionary of lindblad raising operators of type 
        :class:`diracpy.states_operators.qop`. These are the Hermitian
        conjugates of the lindbladlowering terms. These are automatically
        generated when a list of :attr:`jump_ops` are specified on
        intsantiation of a qsys object. Alternatively, this
        dictionary can be manually populated after the object is initialised.
        Used in :class:`quantum_dynamics.lindblad`, 
        :class:`quantum_dynamics.quantum_jumps`, and
        :class:`quantum_dynamics.liouvillian`.
    lindbladgamma : dict 
        Dictionary of coefficients (int, float, complex) of the lindbland 
        operators to apply in the master equation. When generated 
        automatically, the coefficient is 1 for each term since the
        coefficient is combined with :attr:`jump_ops`. When Lindblad
        terms are specified manually after initialisation, it is permitted
        that coefficients different from 1 can be specified here.
        
    """
    
    def __init__(self, hamiltonian_operator=None, H0terms=None, Vterms=None,
                 jump_ops=None, n_int=0, initialstates=None, **kwargs):
        """
        Initialise qsys object.
        
        Initialisation comprises of parsing input hamiltonian (which could be
        an array of terms) and quantum jumps operator(s), building the 
        Hilbert space basis and Hamiltonian matrix in this basis, and defining
        lindblad terms if included.

        Parameters
        ----------
        hamiltonian_operator : :class:`diracpy.states_operators.qop`, optional.
            When the Hamiltonian is represted by a single 
            :class:`diracpy.states_operators.qop` object it passed here. 
            The default is None, in which case the terms should be stated
            by 'H0terms' and 'Vterms'
        H0terms : list or 1D np.array, optional.
            List or 1D np.array of :class:`diracpy.states_operators.qop` 
            objects containing the terms of :math:`H_0` 
            The default is None.
        Vterms : list or 1D np.array, optional.
            List or 1D np.array of :class:`diracpy.states_operators.qop` 
            objects containing the terms of the interaction Hamiltonian 
            :math:`V`. The default is None.
        jump_ops : :class:`diracpy.states_operators.qop`, list or array, optional.
            Open system quantum jump operator or if there is more than one,
            a list or array of quantum jump operators for the system.
            The default is None.
        n_int : int, optional
            Interaction order used when calculating extent of basis defined by
            :func:`build_sys`
            The default is 0.
        initialstates : list of :class:`diracpy.states_operators.ket` objects, optional.
            The initialstate or a set of states to define the basis from.
            See :func:`build_system` for details. Other iterables other than
            list can also be used. The default is None.
        **kwargs : dict
            Opional named argument 'verbose' of type bool. Default is False.
            When 'True' progress is reported during object construction which
            is helpful for extremly large systems.

        Raises
        ------
        TypeError
            Raised when `H0terms` and `Vterms` are given 
            but elements are not of type
            :class:`diracpy.states_operators.qop`.
        NameError
            At least one of `hamiltonian_operator`, 
            `H0terms` or `Vterms` should be given, otherwise
            NameError is raised.
        NameError
            Raised when no `initialstates` are given.
            
        Returns
        -------
        None.

        """
        # hamiltonian operator should be of type dp.qop
        #comment
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
        
        # **kwargs includes keyword argument verbose, which defaults to false
        # this is used to print messages on build progress
        
        
        self.n_int = n_int
        # parse hamiltonians
        self._input_hamiltonian(hamiltonian_operator, H0terms, Vterms)
        # parse jump operators
        self.jump_ops = self._input_jump_ops(jump_ops)
        self.build_sys(initialstates, **kwargs)
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
        
    # Checks at least one of input arguments hamiltonian_operator, H0terms or Vterms are populated.
    # If populated check hamiltonian is a qop or H0terms and/or V0terms are list or array of qops.
    def _input_hamiltonian(self, hamiltonian_operator, H0terms, Vterms):
        if isinstance(hamiltonian_operator, qop):
            self.ham_op = hamiltonian_operator
            self._ham_array_flag = False
        else:
            self.H0terms = []
            self.Vterms = []
            if H0terms == None:
                pass
            elif isinstance(H0terms[0], qop):
                self.H0terms = H0terms
                self._ham_array_flag = True
            else:
                raise TypeError("H0terms should be a list or array of type diracpy.states_operators.qop")
            if Vterms == None:
                pass
            elif isinstance(Vterms[0], qop):
                self.Vterms = Vterms
                self._ham_array_flag = True
            else:
                raise TypeError("Vterms should be a list or array of type diracpy.states_operators.qop")
        try:
            self._ham_array_flag
        except:
            NameError("No hamiltonian has been specified. At least one of hamiltonian_operator, H0terms or Vterms must be given.")
            
    def build_sys(self, initial_states, **kwargs):
        """
        Build basis and Hamiltonian matrix attributes of the qsys object.
        
        Builds the basis first, then the Hamiltonian matrix in this basis.
        For `n_int=0` the basis is comprised of states in `initialstates`.
        For `n_int=1` the basis the `hamiltonian_operator` is applied to each
        state in `initialstates`, and the resultant states are appended to 
        the `initialstates` to give a wider basis.
        Generally, this is repeated `n_int` times to construct a basis of
        states that coherently interact with the `initialstates` up to
        order `n_int` in the `hamilonian_operator`. 
        When `jump_ops` are specified, then each jump_op is applied to each
        of the coherently interacting states. This is repeated until no
        further states are found.

        Parameters
        ----------
        initial_states : list of :class:`diracpy.states_operators.ket` objects, optional.
            The initialstate or a set of states to define the basis from.
            See :func:`build_system` for details. Other iterables other than
            list can also be used. The default is None.
        **kwargs : dict
            Optional named arguement `verbose` which gives information on
            progress of build. Default value is `verbose=False`.

        Raises
        ------
        NameError
            Raised when no `initialstates` are given.

        Returns
        -------
        None.

        """
        # print(initial_states)
        if initial_states == None:
            raise NameError("System not built, no initial basis states given.")
        else:
            self._build(initial_states, **kwargs)
        
    def _build(self, initialstates, build_hmatrix=True, verbose=False):
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
        if verbose: print("building sys basis...")
        t1 = time.time()
        
        self.basis = self._build_coherent_sys(initialstates)
        
        # add states these decay to
        if len(self.jump_ops) == 0:
            basis_incomplete = False
        else:
            basis_incomplete = True
        while basis_incomplete:
            if verbose: print('...appending decay subspaces...')
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
        if verbose: print("...system basis built in {} seconds".format(t2-t1))

        self.adjoint_basis = [bra(state) for state in self.basis]
        self.dim = len(self.basis)
        if build_hmatrix==True:
            if verbose: print("defining hmatrix...")
            t3 = time.time()
            # self.hmatrix = self.make_hmatrix()
            self.make_hmatrix()
            t4 = time.time()
            if verbose: print("...hmatrix evaluated in {} seconds".format(t4-t3))
        
    # builds the basis states coherently coupled to initialstates via hamiltonian  
    def _build_coherent_sys(self, initialstates):
        # build sys when hamiltonian given via arrays H0terms and or Vterms
        if self._ham_array_flag == True:
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
        """
        Build Hamiltonian matrix attribute.
        
        Builds the hamiltonian :class:`diracpy.states_operators` operator
        into a matrix in the basis constructed for the qsys instance.
        The hamiltonian matrix is assigned to the attribute
        :attr:`self.hmatrix <hmatrix>`.

        Returns
        -------
        None.

        """
        if self._ham_array_flag == True:
            matrix_evaluator = self._hmatrix_from_ham_array
        else:
            matrix_evaluator = self._hmatrix_from_ham_op
        self.hmatrix = matrix_evaluator()
        # return matrix_evaluator()
        
    def _hmatrix_from_ham_op(self):
        hmatrix = self.matrix(self.ham_op)
        return hmatrix
            
    def _hmatrix_from_ham_array(self):
        hmatrix = np.zeros([self.dim, self.dim], complex)
        for term_op in self.H0terms:
            term_matrix = self.matrix(term_op)
            hmatrix += term_matrix
        for term_op in self.Vterms:
            term_matrix = self.matrix(term_op)
            hmatrix += term_matrix
        return hmatrix
            
    def matrix(self, operator):
        """
        Build matrix of operator.
        
        Returns a matrix for the given :class:`diracpy.states_operators` 
        operator in the basis constructed for the qsys instance.

        Parameters
        ----------
        operator : :class:`diracpy.states_operators`
            Operator for which matrix representation is found.

        Raises
        ------
        TypeError
            Raised when input is not of type :class:`diracpy.states_operators`.

        Returns
        -------
        matrix_out : numpy.ndarray
            Matrix representation of the input operator.

        """
        if not isinstance(operator, qop):
            raise TypeError("first positional argument should be of type diracpy.states_operators.qop")
        matrix_out = np.zeros([self.dim, self.dim], complex)
        for j, basis_ket in enumerate(self.basis):
            _psi = operator * basis_ket
            for i, basis_bra in enumerate(self.adjoint_basis):
                matrix_out[i,j] = basis_bra * _psi
        return matrix_out
    
    def print_basis(self):
        """
        Print basis states.
        
        Prints the list of basis states defined for qsys instance in
        a readable format.

        Returns
        -------
        None.
        """
        [state.print() for state in self.basis]
        
    def print_ham(self):
        """
        Print hamiltonian matrix.
        
        Prints the hamiltonians hmatrix representation in
        a readable format.

        Returns
        -------
        None.
        """
        self.matprint(self.hmatrix)
                
    def matprint(self, mat, fmt="g"):
        """
        Print matrix.
        
        Prints input matrix 'mat' in
        a readable format.

        Parameters
        ----------
        mat : input matrix
            Matrix to be printed.
        fmt : str, optional
            Format specifier for matrix elements. The default is "g" such
            the the string representation of each matrix element is
            "{:g}".format(<matrix element>).
            See `Format String Syntax <https://docs.python.org/3/library/string.html#formatstrings>`_
            for other formats.

        Returns
        -------
        None.

        """
        col_maxes = [max([len(("{:"+fmt+"}").format(x)) for x in col]) for col in mat.T]
        for x in mat:
            for i, y in enumerate(x):
                print(("{:"+str(col_maxes[i])+fmt+"}").format(y), end="  ")
            print("")
                
    def ham(self, t):
        """
        Return Hamiltonian matrix as a function of t.
        
        Returns hmatrix for this class since the hamiltonian is static.
        This method is defined so that the functionality is extendable to
        time-dependent Hamiltonians.

        Parameters
        ----------
        t : int, float
            Time.

        Returns
        -------
        numpy.ndarray
            The :attr:`hmatrix <hmatrix>`.

        """
        return self.hmatrix
    
    def add_jump_op(self, jump_qop):
        """
        Append jump operator.
        
        Appends `jump_op` to list of jump operators, :attr:`jump_ops`.

        Parameters
        ----------
        jump_qop : :class:`diracpy.states_operators.qop`
            The new jump operator to append.

        Raises
        ------
        TypeError
            Raised when `jump_qop` is not of type 
            :class:`diracpy.states_operators.qop`.

        Returns
        -------
        None.

        """
        if not isinstance(jump_qop, qop):
            raise TypeError("first positional argument should be of type diracpy.states_operators.qop")
        # jump_qop should be the lowering operator associated with the quantum jump
        self.jump_ops.append(jump_qop)
        
    def make_lindblads(self):
        """
        Make lindblad dictionaries.
        
        Defines or redefines the :attr:`lindbladgamma`, :attr:`lindbladraising` 
        and :attr:`lindbladlowering` from :attr:`jump_ops`.
        

        Returns
        -------
        None.
        
        """
        self.lindbladgamma = {}
        self.lindbladraising = {}
        self.lindbladlowering = {}
        for index, jump_qop in enumerate(self.jump_ops):
            self.lindbladgamma[index] = 1
            self.lindbladlowering[index] = self.matrix(jump_qop)
            self.lindbladraising[index] = self.matrix(jump_qop.conj())
            
class qsys_t(qsys):
    r"""Create quantum system for a time-dependent Hamiltonian.
    
    This class constructs quantum system objects for time-dependent
    Hamiltonians. The quantum system consists of the Hamiltonian,
    any Lindblad operators, a basis for the Hilbert space, and the
    Hamiltonian's matrix in this basis. This class inherits from
    :class:`qsys`, and extends its functionality to define a time
    dependent hamiltonian matrix.
    
    Attributes
    ----------
    n_int : int
        The interaction order to which the basis is found
    ham_op : :class:`diracpy.states_operators.qop`, list or 1d np.array of qop
        The hamiltonian used to define system. This can be a single
        :class:`diracpy.states_operators.qop` object, or a list or 1d np.array
        of :class:`diracpy.states_operators.qop` objects
        corresponding to terms that comprise the Hamiltonian.
    jump_ops : list of 1d np.array of :class:`diracpy.states_operators.qop` 
    objects
        Jump operators should specify the lindblad lowering term or terms
        used in the master equation only. The square root of the decay
        coefficient should be included in this term.
    basis : list of :classL`diracpy.states_operators.ket`.
        Basis used for vecrtorisation of the quantum systems Hilbert space.
        All subsequent vectors and matices will be in this basis with same
        component indexing as 'basis'.
    adjoint_basis : list of :classL`diracpy.states_operators.bra`.
        The basis adjoint to :attr:`basis`.
    dim : int
        Dimension of the Hilbert space of the 'qsys' instance.
    lindbladlowering : dict
        Dictionary of lindblad lowering operators of type 
        :class:`diracpy.states_operators.qop`. These are automatically
        generated when a list of :attr:`jump_ops` are specified on
        intsantiation of a qsys object. Alternatively, this
        dictionary can be manually populated after the object is initialised.
        Used in :class:`quantum_dynamics.lindblad`, 
        :class:`quantum_dynamics.quantum_jumps`, and
        :class:`quantum_dynamics.liouvillian`.
    lindbladraising : dict
        Dictionary of lindblad raising operators of type 
        :class:`diracpy.states_operators.qop`. These are the Hermitian
        conjugates of the lindbladlowering terms. These are automatically
        generated when a list of :attr:`jump_ops` are specified on
        intsantiation of a qsys object. Alternatively, this
        dictionary can be manually populated after the object is initialised.
        Used in :class:`quantum_dynamics.lindblad`, 
        :class:`quantum_dynamics.quantum_jumps`, and
        :class:`quantum_dynamics.liouvillian`.
    lindbladgamma : dict 
        Dictionary of coefficients (int, float, complex) of the lindbland 
        operators to apply in the master equation. When generated 
        automatically, the coefficient is 1 for each term since the
        coefficient is combined with :attr:`jump_ops`. When Lindblad
        terms are specified manually after initialisation, it is permitted
        that coefficients different from 1 can be specified here.
    
    """
    
    def __init__(self, static_operator, dynamic_operators, 
                 dynamic_coefficients, **kwargs):
        """
        Initialise qsys_t object.
        
        Initialisation comprises of parsing input hamiltonian (which could be
        an array of terms) and quantum jumps operator(s), building the 
        Hilbert space basis, and defining lindblad terms if included.
        Furthermore, the method :func:`ham` which returns hamiltonian matrix
        as a function of time is initialised from the parameters below.

        Parameters
        ----------
        static_operator : :class:`diracpy.states_operators.qop`.
            Operator describing the static part of the Hamiltonian.
        dynamic_operators : :class:`diracpy.states_operators.qop` or list of qop objects.
            The terms of the Hamiltonian matrix that have dyanmic coeffcients.
            The terms themselves are static - time dependence comes from
            the time dependent factors 'dynamic_coefficients'
        dynamic_coefficients : function or list of function objects.
            Dynamic coefficients of the time dependent part of the Hamiltonian.
        **kwargs : dict
            See list of paramters for :class:`qsys`. From these the
            'hamiltonian_operator', 'H0terms' and 'Vterms' are reduntant in
            for this derived class, but all others optional named arguments
            of :class:`qsys` are valid.

        Returns
        -------
        None.

        """
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
        self._static_operator = static_operator
        self._dynamic_operators = dynamic_operators
        self._dynamic_coefficients = dynamic_coefficients
        self.static_hamiltonian = self._static_operator
        try:
            for op in self._dynamic_operators:
                self.static_hamiltonian = self.static_hamiltonian + op
        except TypeError:
            self.static_hamiltonian = self.static_hamiltonian + self._dynamic_operators
        super().__init__(hamiltonian_operator=self.static_hamiltonian, **kwargs)
        self._ham_matrices()
        
    def _ham_matrices(self):
        self._static_op_matrix = self.matrix(self._static_operator)
        try:
            self._dynamic_op_matrices = [self.matrix(op) for op in self._dynamic_operators]
        except TypeError:
            self._dynamic_op_matrices = [self.matrix( self._dynamic_operators )]
            self._dynamic_coefficients = [self._dynamic_coefficients]
        # Check 'dynamic_coefficients' are correctly specified
        try:
            _testval = 0
            for func in self._dynamic_coefficients:
                _testval += func(0)
        except TypeError:
            raise TypeError("'dynamic_coefficients' should be a scalar function or list of functions of 1 positional argument")
        if len(self._dynamic_coefficients) == len(self._dynamic_op_matrices):
            self._num_dynamic_components = len(self._dynamic_coefficients)
        else:
            raise IndexError("'dynamic_operators' and 'dynamic_coefficients' inputs should be same length")
            
        
    def ham(self, t):
        """
        Time dependent Hamiltonian.
        
        Calculates and returns the matrix for the time dependent
        Hamiltonian at time t.

        Parameters
        ----------
        t : int, float
            time.

        Returns
        -------
        current_hmatrix : numpy.ndarray
            matrix of the time dependent Hamiltonian at time t.

        """
        current_hmatrix = self._static_op_matrix.copy()
        for i in range(self._num_dynamic_components):
            current_hmatrix += self._dynamic_coefficients[i](t) * self._dynamic_op_matrices[i]
        return current_hmatrix
            
        
        
        
        
        
        
        

        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
    
