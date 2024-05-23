#!/usr/bin/env python3
# -*- coding: utf-8 -*-

r"""Defines quantum state spaces and their operators.

This module contains the definition of the :class:`subspace` class, and a
number of classes that inherit from it to define some specific state spaces
together with their operators. These include :class:`two_level_subspace` for
defining the state space of a two level system together with a complete set of
operators, and :class:`fock_subspace` that defines the space of Fock states
together with bosonic number, creation and annihilation operators. See
:ref:`Classes <subspace_classes>` for a full list of subspace classes defined
in this module. User's can define further subspace classes, which should
inherit from :class:`subspace`.

.. _subspace_classes:
    
Classes
-------

    subspace
    two_level_subspace
    three_level_subspace
    fock_subspace
    floquet_subspace
    
Notes
-----
    - State spaces (Hilbert spaces) are the building block of quantum systems, 
      and are spanned by a set of linearly indepedent (typically orthonormal) 
      states. In dirac notation these basis states are distinguished by the 
      value of a designated index. For example, the Fock space is given the 
      index 'n' which can take any integer value, and spin half particle is be 
      given the index 'm' which can take values -1/2 or 1/2. Tensor product 
      spaces can then be formed, where it is clear which index refers to which 
      state space, e.g. :math:`\vert n,m \rangle` label states with the boson 
      number as the first 
      index and the spin by the second index. To reflect this convention in 
      diracypy the position of the state index is an attribute that is common 
      to all statespace classes, i.e. :class:`subspace` or any class that 
      inherits from it, which are all classes within this module. Tensor 
      product spaces are then easily formed by defining new 'subspace'
      objects with a different index position.
    - Once a (possibly infinite) set of state labels is defined for a given 
      subspace, then it is straight forward to define the action of operators 
      on this statespace, e.g. the bosonic creation operator raises the value 
      of the index  by 1, and multiplies the state by the square root of the 
      new index value. In this lies the power of defining state spaces in this 
      abstract way that is natural to dirac notation -- rather than using 
      explicit column vectors for example. Once the state spaces index
      values are at least abstractly defined, then all operators acting on the
      state space can be defined also, using a fixed (and finite) amount of 
      memory.
      
    .. _op_actions:

    - Operators are defined by their actions. In the derived subspace classes
      some specific operators have been defined. To define these, one first
      must define their 'action' on an individual basis state, i.e. one that
      is not a superposition state. This action is then generalised, for 
      example to be able to map over over superposition states
      or to act on other operators to form composite operators,
      by wrapping it in the :class:`diracpy.states_operators.qop class`.
    - Subspace is a misnomer. The term state_space would here be accurate, and
      should be updated in a future version, since each 'subspace' is a state
      space on its own, but when a tensor product space is formed, then it
      is not neccesarily a subspace of the tensor product space.
      
Examples
--------
First, a simple example of a subspace:

>>> import diracpy as dp
>>> tls = dp.two_level_subspace()
>>> psi = tls.sigma_plus * dp.ket(['g'])
>>> print(psi)
ket['e']

This defines a two level system, formed by the state space
:math:`\{ \vert g \rangle, \vert e \rangle \}`, together with the
operators :math:`\sigma_z,\ \sigma_+` and :math:`\sigma_-`. Note here that
the position of the index is by default set to 0. In this example we compute

.. math:: \sigma_+ \vert g \rangle

Next, demonstrate how we can easily extend this to a tensor product of two 
such subpaces:
    
>>> tls1 = dp.two_level_subspace(0)
>>> tls2 = dp.two_level_subspace(1)
>>> psi = tls1.sigma_plus * tls2.sigma_minus * dp.ket(['g','e'])
>>> print(psi)
ket['e','g']

This defines the tensor product state of a pair of two level subspace. We
then calculate

.. math:: \sigma^{(1)}_+ \sigma^{(2)}_+ \vert g,e \rangle

where :math:`\sigma^{(1)}_+` acts on the first two level subspace and
:math:`\sigma^{(2)}_+` acts on the second.

The next example illustrates the need for the operator action methods of the
classes derived from :class:`subspace`. This example is of use to developers
wanting to implement additional subspace classes with different operators than
those defined in the :mod:`diracpy.subspaces` module.
When implementing an operator in diracpy we first define its 'operator action' 
on a suitable basis of the state space on which it acts. For example, the
spin operators :math:`\sigma_z, \ \sigma_+` and :math:`\sigma_-` are defined
in the :class:`two_level_subspace` class. To do this, we have defined a basis
by specifying the index values 'g' and 'e' on which the operators act;

>>> tls.basis
[('e',), ('g',)]

Having defined a statespace basis we can define the action on operators in
terms of this basis. For example, we expect 
:math:`\sigma_z \vert g \rangle = -1/2 \vert g \rangle` and
:math:`\sigma_z \vert e \rangle = 1/2 \vert e \rangle`. The action method
:meth:`two_level_subspace.sigma_z_action` describes this mapping, by matching
the value of the subspace instances index to either 'e' or 'g' and returning
'0.5 * ket_in' for 'e' and '-0.5 * ket_in' for 'g';

>>> ket_e, ket_g = dp.ket(['e']), dp.ket(['g'])
>>> print(tls.sigma_z_action(ket_e))
>>> print(tls.sigma_z_action(ket_g))
0.5 * ket['e']
-0.5 * ket['g']

Because the tls subspace has 'index=0' it matches only the first index to these
values. E.g. the same result is obtained for

>>> psi_1, psi_2 = dp.ket(['e',1,4,'x']), dp.ket(['g',2,3,'y'])
>>> print(tls.sigma_z_action(psi_1))
>>> print(tls.sigma_z_action(psi_2))
0.5 * ket['e', 1, 4, 'x']
-0.5 * ket['g', 2, 3, 'y']

It is straight forward to define the action on these particular basis kets
in this way, but what about generalisation to superposition states? This is
done automatically by turning the action method into a 
:class:diracpy.states_operators.qop object, which inherently defines all the 
rules of how such an operator should behave, e.g. on superposition states,
on other operators, under scalar multiplication, etc. Hence the method
:meth:`two_level_subspace.sigma_z_action` is used to define the
:attr:`two_level_subspace.sigma_z` which is an attribute of the
:class:`two_level_subspace` class. The code to do this, defined in
:meth:`two_level_subspace.__init__` , is simply

>>> self.sigma_z = qop(self.sigma_z_action)

where the default conjugate method of 'self.sigma_z_action' is used since
:math:`\sigma_z` is Hermitian. Finally, we have the operator
:meth:`two_level_subspace.sigma_z_action` define for each subspace instance.
E.g.

>>> print(tls.sigma_z * (ket_e + ket_g))
0.5 * ket['e'] + -0.5 * ket['g']
>>> psi = dp.ket(['e','g'])
>>> print(tls1.sigma_z * tls2.sigma_z * psi)
-0.25 * ket['e', 'g']

The operators :attr:`two_level_subspace.sigma_plus` and
:attr:`two_level_subspace.sigma_minus` are defined in a similar way using
the actions methods :attr:`two_level_subspace.sigma_plus_action` and
:attr:`two_level_subspace.sigma_minus_action`.


"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on Tue Feb 15 15:44:22 2022
# @author: benjaminyuen

# Still to implement/change
# replace **kwargs with index=0 in subspace.__init__()

from diracpy.states_operators import ket
#from diracpy.states import bra
from diracpy.states_operators import qop
from numpy import sqrt

class subspace:
    """
    Parent class to define shared properties of subspaces.
    
    All subspaces require an index which designates the position of the state
    label for the subspace. This is a standard convention one automatically
    chooses when writing states in dirac notation, but must be specified in 
    diracpy.
    
    Attributes
    ----------
    index : int
        The index of the state label for this subspace. The default is 0.
    basis : list
        An empty list to (optionally) contain basis states for the susbpace.
        
    """
    
    def __init__(self, index=0):
        """
        Construct subspace object.
        
        Constructor for subspace class.

        Parameters
        ----------
        index : int
            index of state label. Default 0.

        Returns
        -------
        None.

        """
        self.index = index
        self.basis = []
        
    def modify_bstate(self, bstate, new_label):
        """
        Modify index value of subspace instance's index.
        
        Operators in derived subspace classes will modify the index
        of the subspace instance.
        State vectors (:class:`diracpy.states_operators.bra` and 
        :class:`diracpy.states_operators.ket`) can be labelled by list 
        (or iterable) of multiple indices. This method selects and modifies
        the index for a particular subspace instance.
        
        Parameters
        ----------
        bstate : list, tuple or string.
            An iterable that gives the state index values of a basis state.
        new_label : int or string
            The new value of the state index. By convention this should
            be of same type as the old index value, i.e. bstate[self.index].

        Returns
        -------
        list, tuple or string.
            The modified iterable of state indices.

        """
        # produces a copy of a basis state tuple with the element at index self.index changed to new_label
        # modifies bstate made by concatenation
        return bstate[:self.index] + (new_label,) + bstate[self.index+1:]
        
        
class two_level_subspace(subspace):
    r"""
    Create two level system subspace object.
    
    Creates a subspace object the represents a two level system, together with
    a complete set (i.e. Lie algebra closes) of operators that act on a two
    level system. This class inherits from :class:`subspace`.
    
    Attributes
    ----------
    index : int
        Inherited from :class:`subspace`.
    basis : list of tuple
        Specifies the set of basis states.
        A list of tuples, where each tuple contains the basis state index.
        In this class the indices take string values 'e' or 'g'.
    sigma_z_dict : dict
        Dictionary containing 'spin projections' for the states 'e' and 'g'.
        Default is {'e':0.5, 'g':-0.5}.
        This can be modified after construction 
        (see :ref:`example <tls_example1>`).
    sigma_z : :class:`diracpy.states_operators.qop`
        The spin projection operator, :math:`\sigma_z`.
    sigma_plus : :class:`diracpy.states_operators.qop`
        The raising operator, :math:`\sigma_+`.
    sigma_minus : :class:`diracpy.states_operators.qop`
        The lowering operator, :math:`\sigma_-`.
    
    Examples
    --------
    .. _tls_example1:
    
    >>> tls = dp.two_level_subspace()
    >>> states = [dp.ket(['e']), dp.ket(['g'])]
    >>> [print(tls.sigma_z * state) for state in states];
    0.5 * ket['e']
    -0.5 * ket['g']
    
    Now modify :attr:`sigma_z_dict` so that the ground state has energy zero and
    the excited state energy is 1.
    
    >>> tls.sigma_z_dict = {'e':1, 'g':0}
    >>> [print(tls.sigma_z * state) for state in states];
    ket['e']
    ket[]
    
    """
    
    # need to build index into ket referencing so that it can handle multi-indexed basis states
    def __init__(self, **kwargs):
        """
        Construct object two_level_subspace object.

        Construct object to represent a two level system state space, spanned
        by `{diracpy.ket(['e']), diracpy.ket(['g'])}`, together with operators
        on this state space.

        Parameters
        ----------
        **kwargs : dict, optional.
            Optional named argument `index` to specify :attr:`index`. Default
            value of `index` is zero.

        Returns
        -------
        None.

        """
        super().__init__(**kwargs)
        self.basis = [('e',),('g',)]
        self.sigma_z_dict = {'e':0.5, 'g':-0.5}
        self.sigma_z = qop(self.sigma_z_action)
        self.sigma_plus = qop(self.sigma_plus_action, self.sigma_minus_action)
        self.sigma_minus = qop(self.sigma_minus_action, self.sigma_plus_action)

    def sigma_z_action(self, ket_in):
        r"""
        Define action for :attr:`sigma_z` operator.
        
        Defines mapping of ket_in onto ket_out under the action of the
        :math:`\sigma_z` operator where ket_in and ket_out are individual
        basis states, i.e. not superposition states.
        If ket_in is in state 'g' (for the subspace instance index),
        then the scalar multiplied is :attr:`sigma_z_dict`.'g',
        which has default value `-0.5`. If ket_in is in state 'e', then
        then the scalar multiplier is :attr:`sigma_z_dict`.'e',
        which has default value `0.5`. 
        See :ref:`note <op_actions>`  on operator action
        for general description.

        Parameters
        ----------
        ket_in : :class:`diracpy.states_operators.ket`
            A ket corresponding to an indivual basis state. It can be a tensor
            product state, but not a superposition.

        Returns
        -------
        ket_out : :class:`diracpy.states_operators.ket`
            output state :math:`\sigma_z \vert \psi \rangle`.

        """
        # assumes that the input ket_in is for a single basis state but with arbitrary cnum.
        try:
            basisstate_in = list(ket_in.vec.keys())[0] # state index values
            c_in = ket_in.vec[basisstate_in]
            c_out = self.sigma_z_dict[basisstate_in[self.index]] * c_in
            ket_out = ket(basisstate_in, c_out)
        except IndexError:
            ket_out = ket()
        return ket_out
    
    def sigma_minus_action(self, ket_in):
        r"""
        Define action for :attr:`sigma_minus` operator.
        
        This action defines the mapping 
        :math:`\vert g\rangle \mapsto \vert e \rangle`,
        :math:`\vert e \rangle \mapsto 0`.

        Parameters
        ----------
        ket_in : :class:`diracpy.states_operators.ket`
            Input ket.

        Returns
        -------
        ket_out : :class:`diracpy.states_operators.ket`
            Output ket.

        """
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
        r"""
        Define action for :attr:`sigma_plus` operator.
        
        This action defines the mapping 
        :math:`\vert g\rangle 0`,
        :math:`\vert e \rangle \mapsto \vert g \rangle`.

        Parameters
        ----------
        ket_in : :class:`diracpy.states_operators.ket`
            Input ket.

        Returns
        -------
        ket_out : :class:`diracpy.states_operators.ket`
            Output ket.

        """
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
    r"""
    Create three level statespace and associated operators.
    
    Creates subspace objects that describe a three level lambda system with
    two ground states 
    :math:`\vert g_1 \rangle` and :math:`\vert g_2 \rangle`
    and one shared excited state :math:`\vert e \rangle`.
    Inherits from :class:`subspace` class.
    
    Attributes
    ----------
    basis : list of tuples.
        Define the basis state index values for this state space class.
    sigma_z_dict : dict
        Defines a lookup table for values returned by the 'spin projection'
        like operator :attr:`sigma_z`. This can be modified after construction
        in the same way as the two level system
        :ref:`example <tls_example1>`.
    sigma_z : :class:`diracpy.states_operators.qop`
        Spin projection operator. Maps basis states back to themselves
        and multiplies them by the coefficients defined by 
        :attr:`sigma_z_dict`.
    sigma_minus1 : :class:`diracpy.states_operators.qop`
        The lowering operator :math:`\vert g_1 \rangle \langle e \vert`
    sigma_plu1 : :class:`diracpy.states_operators.qop`
        The raising operator :math:`\vert e \rangle \langle g_1 \vert`
    sigma_minus2 : :class:`diracpy.states_operators.qop`
        The lowering operator :math:`\vert g_2 \rangle \langle e \vert`
    sigma_plus1 : :class:`diracpy.states_operators.qop`
        The raising operator :math:`\vert e \rangle \langle g_2 \vert`
    """
    
    def __init__(self, **kwargs):
        """
        Construct state space.
        
        Constructs the three level system state space spanned by
        `dp.ket(['g_1'])`, `dp.ket(['g_2'])`, and `dp.ket(['e'])`, together
        with associated operators.

        Parameters
        ----------
        **kwargs : dict
            optional named argument index of type int that specifies the index
            position for state space instance. Defaults to 0.

        Returns
        -------
        None.
        """
        super().__init__(**kwargs)
        self.basis = [('e',),('g1',),('g2',)]
        self.sigma_z_dict = {'g1':-0.5, 'g2':0.5, 'e':0}
        self.sigma_z = qop(self.sigma_z_action)
        self.sigma_minus1 = qop(self.sigma_minus1_action, self.sigma_plus1_action)
        self.sigma_plus1 = qop(self.sigma_plus1_action, self.sigma_minus1_action)
        self.sigma_minus2 = qop(self.sigma_minus2_action, self.sigma_plus2_action)
        self.sigma_plus2 = qop(self.sigma_plus2_action, self.sigma_minus2_action)
    
    def sigma_minus1_action(self, ket_in):
        r"""
        Define action of :attr:`sigma_minus1`.
        
        Define lowering operator action 
        :math:`\vert g_1 \rangle \langle e \vert` on basis states defined.

        Parameters
        ----------
        ket_in : :class:`diracpy.states_operators.ket`
            Input ket.

        Returns
        -------
        ket_out : :class:`diracpy.states_operators.ket`
            Output ket.
        """
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
        r"""
        Define action of :attr:`sigma_plus1`.
        
        Define raising operator action 
        :math:`\vert e \rangle \langle g1 \vert` on basis states defined.

        Parameters
        ----------
        ket_in : :class:`diracpy.states_operators.ket`
            Input ket.

        Returns
        -------
        ket_out : :class:`diracpy.states_operators.ket`
            Output ket.
        """
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
        r"""
        Define action of :attr:`sigma_minus2`.
        
        Define lowering operator action 
        :math:`\vert e \rangle \langle g2 \vert` on basis states defined.

        Parameters
        ----------
        ket_in : :class:`diracpy.states_operators.ket`
            Input ket.

        Returns
        -------
        ket_out : :class:`diracpy.states_operators.ket`
            Output ket.
        """
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
        r"""
        Define action of :attr:`sigma_plus2`.
        
        Define raising operator action 
        :math:`\vert e \rangle \langle g2 \vert` on basis states defined.

        Parameters
        ----------
        ket_in : :class:`diracpy.states_operators.ket`
            Input ket.

        Returns
        -------
        ket_out : :class:`diracpy.states_operators.ket`
            Output ket.
        """
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
        r"""
        Define action of :attr:`sigma_z`.
        
        Define 'spin projection' operator action which with the default
        coefficients of :attr:`sigma_z_dict` takes the form
        :math:`0 \vert e \rangle \langle e \vert - 0.5 \vert g_1 \rangle 
        \langle g_1 \vert + 0.5 \vert g_2 \rangle \langle g_2 \vert`.

        Parameters
        ----------
        ket_in : :class:`diracpy.states_operators.ket`
            Input ket.

        Returns
        -------
        ket_out : :class:`diracpy.states_operators.ket`
            Output ket.
        """
        try:
            basisstate_in = list(ket_in.vec.keys())[0]
            key = basisstate_in[self.index]
            ket_out = self.sigma_z_dict[key] * ket_in
        except IndexError:
            ket_out = ket()
        return ket_out
        
class fock_subspace(subspace):
    r"""
    Create Fock state space and associated operators.
    
    Creates subspace objects that describe Fock space spanned by
    :math:`\vert n \rangle \forall n \in \mathbb N`. The number, creation
    and annihilation operators 
    :math:`N, \ a^{\dagger}` and :math:`a`, are also defined. Inherits
    from :class:`subspace` class.
    
    Attributes
    ----------
    basis : list of tuples.
        Defaults to an empty list since the infinte dimensional state space
        is defined abstractly.
    a : :class:`diracpy.states_operators.qop`
        Bosonic annihilation operator.
    adag : :class:`diracpy.states_operators.qop`
        Bosonic creation operator.
    n : :class:`diracpy.states_operators.qop`
        Bosonic number operator equivalent to 'adag * a'.
    """
    
    def __init__(self, **kwargs):
        """
        Construct Fock space instance.
        
        Creates a Fock state space instance, together with opertors.

        Parameters
        ----------
        **kwargs : dict
            Optional named arguments of 'index', the state space's 
            index position which defaults to 0. If 'max_n' is specified then a
            list of basis states in made up to 'max_n'. The need for this
            is superceded by the automatic basis state generation perfomed in
            :mod:`diracpy.quantum_systems`.

        Returns
        -------
        None.

        """
        super().__init__(**kwargs)
#        self.index = 0
        if 'max_n' in kwargs:
            self.basis = [(n,) for n in range(kwargs['max_n']+1)]
        self.a = qop(self.a_action, self.adag_action)
        self.adag = qop(self.adag_action, self.a_action)
        self.n = qop(self.n_action)
        
    def a_action(self, ket_in):
        r"""
        Define action of :attr:`a`.
        
        Define annihilation operator action which takes the form
        :math:`a \vert n \rangle = \sqrt n \vert n-1 \rangle \forall n>0`
        and `a \vert 0 \rangle = 0`.

        Parameters
        ----------
        ket_in : :class:`diracpy.states_operators.ket`
            Input ket.

        Returns
        -------
        ket_out : :class:`diracpy.states_operators.ket`
            Output ket.
        """
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
        r"""
        Define action of :attr:`adag`.
        
        Define creation operator action which takes the form
        :math:`a^{\dagger} \vert n \rangle = \sqrt{n+1} \vert n+1 \rangle 
        \forall n>=0`.


        Parameters
        ----------
        ket_in : :class:`diracpy.states_operators.ket`
            Input ket.

        Returns
        -------
        ket_out : :class:`diracpy.states_operators.ket`
            Output ket.
        """
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
        r"""
        Define action of :attr:`n`.
        
        Define number operator action which takes the form
        :math:`N \vert n \rangle = n \vert n \rangle`.

        Parameters
        ----------
        ket_in : :class:`diracpy.states_operators.ket`
            Input ket.

        Returns
        -------
        ket_out : :class:`diracpy.states_operators.ket`
            Output ket.
        """
        try:
            basisstate_in = list(ket_in.vec.keys())[0]
            n_in = basisstate_in[self.index]
            ket_out = n_in * ket_in
        except IndexError:
            ket_out = ket()
        return ket_out
    
class floquet_subspace(subspace):
    r"""
    Create Floquet state space and associated operators.
    
    Creates subspace objects that describe the state space
    :math:`\{ \vert n \rangle \vert \forall n \in \mathbb Z \}`, together with
    the number, raising and lowering operators,
    :math:`N, \ a^{\dagger}` and :math:`a`. This is used to desribe
    highly excited bosonic states where :math:`\sqrt{n+1}\approx \sqrt{n}`.
    
    Attributes
    ----------
    basis : list of tuples.
        Defaults to an empty list since the infinte dimensional state space
        is defined abstractly.
    a : :class:`diracpy.states_operators.qop`
        Lowering operator.
    adag : :class:`diracpy.states_operators.qop`
        Raising operator.
    n : :class:`diracpy.states_operators.qop`
        Number operator. Essential since not equivalent to
        :math:`a^{dagger} a`.
    """
    
    def __init__(self, **kwargs):
        """
        Construct Floquet space instance.
        
        Creates a Floquet state space instance, together with opertors.

        Parameters
        ----------
        **kwargs : dict
            Optional named arguments of 'index', the state space's 
            index position which defaults to 0. If 'max_n' is specified then a
            list of basis states in made up to 'max_n'. The need for this
            is superceded by the automatic basis state generation perfomed in
            :mod:`diracpy.quantum_systems`.

        Returns
        -------
        None.

        """
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
        r"""
        Define action of :attr:`a`.
        
        Define lowering operator action which takes the form
        :math:`a \vert n \rangle = \vert n-1 \rangle \forall n>\mathbb Z`.

        Parameters
        ----------
        ket_in : :class:`diracpy.states_operators.ket`
            Input ket.

        Returns
        -------
        ket_out : :class:`diracpy.states_operators.ket`
            Output ket.
        """
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
        r"""
        Define action of :attr:`adag`.
        
        Define raising operator action which takes the form
        :math:`a \vert n \rangle = \vert n+1 \rangle \forall n>\mathbb Z`.

        Parameters
        ----------
        ket_in : :class:`diracpy.states_operators.ket`
            Input ket.

        Returns
        -------
        ket_out : :class:`diracpy.states_operators.ket`
            Output ket.
        """
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
        r"""
        Define action of :attr:`n`.
        
        Define number operator action which takes the form
        :math:`N \vert n \rangle = n \vert n \rangle`.

        Parameters
        ----------
        ket_in : :class:`diracpy.states_operators.ket`
            Input ket.

        Returns
        -------
        ket_out : :class:`diracpy.states_operators.ket`
            Output ket.
        """
        try:
            basisstate_in = list(ket_in.vec.keys())[0]
            n_in = basisstate_in[self.index]
            ket_out = n_in * ket_in
        except IndexError:
            ket_out = ket()
        return ket_out