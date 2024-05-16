"""Defines quantum state vectors 'bra' and 'ket' and quantum operators.

This module contains classes to define 'bra' and 'ket' objects as well as
'qop' (quantum operator) objects that act on these. These three classes
encode the required mathematical actions of such objects within a Hilbert
space. E.g. 'ket' objects satisy the axioms of a complex vector space, and
'bra' objects are there duals. Inner and outer products are defined. 'qop'
objects map 'ket' objects to new 'ket' objects, and 'bra' objects to new 'bra'
objects, and can be combined linearly as required.

Classes
-------
    qvec
    ket
    bra
    qop

Notes
-----
    - Objects cannot be pickled due to presence of lambda functions.
    - Partial contraction of bra and ket (e.g. for partial trace) not yet
    implemented.
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#C reated on Tue Feb 15 15:33:45 2022

# @author: benjaminyuen

# Still To implement
# problem pickling __add__(), probably due to lambda functions.
# partial contraction of bras/kets with kets/bras of a factor space.
# new_instance methods of qvec, ket and bra should be private, _new_instance.
# Constructor for qvec class should be rewritten, utilising private classes.

from numpy import conj as cconj
from copy import deepcopy
#from diracpy.operators import qop

#%% qvec class defintion
class qvec:
    """
    Parent class to define general properties of bra and ket vectors.
    
    This class defines attributes and methods common to both types of vectors.
    The bra and ket classes below are child classes of qvec.
    
    Attributes
    ----------
    type : str
        a string that will be either 'bra' or 'ket' in child classes
        (default '')
    vec : dict
        a dictionary whose keys are a list of basis state indices and values
        are the complex scalar coefficients of each basis state (default {}).
        The qvec represents the linear superposition of these basis states
        using the coefficients given here.
    """
    
    def __init__(self, basisstate = None, cnum=1, **kwargs):
        # **kwargs should go - but need to test with unit_tests.py first
        """
        Initialize a new instance of a qvec.
        
        Parameters
        ----------
        basisstate : list, tuple, string or a qvec, optional
            An iterable that gives the state index values of a basis state. 
            The default is None which yields the zero qvec. If type str given
            then state indices are each character.
        cnum : int, float or complex, or iterable of these, optional.
            The scalar multiple of the given basis state. The default is 1.
            If iterable, then these are the coeffients of each basis state.
        **kwargs : dict, optional
            Extra arguments to qvec.

        Returns
        -------
        None.
        """
        self.type = ''
        self.vec = {}
        # populate vec with basis states and cnum's
        try:
            # when basisstate is already a qvec, ket or bra
            self.vec = basisstate.vec
        except AttributeError:
            # when basisstate is an iterable or an iterable of iterables
            try:
                # when cnum is an interable corresponding to several state...
                for i, c in enumerate(cnum):
                    self.vec[ tuple( basisstate[i] ) ] = c
            except TypeError:
                # when cnum is not iterable...
                if basisstate != None: # if zero then empty state is created, i.e. self.vec= {}.
                    self.vec[ tuple( basisstate )] = cnum
                
    def new_instance(self, basisstate = None, cnum=1):
        """
        Create new instance of a qvec.
        
        Creates a new qvec with for the basis state given by basistate and
        the cnum specified. If non basis state is given then this produces 
        the zero qvec.
        
        Parameters
        ----------
        basisstate : list, tuple, str or qvec, optional
            The index values for the basis state, or a qvec. 
            The default is None.
        cnum : int, float, complex, optional
            The scalar multiplier of the given basis state. The default is 1.

        Returns
        -------
        qvec
            new qvec.
        """
        # This should be turned into a private method
        return qvec(basisstate, cnum)
            
    def copy(self):
        """
        Return a deepcopy of the qvec object.
        
        Returns
        -------
        qvec
            copy of self.
        """
        return deepcopy(self)
    
    def c(self, basisstate):
        """
        Return qvec coefficient.
        
        gets coefficient of specified basistate from the qvec object.
        If the basisstate is not in the qvec then zero is returned.

        Parameters
        ----------
        basisstate : list, tuple, string
            index values of the basisstate to get coefficient of.

        Returns
        -------
        coefficient : int, float or complex
            The coefficient of specified basistate
        """
        try:
            coefficient = self.vec[tuple(basisstate)]
        except KeyError:
            coefficient = 0
        return coefficient
    
    def update(self, basisstate, cnum=1):
        """
        Update the cnum for specified basistate with coefficient cnum.

        Parameters
        ----------
        basisstate : int, list, or tuple
            index values of basisstate.
        cnum : int, float or complex, optional
            Replacement coefficient value. The default is 1.

        Returns
        -------
        None.
        """
        self.vec[tuple(basisstate)] = cnum
        
    def print(self):
        """
        Print qvec coefficients and basisstates.

        Returns
        -------
        None.
        """
        if self.vec:
            [print(value, ' * ', self.type, key) for key, value in self.vec.items()]
        else:
            print(self.type+'[]')
            
    def __str__(self):
        """
        Implement str method for qvecs.

        Returns
        -------
        output : str
            returns string in format 
            cnum1 * basisstate1 + cnum2 * basisstate2 ...
        """
        output = ''
        if self.vec:
            for state, cnum in self.vec.items():
                if cnum == 0:
                    pass
                elif cnum == 1:
                    output += self.type + str(list(state)) + ' + '
                else:
                    output += str(cnum) + ' * ' + self.type + str(list(state)) + ' + '
            output = output[:-3]
        else:
            output = self.type+'[]'
        return output
    
    def __neg__(self):
        """
        Negation operator method.

        Returns
        -------
        qvec
            -1 * self.
        """
        return self.__mul__(-1)
    
    def __mul__(self, scalar):
        """
        Left multiplication.
        
        multiplication of qvecs from left by a scalar, scalar * qvec.
        multiplies all basisstates coefficients of this qvec by the scalar

        Parameters
        ----------
        scalar : int, float, or complex
            multiplier.

        Returns
        -------
        output : qvec
            multiplied qvec.
        """
        if isinstance(scalar, (int, float, complex)):
            if scalar == 0:
                output = self.new_instance()
            else:
                basisstates = [state for state in self.vec]
                cnums = [scalar * cnum for cnum in self.vec.values()]
                output = self.new_instance(basisstate = basisstates, cnum = cnums)
        else:
            output = NotImplemented
        return output
    
    def __rmul__(self, scalar):
        """
        Right multiplication.
        
        multiplication of qvecs from the right, qvec * scalar.

        Parameters
        ----------
        scalar : int, float or complex
            multiplier.

        Returns
        -------
        output : qvec
            right multiplied qvec.
        """
        if isinstance(scalar, (int, float, complex)):
            output = self.__mul__(scalar)
        else:
            output = NotImplemented
        return output
    
    def __truediv__(self, scalar):
        """
        Division.
        
        Division operator for qvec / scalar. Divides each basisstate
        coefficient by the scalar.

        Parameters
        ----------
        scalar : int, float, or complex
            divisor.

        Returns
        -------
        output : qvec
            divided qvec.
        """
        if isinstance(scalar,(int, float, complex)):
            basisstates = [state for state in self.vec]
            cnums = [cnum / scalar for cnum in self.vec.values()]
            output = self.new_instance(basisstate = basisstates, cnum = cnums)
        else:
            output = NotImplemented
        return output
    
    def __add__(self, other):
        """
        Left addition.
        
        Addition operator adds a qvec to another qvec of the same type.

        Parameters
        ----------
        other : qvec
            must be of same type, e.g. bra or ket.

        Returns
        -------
        output : qvec
            qvec of same type.
        """
        if type(self) == type(other):
            new_qvec = self.copy()
            for other_bstate in other.vec:
                try:
                    new_qvec.vec[other_bstate] += other.vec[other_bstate]
                except KeyError:
                    new_qvec.vec[other_bstate] = other.vec[other_bstate]
                if new_qvec.vec[other_bstate] == 0:
                    new_qvec.vec.pop(other_bstate)
            output = new_qvec
        else:
            output = NotImplemented
        return output
    
    def __radd__(self, other):
        """
        Right addition.
        
        Right addition operator for qvecs

        Parameters
        ----------
        other : qvec
            A qvec as same type.

        Returns
        -------
        qvec
            A qvec of same type.

        """
        return self.__add__(other)
    
    def __sub__(self, other):
        """
        Left subtraction.
        
        subtraction operators for two qvecs of the same type

        Parameters
        ----------
        other : qvec
            qvec of same type.

        Returns
        -------
        output : qvec
            qvec of same type.
        """
        if type(self) == type(other):
            new_qvec = self.copy()
            for other_bstate in other.vec:
                try:
                    new_qvec.vec[other_bstate] -= other.vec[other_bstate]
                except KeyError:
                    new_qvec.vec[other_bstate] = -other.vec[other_bstate]
                if new_qvec.vec[other_bstate] == 0:
                    new_qvec.vec.pop(other_bstate)
            output = new_qvec
        else:
            output = NotImplemented
        return output
        
    def __eq__(self, other):
        """
        Equality.
        
        equality operator for two vectors of the same type, i.e. qvec1 == qvec2
        Compares basisstate coefficients to see if qvecs are equal

        Parameters
        ----------
        other : qvec
            qvec of same type.

        Returns
        -------
        bool_out : bool
            True if qvec vec dictionaries are equal.
        """
        bool_out = (self.type == other.type and self.vec == other.vec)
        return bool_out

#%% ket class defintion
class ket(qvec):
    """
    Creates an object to represent a ket.
    
    This is a child class of qvec where many of the vector space and inner 
    operations on kets are defined. Methods that are specific to ket's are 
    defined here.
    
    Attributes
    ----------
    type : str
        value is 'ket'
    vec : dict
        inherited from qvec. All ket's are represented by a dictionary that
        gives the coefficients of each basis state that makes up the object.
        The the set of indices of each basisstate are the keys of this
        dictionary.
     
    """
    
    def __init__(self, basisstate = None, cnum=1):
        """
        Initialize a new instance of ket.
        
        Constructor for ket's inherits from qvecs. type is set to 'ket'.

        Parameters
        ----------
        basisstate : list, tuple, string or a qvec, optional
            indices or list of indices for basis state(s). The default is None.
        cnum : int, float or complex, or iterable of these, optional
            complex coefficients for basis states. The default is 1.

        Returns
        -------
        None.

        """
        super().__init__(basisstate, cnum)
        self.type = 'ket'
        
    def new_instance(self, basisstate = None, cnum=1):
        """
        Create a new instance of a ket.

        Parameters
        ----------
        basisstate : list, tuple, string or a qvec, optional
            indices or list of indices for basis state(s). The default is None.
        cnum : int, float or complex, or iterable of these, optional
            complex coefficients for basis states. The default is 1.

        Returns
        -------
        ket
            Returns new ket.

        """
        # This should be turned into a private method
        return ket(basisstate, cnum)
    
    def __mul__(self, other):
        """
        Left multiplication.
        
        Multiplication rule of a scalar or a bra by the ket.
        Scalar multiplication inherits fro qvec class. Multiplication a bra
        is the outer product which gives an operator (qop object)

        Parameters
        ----------
        other : int, float, complex, bra.
            multiplier.

        Returns
        -------
        output : ket or qop
            If multiplier is a scalar then a ket is returnes. If multiplier is
            a bra then a qop is returned.

        """
        if isinstance(other, (int, float, complex)):
            output = super().__mul__(other)
            # if other == 0:
            #     output = self.new_instance()
            # else:
            #     basisstates = [state for state in self.vec]
            #     cnums = [other * cnum for cnum in self.vec.values()]
            #     output = self.new_instance(basisstate = basisstates, cnum = cnums)
        elif isinstance(other, bra):
            # outer product between a ket and a bra
            action = lambda ket_in : (other * ket_in) * self
            conj_action = lambda ket_in : (self.conj() * ket_in) * other.conj()
            output = qop(action, conj_action)
        else:
            output = NotImplemented
        return output
        
    def __rmul__(self, other):
        """
        Right multiplication.
        
        Right multiplication be either a scalar or a bra. When the multiplier
        is a scalar then multiplication rule inherits from qvec. If multiplier
        is bra then the inner product is given, returning a scalar.

        Parameters
        ----------
        other : int, float, complex, bra
            multiplier.

        Returns
        -------
        output : ket or int, float, complex
            multiple of ket or inner product.

        """
        if isinstance(other, (int, float, complex)):
            output = super().__rmul__(other)
        elif isinstance(other, bra):
            output = 0
            for basis_bra, c_bra in other.vec.items():
                for basis_ket, c_ket in self.vec.items():
                    if basis_bra == basis_ket:
                        output += c_bra * c_ket
        else:
            output = NotImplemented
        return output
    
    def conj(self):
        """
        Hermitian conjugate.
        
        Defines conjugation of the ket - turning the ket into a bra

        Returns
        -------
        conj_bra : bra
            The conjugate bra to this ket.

        """
        basisstates = list(self.vec.keys())
        cnums = [cconj(value) for value in self.vec.values()]
        conj_bra = bra(basisstates, cnums)
        return conj_bra
    
#%% bra
class bra(qvec):
    """
    Creates an object to represent a ket.
    
    This is a child class of qvec where many of the vector space and inner 
    operations on kets are defined. Methods that are specific to bra's are 
    defined here.
    
    Attributes
    ----------
    type : str
        value is 'bra'
    vec : dict
        inherited from qvec. All ket's are represented by a dictionary that
        gives the coefficients of each basis state that makes up the object.
        The the set of indices of each basisstate are the keys of this
        dictionary.
        
    """
    
    def __init__(self, basisstate = None, cnum=1):
        """
        Initialize a new instance of bra.

        Parameters
        ----------
        basisstate : list, tuple, string or a qvec, optional
            indices or list of indices for basis state(s). The default is None.
        cnum : int, float or complex, or iterable of these, optional
            complex coefficients for basis states. The default is 1.

        Returns
        -------
        None.

        """
        super().__init__(basisstate, cnum)
        self.type = 'bra'
        
    def new_instance(self, basisstate = None, cnum=1):
        """
        Return new instance of a bra.

        Parameters
        ----------
        basisstate : list, tuple, string or a qvec, optional
        cnum : int, float or complex, or iterable of these, optional
        
        Returns
        -------
        bra
            bra object for given basisstate and cnum.

        """
        # This should be turned into a private method
        return bra(basisstate, cnum)
        
    def __mul__(self, other):
        """
        Left multiplication.

        Parameters
        ----------
        other : int, float, complex, ket, qop
            The multiplier.

        Returns
        -------
        output : bra, int, float, complex, 
            if other is a scalar or a qop then return type is bra
            if other is ket then a scalar is returned
        """
        if isinstance(other, (int, float, complex)):
            output = super().__mul__(other)
        elif isinstance(other, ket):
            output = 0
            for basis_bra, c_bra in other.vec.items():
                for basis_ket, c_ket in self.vec.items():
                    if basis_bra == basis_ket:
                        output += c_bra * c_ket
        elif isinstance(other, qop):
            output_conj = other.conj() * self.conj()
            output = output_conj.conj()
        else:
            output = NotImplemented
        return output
    
    def conj(self):
        """
        Hermitian conjugate.

        Returns
        -------
        conj_ket : ket
            The conjugate ket to this bra.

        """
        basisstates = list(self.vec.keys())
        cnums = [cconj(value) for value in self.vec.values()]
        conj_ket = ket(basisstates, cnums)
        return conj_ket
    
class qop:
    def __init__(self, action, conj_action = None):
        # action should be a function which takes as input a ket and returns a ket
        # the action function should correspond to the linear operators action on a state vector
        self.action = action
        if conj_action == None:
            # if no conj_action give, assume the operator is Hermitian
            self.conj_action = self.action
        else:
            self.conj_action = conj_action
            
#    def check_empty_ket(self, action, ket_in):
#        if bool(ket_in.vec):
#            ket_out = action(ket_in)
#        else:
#            ket_out = ket()
#        return ket_out
            
        
    def __mul__(self, other):
        if isinstance(other, ket):
            # Operator action on state vectors
            ket_out = ket()
            for bstate, c in other.vec.items():
                _state_in = ket(bstate, c)
                ket_out = ket_out + self.action(_state_in)
            output = ket_out
        elif isinstance(other, (int, float, complex)):
            # multiplication of operator with a scalar
#            new_action = lambda ket_in : self.action(ket_in) * other
#            new_conj_action = lambda ket_in : self.conj_action(ket_in) * other
            def new_action(ket_in):
                return self.action(ket_in) * other
            def new_conj_action(ket_in):
                return self.conj_action(ket_in) * cconj(other)
            qop_out = qop(new_action, new_conj_action)
            output = qop_out
        elif isinstance(other, (qop,)):
            # multiplication of two operators
            new_action = lambda ket_in : self.action( other.action(ket_in) )
            new_conj_action = lambda ket_in : other.conj_action(self.conj_action(ket_in))
            qop_out = qop(new_action, new_conj_action)
            output = qop_out
        else:
            # multiplication by any other object not implemented
            output = NotImplemented
        return output
        
    def __rmul__(self, other):
        if isinstance(other, (int, float, complex)):
            # multiplication scalar * operator
            new_action = lambda ket_in : self.action(ket_in) * other
            new_conj_action = lambda ket_in : self.conj_action(ket_in) * cconj(other)
            qop_out = qop(new_action, new_conj_action)
            output = qop_out
        else:
            # no other right multiplications implemented
            # Right multiplication of a state vector of type bra covered by bra __mul__ method.
            output = NotImplemented
        return output
    
    def __truediv__(self, other):
        if isinstance(other, (int, float, complex)):
            # division of operator by a scalar
            new_action = lambda ket_in : self.action(ket_in) / other
            new_conj_action = lambda ket_in : self.conj_action(ket_in) / cconj(other)
            qop_out = qop(new_action, new_conj_action)
            output = qop_out
        else:
            output = NotImplemented
        return output

    # test comment
    def __add__(self, other):
        if isinstance(other, qop):
            # Addition of two operators
            new_action = lambda ket_in : self.action(ket_in) + other.action(ket_in)
            new_conj_action = lambda ket_in : self.conj_action(ket_in) + other.conj_action(ket_in)
            qop_out = qop(new_action, new_conj_action)
            output = qop_out
        else:
            # Only two operators may be added. All other types not implemented
            output = NotImplemented
        return output
    
    def __sub__(self, other):
        if isinstance(other, qop):
            # Addition of two operators
            new_action = lambda ket_in : self.action(ket_in) - other.action(ket_in)
            new_conj_action = lambda ket_in : self.conj_action(ket_in) - other.conj_action(ket_in)
            qop_out = qop(new_action, new_conj_action)
            output = qop_out
        else:
            # Only two operators may be added. All other types not implemented
            output = NotImplemented
        return output
    
    def __neg__(self):
        return self.__mul__(-1)
    
    def conj(self):
        # Hermitian conjugation
        conj_op = qop(self.conj_action, self.action)
        return conj_op