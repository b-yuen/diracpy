#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 15 15:33:45 2022

@author: benjaminyuen
"""

# Still To implement
# problem pickling __add__(), probably due to lambda functions.
# partial contraction of bras/kets with kets/bras of a factor space

from numpy import conj as cconj
from copy import deepcopy
#from diracpy.operators import qop


class qvec:
    """
    Parent class to define general properties of bra and ket vectors
    
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
        
    Methods
    -------
    new_instance(basisstate=None, cnum=1)
        returns a new qvec object
    copy()
        returns a deepcopy of the qvec object
    c(basistate)
        gives the coefficient of the given basis state for this qvec
    update(basisstate, cnum=1)
        updates the coefficent of the given basis state for this qvec
    print()
        prints each component of the qvec
    __str__()
        defines output of str(qvec) and print(qvec)
    __neg__()
        defines negation of a qvec
    __mul__(scalar)
        defines scalar multiplication of a qvec
    __rmul__(scalar)
        defines multiplication by a scalar from the right
    __truediv__(scalar)
        defines scalar devision
    __add__(other) 
        defines addition qvec addition from left
    __radd__(other)
        defines addition qvec addition from the right
    __sub__(other)
        defines subtraction of qvecs
    __eq__(other)
        defines equality qvecs
    """
    
    def __init__(self, basisstate = None, cnum=1, **kwargs):
        # **kwargs should go - but need to test with unit_tests.py first
        """
        Parameters
        ----------
        basisstate : list, tuple, string or a qvec, optional
            An iterable that gives the state index values of a basis state. 
            The default is None which yields the zero qvec. If type str given
            then state indices are each character.
        cnum : int, float or complex, optional
            The scalar multiple of the given basis state. The default is 1.
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
        Creates a new qvec with for the basis state given by basistate and
        the cnum specified.
        
        If non basis state is given then this produces the zero qvec.
        
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
        
        return qvec(basisstate, cnum)
            
    def copy(self):
        """
        Returns
        -------
        qvec
            copy of self.
        """
        
        return deepcopy(self)
    
    def c(self, basisstate):
        """
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
        update the cnum for specified basistate with coefficient cnum.

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
        prints qvec coefficients and basisstates

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
        str method for qvecs

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
        negation operator method

        Returns
        -------
        qvec
            -1 * self.

        """
        return self.__mul__(-1)
    
    def __mul__(self, scalar):
        """
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
        Addition operator to adds a qvec to another qvec of the same type.

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

class ket(qvec):
    def __init__(self, basisstate = None, cnum=1):
        super().__init__(basisstate, cnum)
        self.type = 'ket'
        
    def new_instance(self, basisstate = None, cnum=1):
        return ket(basisstate, cnum)
    
    def __mul__(self, other):
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
        basisstates = list(self.vec.keys())
        cnums = [cconj(value) for value in self.vec.values()]
        conj_bra = bra(basisstates, cnums)
        return conj_bra
    
class bra(qvec):
    def __init__(self, basisstate = None, cnum=1):
        super().__init__(basisstate, cnum)
        self.type = 'bra'
        
    def new_instance(self, basisstate = None, cnum=1):
        return bra(basisstate, cnum)
        
    def __mul__(self, other):
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