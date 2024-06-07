#!/usr/bin/env python3
# -*- coding: utf-8 -*-

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
# Error handling for when incorrect types are given (not that magic methods
# require NotImplemented to be returned instead, as coded below).

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
                if basisstate != None and len(tuple(basisstate)) > 0:
                    self.vec[ tuple( basisstate )] = cnum
                else:
                    # zero qvec is created, i.e. self.vec= {}
                    pass
                
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
    
    This is a child class of :class:`qvec` where many of the vector space and inner 
    operations on kets are defined. Methods that are specific to ket's are 
    defined here.
    
    Attributes
    ----------
    type : str
        value is 'ket'
    vec : dict
        inherited from :class:`qvec`. All ket's are represented by a dictionary that
        gives the coefficients of each basis state that makes up the object.
        The the set of indices of each basisstate are the keys of this
        dictionary.
     
    """
    
    def __init__(self, basisstate = None, cnum=1):
        """
        Initialize a new instance of ket.
        
        Constructor for ket's inherits from :class:`qvecs`. type is set to 'ket'.

        Parameters
        ----------
        basisstate : list, tuple, string or a :class:`qvec`, optional
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
        Scalar multiplication inherits from :class:`qvec` class. Multiplication
        with a :class:`bra` is the outer product which gives an operator 
        (:class:`qop` object)

        Parameters
        ----------
        other : int, float, complex, bra.
            multiplier.

        Returns
        -------
        output : :class:`ket` or :class:`qop` object.
            If multiplier is a scalar then a ket is returnes. If multiplier is
            a :class:`bra` then a :class:`qop` is returned.

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
        
        Right multiplication be either a scalar or a :class:`bra`. When the multiplier
        is a scalar then multiplication rule inherits from :class:`qvec`. If multiplier
        is :class:`bra` then the inner product is given, returning a scalar.

        Parameters
        ----------
        other : int, float, complex, :class:`bra`
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
        
        Defines conjugation of the :class:`ket` object - turning the ket into 
        a bra

        Returns
        -------
        conj_bra : :class:`bra`
            The conjugate bra to this ket.

        """
        basisstates = list(self.vec.keys())
        cnums = [cconj(value) for value in self.vec.values()]
        conj_bra = bra(basisstates, cnums)
        return conj_bra
    
#%% bra
class bra(qvec):
    """
    Creates an object to represent a bra.
    
    This is a child class of :class:`qvec` where many of the vector space and 
    inner operations on bra's are defined. Methods that are specific to bra's 
    are defined here.
    
    Attributes
    ----------
    type : str
        value is 'bra'
    vec : dict
        inherited from :class:`qvec`. All ket's are represented by a dictionary that
        gives the coefficients of each basis state that makes up the object.
        The the set of indices of each basisstate are the keys of this
        dictionary.
        
    """
    
    def __init__(self, basisstate = None, cnum=1):
        """
        Initialize a new instance of bra.

        Parameters
        ----------
        basisstate : list, tuple, string or a :class:`qvec`, optional
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
        Return new instance of a :class:`bra`.

        Parameters
        ----------
        basisstate : list, tuple, string or a :class:`qvec` object, optional
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
    """
    Quantum operators.
    
    Defines quantum operator (qop) objects. qop objects act on ket and bra
    objects using __mul__ and __rmul__ magic mathods (i.e. * operator). qop
    objects can also be added and subtracted from each other, multiplied or
    divided by scalars and multiplied by other qop objects. All such operations
    form new qop objects.
    
    Attributes
    ----------
    action : callable[[ket], ket]
        Function that defines the logical operation of the operator on kets.
    conj_action : callable[[]ket, ket], optional
        Function that defines logical operation of Hermitian conjugate of the
        operator on kets. If not defined then operator assumed to be Hermitian.
    """
    
    def __init__(self, action, conj_action = None):
        """
        Inititialize quantum operator.
        
        Quantum operators act on kets, mapping them onto a new ket and complex
        coefficient. The operation on kets is defined by the action function
        that must be passed on initialization. The action function must 
        therefore be a function whose argument is a ket and returns a ket and
        corresponds to the linear operatos action on a state vector. The
        action function must be able to determine the output ket, together with
        any multiplicative factors, from the basis index values of the input 
        ket. The conjugate action works in the same way, but defines how the
        Hermitian conjugate of this operator acts on kets. The conjugate action
        is used to determine how this operator acts on bras to its left.

        Parameters
        ----------
        action : callable[[ket], ket]
            Function that defines the logical operation of the operator on 
            kets.
        conj_action : callable[[]ket, ket], optional
            Function that defines action of the Hermitian conjugate of this
            operator

        Returns
        -------
        None.

        """
        self.action = action
        if conj_action == None:
            # if no conj_action given, assume the operator is Hermitian
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
        """
        Left multiplication.
        
        Defines multiplication from the left on a scalar, ket, or operator.
        When other is a scalar, then a new action is defined that includes
        multiplication by this scalar.

        Parameters
        ----------
        other : int, float, complex, ket, qop.
            scalar, ket or qop the operator acts on.
            If other is a scalar, then a new action is defined that includes
            multiplication by this scalar.
            If other is a ket, then action is applied to determine the ket that
            should be returned.
            If other is an operator, then a new composite action is defined.

        Returns
        -------
        output : qop, ket
            If other is a scalar or qop, then a qop is returned
            If other is a ket then a ket is returned.
        """
        if isinstance(other, ket):
            # Operator action on state vectors
            ket_out = ket()
            for bstate, c in other.vec.items():
                _state_in = ket(bstate, c)
                ket_out = ket_out + self.action(_state_in)
            output = ket_out
        elif isinstance(other, (int, float, complex)):
            # multiplication of operator with a scalar
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
        """
        Right multiplication.
        
        Defines multiplication to the right, on a scalar. This must be defined
        since a scalars __mul__ method will not include instances where other
        is a qop. Right multiplication on bra and qop are determined from
        other.__mul__ .

        Parameters
        ----------
        other : int, float, complex
            Scalar multiplier.

        Returns
        -------
        output : qop
            A new qop with action modified to include multiplication by given
            scalar.
        """
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
        """
        Scalar division.
        
        Defines of an operator division by a scalar, by modifying the operator
        action function appropriately.

        Parameters
        ----------
        other : int, float, complex
            Scalar divisor.

        Returns
        -------
        output : qop
            A qop with modified action (and conj_action) that accounts for the
            specified scalar division.

        """
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
        """
        Left Addtion.
        
        Addition operation between an operator and another operator to its 
        right.

        Parameters
        ----------
        other : qop
            The quantum operator to add to.

        Returns
        -------
        output : qop
            A new operator whose action is the sum of the outputs of 
            self.action and other.action.
        """
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
        """
        Left subtraction.
        
        Defines operator subtraction of another qop from this qop.

        Parameters
        ----------
        other : qop
            qop to subtract.

        Returns
        -------
        output : qop
            new qop with action that returns the output of this operators
            action minus the other operators action.
        """
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
        """
        Negation.
        
        Defines the negative of this qop.

        Returns
        -------
        qop
            A new qop that is -1 * this qop.

        """
        return self.__mul__(-1)
    
    def conj(self):
        """
        Hermition conjugate.
        
        Defines Hemrition conjugation of this operator.

        Returns
        -------
        conj_op : qop
            A new qop with action and conj_action swaped.
        """
        # Hermitian conjugation
        conj_op = qop(self.conj_action, self.action)
        return conj_op
    
    
    
    
    
    