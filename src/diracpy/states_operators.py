#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 15 15:33:45 2022

@author: benjaminyuen
"""

# Still To implement
# operator subtraction
# state subtraction
# superposition of qvecs formed using qvec or derived classes
# problem pickling __add__(), probably due to lambda functions.
# partial contraction of bras/kets with kets/bras of a factor space

from numpy import conj as cconj
from copy import deepcopy
#from diracpy.operators import qop

# need to include scalar division in qvec and qop

class qvec:
    def __init__(self, basisstate = None, cnum=1, **kwargs):
#        if kwargs.get('vector_type') == 'bra':
#            self.type = 'bra'
#        else:
#            self.type = 'ket'
        self.type = ''
        self.vec = {}
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
        return qvec(basisstate, cnum)
            
    def copy(self):
        return deepcopy(self)
    
    def c(self, basisstate):
        try:
            coefficient = self.vec[tuple(basisstate)]
        except KeyError:
            coefficient = 0
        return coefficient
    
    def update(self, basisstate, cnum=1):
        self.vec[tuple(basisstate)] = cnum
        
    def print(self):
        if self.vec:
            [print(value, ' * ', self.type, key) for key, value in self.vec.items()]
        else:
            print('')
            
    def __str__(self):
        output = ''
        for state, cnum in self.vec.items():
            if cnum == 0:
                pass
            elif cnum == 1:
                output += self.type + str(list(state)) + ' + '
            else:
                output += str(cnum) + ' * ' + self.type + str(list(state)) + ' + '
        output = output[:-3]
        return output
    
    def __mul__(self, scalar):
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
        if isinstance(scalar, (int, float, complex)):
            output = self.__mul__(scalar)
        else:
            output = NotImplemented
        return output
    
    def __div__(self, scalar):
        if isinstance(scalar,(int, float, complex)):
            basisstates = [state for state in self.vec]
            cnums = [cnum / scalar for cnum in self.vec.values()]
            output = self.new_instance(basisstate = basisstates, cnum = cnums)
        else:
            output = NotImplemented
        return output
    
    def __add__(self, other):
        if type(self) == type(other):
            new_qvec = self.copy()
            for other_bstate in other.vec:
                try:
                    new_qvec.vec[other_bstate] += other.vec[other_bstate]
                except KeyError:
                    new_qvec.vec[other_bstate] = other.vec[other_bstate]
            output = new_qvec
        else:
            output = NotImplemented
        return output
    
    def __radd__(self, other):
        return self.__add__(other)
    
    def __eq__(self, other):
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
            if other == 0:
                output = self.new_instance()
            else:
                basisstates = [state for state in self.vec]
                cnums = [other * cnum for cnum in self.vec.values()]
                output = self.new_instance(basisstate = basisstates, cnum = cnums)
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
    def __init__(self, action, conj_action = 0):
        # action should be a function which takes as input a ket and returns a ket
        # the action function should correspond to the linear operators action on a state vector
        self.action = action
        if conj_action == 0:
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
                return self.conj_action(ket_in) * other
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
            new_conj_action = lambda ket_in : self.conj_action(ket_in) * other
            qop_out = qop(new_action, new_conj_action)
            output = qop_out
        else:
            # no other right multiplications implemented
            # Right multiplication of a state vector of type bra covered by bra __mul__ method.
            output = NotImplemented
        return output
    
    def __div__(self, other):
        if isinstance(other, (int, float, complex)):
            # division of operator by a scalar
            new_action = lambda ket_in : self.action(ket_in) / other
            new_conj_action = lambda ket_in : self.conj_action(ket_in) / other
            qop_out = qop(new_action, new_conj_action)
            output = qop_out
        else:
            output = NotImplemented
        return output
    
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
    
    def conj(self):
        # Hermitian conjugation
        conj_op = qop(self.conj_action, self.action)
        return conj_op