#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import diracpy as dp


class unittest_base:
    def __init__(self):
        self.empty_ket = dp.ket([])
        self.empty_bra = dp.bra([])
        
class unittest_base_vector(unittest_base):
    def __init__(self):
        unittest_base.__init__(self)
        self.bra_u = dp.bra([1])
        self.bra_v = dp.bra([2])
        self.bra_w = dp.bra([3])
        self.ket_u = dp.ket([1])
        self.ket_v = dp.ket([2])
        self.ket_w = dp.ket([3])
        self.bra_id_add = self.empty_bra
        self.ket_id_add = self.empty_ket
        
    def com_test_bra(self):
        LHS = self.bra_u + self.bra_v
        RHS = self.bra_v + self.bra_u
        return LHS == RHS
    
    def com_test_ket(self):
        LHS = self.ket_u + self.ket_v
        RHS = self.ket_v + self.ket_u
        return LHS == RHS
   
    def assoc_test_bra(self):
        LHS = self.bra_u + (self.bra_v + self.bra_w)
        RHS = (self.bra_u + self.bra_v) + self.bra_w
        return LHS == RHS
 
    def assoc_test_ket(self):
        LHS = self.ket_u + (self.ket_v + self.ket_w)
        RHS = (self.ket_u + self.ket_v) + self.ket_w
        return LHS == RHS

    def id_test_bra(self):
        LHS = self.bra_u + self.bra_id_add
        RHS = self.bra_u
        return LHS == RHS

    def id_test_ket(self):
        LHS = self.ket_u + self.ket_id_add
        RHS = self.ket_u
        return LHS == RHS
    
    def inv_test_bra(self):
        LHS = self.bra_u + (-self.bra_u)
        RHS = self.bra_id_add
        return LHS == RHS
    
    def inv_test_ket(self):
        LHS = self.ket_u + (-self.ket_u)
        RHS = self.ket_id_add
        return LHS == RHS  
    
    def sca_mul_test_bra(self):
        LHS = (2)**(1/2) *((8)**(1/2)*self.bra_u)
        RHS = ((2)**(1/2) *(8)**(1/2))*self.bra_u
        return LHS == RHS

    def sca_mul_test_ket(self):
        LHS = (2)**(1/2) *((8)**(1/2)*self.ket_u)
        RHS = ((2)**(1/2) *(8)**(1/2))*self.ket_u
        return LHS == RHS
        
    def id_sca_mul_test_bra(self):
        LHS = 1* self.bra_u
        RHS = self.bra_u
        return LHS == RHS
        
    def id_sca_mul_test_ket(self):
        LHS = 1* self.ket_u
        RHS = self.ket_u
        return LHS == RHS
    
    def dist_vec_add_test_bra(self):
        LHS = (2)**(1/2) * (self.bra_u + self.bra_v)
        RHS = (2)**(1/2) * self.bra_u + (2)**(1/2) * self.bra_v
        return LHS == RHS
    
    def dist_vec_add_test_ket(self):
        LHS = (2)**(1/2) * (self.ket_u + self.ket_v)
        RHS = (2)**(1/2) * self.ket_u + (2)**(1/2) * self.ket_v
        return LHS == RHS
    
    def dist_sca_add_test_bra(self):
        LHS = ((2)**(1/2) + 2) * self.bra_u
        RHS = (2)**(1/2) * self.bra_u + 2 * self.bra_u
        return LHS == RHS
    
    def dist_sca_add_test_ket(self):
        LHS = ((2)**(1/2) + 2) * self.ket_u
        RHS = (2)**(1/2) * self.ket_u + 2 * self.ket_u
        return LHS == RHS
    
class unittest_vector(unittest_base_vector):
    def __init__(self):
        super().__init__()
        
    def unit_test_vec(self, end_test = False):
        com_result, com_message = self.com_test()
        asoc_result, asoc_message = self.asoc_test()
        id_result, id_message = self.id_test()
        inv_result, inv_message = self.inv_test()
        sca_mul_result, sca_mul_message = self.sca_mul_test()
        id_scal_mul_result, id_scal_mul_message = self.id_scal_mul_test()
        dist_sca_add_result, dist_sca_add_message = self.dist_sca_add_test()
        dist_vec_add_result, dist_vec_add_message = self.dist_vec_add_test()
        
        total_results = [com_result, asoc_result, id_result, inv_result, sca_mul_result, 
                         id_scal_mul_result, dist_sca_add_result, dist_vec_add_result]  
        total_message = np.array([com_message, asoc_message, id_message, inv_message, sca_mul_message, 
                         id_scal_mul_message, dist_sca_add_message, dist_vec_add_message])
        
        if np.prod(total_results):
            return True, "All Vector Space Tests Passed"
        else:
            if not end_test:
                return False, total_message[np.invert(total_results)]
            else:
                output = total_message[np.invert(total_results)].tolist()
                print(*output, sep="\n")
        
    def com_test(self):
        bra = self.com_test_bra()
        ket = self.com_test_ket()
        total = bra * ket
        if total:
            return True, "Pass"
        elif bra == ket:
            return False, "Vector Space Commutivity Fail:both bra and ket"
        elif not bra:
            return False, "Vector Space Commutivity Fail:bra"
        elif not ket:
            return False, "Vector Space Commutivity Fail:ket"
    
    def asoc_test(self):
        bra = self.assoc_test_bra()
        ket = self.assoc_test_ket()
        total = bra * ket
        if total:
            return True, "Pass"
        elif bra == ket:
            return False, "Vector Space Associativity Fail:both bra and ket"
        elif not bra:
            return False, "Vector Space Associativity Fail:bra"
        elif not ket:
            return False, "Vector Space Associativity Fail:ket"
    
    def id_test(self):
        bra = self.id_test_bra()
        ket = self.id_test_ket()
        total = bra * ket
        if total:
            return True, "Pass"
        elif bra == ket:
            return False, "Vector Space Identity Fail:both bra and ket"
        elif not bra:
            return False, "Vector Space Identity Fail:bra"
        elif not ket:
            return False, "Vector Space Identity Fail:ket"
        
    def inv_test(self):
        bra = self.inv_test_bra()
        ket = self.inv_test_ket()
        total = bra * ket
        if total:
            return True, "Pass"
        elif bra == ket:
            return False, "Vector Space Inverse Fail:both bra and ket"
        elif not bra:
            return False, "Vector Space Inverse Fail:bra"
        elif not ket:
            return False, "Vector Space Inverse Fail:ket"
    
    def sca_mul_test(self):
        bra = self.sca_mul_test_bra()
        ket = self.sca_mul_test_ket()
        total = bra * ket
        if total:
            return True, "Pass"
        elif bra == ket:
            return False, "Vector Space Scalar Multiplication Fail:both bra and ket"
        elif not bra:
            return False, "Vector Space Scalar Multiplication Fail:bra"
        elif not ket:
            return False, "Vector Space Scalar Multiplication Fail:ket"
    
    def id_scal_mul_test(self):
        bra = self.id_sca_mul_test_bra()
        ket = self.id_sca_mul_test_ket()
        total = bra * ket
        if total:
            return True, "Pass"
        elif bra == ket:
            return False, "Vector Space Identity Scalar Multiplication Fail:both bra and ket"
        elif not bra:
            return False, "Vector Space Identiy Scalar Multiplication Fail:bra"
        elif not ket:
            return False, "Vector Space Identity Scalar Multiplication Fail:ket"
        
    def dist_sca_add_test(self):
        bra = self.dist_sca_add_test_bra()
        ket = self.dist_sca_add_test_ket()
        total = bra * ket
        if total:
            return True, "Pass"
        elif bra == ket:
            return False, "Vector Space Distributivity of Scalar Multiplication wrt Vector Addition Fail:both bra and ket"
        elif not bra:
            return False, "Vector Space Distributivity of Scalar Multiplication wrt Vector Addition Fail:bra"
        elif not ket:
            return False, "Vector Space Distributivity of Scalar Multiplication wrt Vector Addition Fail:ket"
    
    def dist_vec_add_test(self):
        bra = self.dist_vec_add_test_bra()
        ket = self.dist_vec_add_test_ket()
        total = bra * ket
        if total:
            return True, "Pass"
        elif bra == ket:
            return False, "Vector Space Distributivity of Scalar Multiplication wrt Scalar Addition Fail:both bra and ket"
        elif not bra:
            return False, "Vector Space Distributivity of Scalar Multiplication wrt Scalar Addition Fail:bra"
        elif not ket:
            return False, "Vector Space Distributivity of Scalar Multiplication wrt Scalar Addition Fail:ket"

class unittest_innerproduct(unittest_base):
    def __init__(self):
        super().__init__()
        self.ket_u = 1/(2**(0.5)) * dp.ket([2]) + 1/(8**(0.5)) * dp.ket([4])
        self.bra_u = 1/(2**(0.5)) * dp.bra([2]) + 1/(8**(0.5)) * dp.bra([4])
        self.bra_v = (1 + 2**(0.5)) * dp.bra([6]) - 4*1j * dp.bra([7])
        self.ket_v = (1 + 2**(0.5)) * dp.ket([6]) + 4*1j * dp.ket([7])
    
    def unit_test_ip(self, end_test = False):
        zero_r,zero_m = self.zero_prod()
        norm_r,norm_m = self.norm_test()
        ortho_r, ortho_m = self.ortho_test()
        sym_r,sym_m = self.sym_test()
        add_ax_r,add_ax_m = self.add_axiom_test()
        add_form_r,add_form_m = self.add_form_test()
        sca_mul_lef_r,sca_mul_lef_m = self.sca_mul_axiom_left_test()
        mult_sp_r,mult_sp_m = self.mult_subspaces_test()
        
        total_r = [ zero_r, norm_r, ortho_r, sym_r, add_ax_r, add_form_r,sca_mul_lef_r,
                   mult_sp_r]
        total_m = np.array([ zero_m, norm_m, ortho_m, sym_m, add_ax_m, add_form_m,sca_mul_lef_m,
                   mult_sp_m])
        if np.prod(total_r):
            return True, "All Inner Product Tests Passed"
        else:
            if not end_test:
                return False, total_m[np.invert(total_r)]
            else:
                output = total_m[np.invert(total_r)].tolist()
                print(*output, sep="\n")
    
    def zero_prod(self):
        ket = self.ket_u
        bra = self.bra_u
        LHS = self.empty_bra * ket
        RHS = bra * self.empty_ket
        if LHS == RHS == 0:
            return True, "Pass"
        else:
            return False, "The Inner Product with the Zero Vector is zero: Fail"
    
    def norm_test(self):
        bra = self.bra_v
        ket = self.ket_v
        prod = bra * ket
        prod_nonnegative = np.real(prod) >= 0
        prod_real = np.imag(prod) == 0
        prod_tot = prod_real * prod_nonnegative
        if prod_tot:
            return True, "Pass"
        elif prod_real == prod_nonnegative:
            return False, "Norm Neither Real or Nonnegative:Fail"
        elif not prod_real:
            return False, "Norm is Real: Fail"
        elif not prod_nonnegative:
            return False, "Norm is Negative: Fail"
    
    def ortho_test(self):
        bk_orth = dp.bra([1]) * dp.ket([1])
        bk_orth_check = bk_orth == 1
        
        bk_non_orth = dp.bra([1]) * dp.ket([2])
        bk_non_orth_check = bk_non_orth == 0
        
        total = bk_orth_check * bk_non_orth_check
        
        if total:
            return True, "Pass"
        elif not bk_orth_check: 
            return False, "Orthogonal Vectors Don't Inner Product to One: Fail"
        elif not bk_non_orth_check:
            return False, "Different Vectors Don't Inner Product to Zero: Fail"
    def sym_test(self):
        LHS = self.bra_u*self.ket_v
        RHS = self.bra_v*self.ket_u
        if LHS == RHS:
            return True, "Pass"
        else:
            return False, "Symetric axiom of the innerproduct test: Fail"
        
    def add_axiom_test(self):
        #Might be reduntant
        LHS = (4 + 5j)*dp.bra([4]) * dp.ket([4])
        RHS = (4)*dp.bra([4])*dp.ket([4]) + (5j)*dp.bra([4])*dp.ket([4])
        if LHS == RHS:
            return True, "Pass"
        else:
            return False, "Additivty axiom of the innerproduct test: Fail"
        
    def add_form_test(self):
        LHS = (3 - 4j)*dp.bra([1]) * (3 + 4j)*dp.ket([1])
        RHS = (3*3)*dp.bra([1])*dp.ket([1]) + 2* np.real((3*4j)*dp.bra([1])*dp.ket([1])) + (4j*(-4j))*dp.bra([1])*dp.ket([1])
        if LHS == RHS:
            return True, "Pass"
        else:
            return False, "Additive Formula Test: Fail"
    def sca_mul_axiom_left_test(self):   
        LHS = (5*dp.bra([2]))*dp.ket([2])
        RHS = 5*dp.bra([2])*dp.ket([2])
        if LHS == RHS:
            return True, "Pass"
        else:
            return False, "Scalar multiplacation in the Bra vector in the innerproduct test: Fail"  
            
    def mult_subspaces_test(self):
        bra1u = dp.bra([2,3j])
        ket1u = dp.ket([2,3j])
        bra1v = dp.bra([3,3j])
        ket1v = dp.ket([3,3j])
        
        bra2 = (1+2)*dp.bra([2,3])
        bra2a = dp.bra([2,3])
        bra2b = 2*dp.bra([2,3])
        ket2 = dp.ket([2,3])
        bra3 = 2*dp.bra([2,3j])
        
        com = bra1u * ket1v == bra1v * ket1u
        add = bra2 * ket2 == bra2a*ket2 + bra2b*ket2
        homo = bra3 * ket1u == 2*bra1u*ket1u
        total = com * add * homo
        if total:
            return True, "Pass"
        elif com == add == homo:
            return False, "Multiple Subspaces all: Fail"
        elif not com and com == add:
            return False, "Multiple Subspaces Commutivity and Addition: Fail"
        elif not com and com == homo:
            return False, "Multiple Subspaces Commutitity and Homogeneity: Fail"
        elif not com: 
            return False, "Multiple Subspaces Commutivity: Fail"
        elif not add and add == homo:
            return False, "Multiple Subspaces Addition and Homogeneity: Fail"
        elif not add:
            return False, "Multiple Subspaces Addition: Fail"
        elif not homo:
            return False, "Multiple Subspaces Homogeneity: Fail"

class unittest_qop(unittest_base):
    def __init__(self):
        super().__init__()
        cav = dp.fock_subspace(index = 0)
        atom = dp.two_level_subspace(index = 1)
        self.example_ground_ket = dp.ket([1,'g'])
        self.example_excited_ket = dp.ket([0,'e'])
        self.a = cav.a
        self.adag = cav.adag
        self.sig_p = atom.sigma_plus
        self.sig_m = atom.sigma_minus
    
    def unit_test_qop(self, end_test = False):
        Focka_r, Focka_m = self.Fock_a_op_test()
        Fockadag_r, Fockadag_m = self.Fock_adag_op_test()
        TwoLevelsp_r, TwoLevelsp_m = self.Two_Level_sp_op_test()
        TwoLevelsm_r, TwoLevelsm_m = self.Two_Level_sm_op_test()
        #Floqet1_r, Floqet1_m = self.Floquet_1_op_test()
        
        Conj_a_r, Conj_a_m = self.Conjugate_a_Test()
        Conj_adag_r, Conj_adag_m = self.Conjugate_adag_Test()
        Conj_sp_r, Conj_sp_m = self.Conjugate_sp_Test()
        Conj_sm_r, Conj_sm_m = self.Conjugate_sm_Test()
        #Conj_Floqet1_r, Conj_Floqet1_m = self.Floqet_conj_test()
        
        Sca_mul_r,Sca_mul_m = self.Scalar_mult_op_test()
        Sca_div_r,Sca_div_m = self.Scalar_div_op_test()
        Op_add_r, Op_add_m  = self.Op_add_test()
        Op_minus_r, Op_minus_m  = self.Op_minus_test()
        Op_Same_mul_r, Op_Same_mul_m  = self.Op_Same_Mul_test()
        Op_Diff_mul_r, Op_Diff_mul_m  = self.Op_Diff_Mul_test()
        
        total_r = [Focka_r, Fockadag_r, TwoLevelsp_r, TwoLevelsm_r, 
                   Conj_a_r, Conj_adag_r, Conj_sp_r, Conj_sm_r,
                   Sca_mul_r, Sca_div_r, Op_add_r, Op_minus_r, Op_Same_mul_r, Op_Diff_mul_r]
        total_m = np.array([Focka_m, Fockadag_m, TwoLevelsp_m, TwoLevelsm_m, 
                   Conj_a_m, Conj_adag_m, Conj_sp_m, Conj_sm_m,
                   Sca_mul_m, Sca_div_m, Op_add_m, Op_minus_m, Op_Same_mul_m, Op_Diff_mul_m])
        if np.prod(total_r):
            return True,  "All QOP Tests Passed"
        else:
            if not end_test:
                return False, total_m[np.invert(total_r)]
            else:
                output = total_m[np.invert(total_r)].tolist()
                print(*output, sep="\n")
        
    def Fock_a_op_test(self):
        op = self.a
        bra_LHS = dp.bra([3])* op
        bra_RHS = 2.0 * dp.bra([4])
        bra = bra_LHS == bra_RHS
        
        ket_LHS = op * dp.ket([4])
        ket_RHS = 2.0 * dp.ket([3])
        ket = ket_LHS == ket_RHS
        
        if bra*ket:
            return True, "Pass"
        elif bra == ket:
            return False, "Fock Space 'a' Operators Act incorrectly on both bra and ket: Fail"
        elif not ket:
            return False, "Fock Space 'a' Operators Acting on Ket: Fail"
        elif not bra:
            return False, "Fock Space 'a' Operators Acting on Bra: Fail"
        
    def Fock_adag_op_test(self):
        op = self.adag
        bra_LHS = dp.bra([4])* op
        bra_RHS = 2.0 * dp.bra([3])
        bra = bra_LHS == bra_RHS
        
        ket_LHS = op * dp.ket([3])
        ket_RHS = 2.0 * dp.ket([4])
        ket = ket_LHS == ket_RHS
        
        if bra*ket:
            return True, "Pass"
        elif bra == ket:
            return False, "Fock Space 'a dagger' Operators Act incorrectly on both bra and ket: Fail"
        elif not ket:
            return False, "Fock Space 'a dagger' Operators Acting on Ket: Fail"
        elif not bra:
            return False, "Fock Space 'a dagger' Operators Acting on Bra: Fail"
    
    def Two_Level_sp_op_test(self):
        op = self.sig_p
        bra_LHS = dp.bra([1,'e']) * op
        bra_RHS = dp.bra([1,'g'])
        bra = bra_LHS == bra_RHS
        
        ket_LHS = op * dp.ket([1,'g'])
        ket_RHS = dp.ket([1,'e'])
        ket = ket_LHS == ket_RHS
        
        if bra*ket:
            return True, "Pass"
        elif bra == ket:
            return False, "Two Level 'sigma plus' Operators Act incorrectly on both bra and ket: Fail"
        elif not ket:
            return False, "Two Level 'sigma plus' Operators Acting on Ket: Fail"
        elif not bra:
            return False, "Two Level 'sigma plus' Operators Acting on Bra: Fail"
        
    def Two_Level_sm_op_test(self):
        op = self.sig_p
        bra_LHS = dp.bra([1,'e']) * op
        bra_RHS = dp.bra([1,'g'])
        bra = bra_LHS == bra_RHS
        
        ket_LHS = op * dp.ket([1,'g'])
        ket_RHS = dp.ket([1,'e'])
        ket = ket_LHS == ket_RHS
        
        if bra*ket:
            return True, "Pass"
        elif bra == ket:
            return False, "Two Level 'sigma minus' Operators Act incorrectly on both bra and ket: Fail"
        elif not ket:
            return False, "Two Level 'sigma minus' Operators Acting on Ket: Fail"
        elif not bra:
            return False, "Two Level 'sigma minus' Operators Acting on Bra: Fail"
    
    def Floquet_1_op_test(self):
        #NOT IMPLEMENTED
        pass
    
    def Conjugate_a_Test(self):
        op = self.a
        conjop = self.adag
        LHS = op.conj() * dp.ket([2])
        RHS = conjop * dp.ket([2])
        if LHS == RHS:
            return True, "Pass"
        else:
            return False, "The Conjugate of 'a' is not 'a dagger': Fail"
    
    def Conjugate_adag_Test(self):
        op = self.adag
        conjop = self.a
        LHS = op.conj() * dp.ket([2])
        RHS = conjop * dp.ket([2])
        if LHS == RHS:
            return True, "Pass"
        else:
            return False, "The Conjugate of 'a' is not 'a dagger': Fail"
    
    def Conjugate_sp_Test(self):
        op = self.sig_p
        conjop = self.sig_m
        LHS = op.conj() * dp.ket([2,'e'])
        RHS = conjop * dp.ket([2,'e'])
        if LHS == RHS:
            return True, "Pass"
        else:
            return False, "The Conjugate of 'sigma plus' is not 'sigma minus': Fail"
    
    def Conjugate_sm_Test(self):
        op = self.sig_m
        conjop = self.sig_p
        LHS = op.conj() * dp.ket([2,'g'])
        RHS = conjop * dp.ket([2,'g'])
        if LHS == RHS:
            return True, "Pass"
        else:
            return False, "The Conjugate of 'sigma minus' is Not 'sigma plus': Fail"
        
    def Floqet_conj_test(self):
        #NOT IMPLEMENTED
        pass
    
    def Scalar_mult_op_test(self):
        ket = dp.ket([1,'g'])
        LHS_a = (4* self.a) * ket
        RHS_a = 4 * ( self.a * ket)
        a_result = LHS_a == RHS_a
        
        LHS_sp = (4 * self.sig_p) * ket
        RHS_sp = 4 * (self.sig_p * ket)
        sp_result = LHS_sp == RHS_sp
        
        if a_result*sp_result:
            return True, "Pass"
        elif a_result == sp_result:
            return False, "Scalar Multiplacation Works For Neither 'a' nor 'sigma plus': Fail"
        elif not a_result:
            return False, "Scalar Multiplication Doesn't Work for 'a': Fail"
        elif not sp_result:
            return False, "Scalar Multiplication Doesn't Work for 'sigma plus': Fail"
        
    def Scalar_div_op_test(self):
        ket = dp.ket([1,'g'])
        LHS_a = (self.a / 4) * ket
        RHS_a = ( self.a * ket) / 4
        a_result = LHS_a == RHS_a
        
        LHS_sp = (self.sig_p / 4) * ket
        RHS_sp =(self.sig_p * ket)/4 
        sp_result = LHS_sp == RHS_sp
        
        if a_result*sp_result:
            return True, "Pass"
        elif a_result == sp_result:
            return False, "Scalar Division Works For Neither 'a' Nor 'sigma plus': Fail"
        elif not a_result:
            return False, "Scalar Division Doesn't Work for 'a': Fail"
        elif not sp_result:
            return False, "Scalar Division Doesn't Work for 'sigma plus': Fail"
    
    def Op_add_test(self):
        ket = dp.ket([1,'g'])
        LHS_a = (self.a + self.a) * ket
        RHS_a = self.a * ket + self.a * ket
        a_result = LHS_a == RHS_a
        
        LHS_sp = (self.sig_p + self.sig_p) * ket
        RHS_sp = self.sig_p * ket + self.sig_p * ket
        sp_result = LHS_sp == RHS_sp
        
        if a_result*sp_result:
            return True, "Pass"
        elif a_result == sp_result:
            return False, "Operator Addition Works For Neither 'a' Nor 'sigma plus': Fail"
        elif not a_result:
            return False, "Operator Addition Doesn't Work for 'a': Fail"
        elif not sp_result:
            return False, "Operator Addition Doesn't Work for 'sigma plus': Fail"
    
    def Op_minus_test(self):
        ket = dp.ket([1,'g'])
        LHS_a = (self.a - self.a) * ket
        RHS_a = self.a * ket - self.a * ket
        a_result = LHS_a == RHS_a
        
        LHS_sp = (self.sig_p - self.sig_p) * ket
        RHS_sp = self.sig_p * ket - self.sig_p * ket
        sp_result = LHS_sp == RHS_sp
        
        if a_result*sp_result:
            return True, "Pass"
        elif a_result == sp_result:
            return False, "Operator Subtraction Works For Neither 'a' Nor 'sigma plus': Fail"
        elif not a_result:
            return False, "Operator Subtraction Doesn't Work for 'a': Fail"
        elif not sp_result:
            return False, "Operator Subtraction Doesn't Work for 'sigma plus': Fail"
    
    def Op_Diff_Mul_test(self):
        LHS = self.a * self.sig_p* dp.ket([1,'g'])
        RHS_1 = self.a * dp.ket([1,'e'])
        RHS_2 = self.sig_p * dp.ket([0,'g'])
        RHS_3 = dp.ket([0,'e'])
        total = LHS == RHS_1 == RHS_2 == RHS_3
        if total:
            return True, "Pass"
        else:
            return False, "Multiplication of Different Operators: Fail"
    
    def Op_Same_Mul_test(self):
        LHS = self.adag * self.adag* dp.ket([0,'g'])
        RHS_1 = self.adag * dp.ket([1,'g'])
        RHS_2 = 1.4142135623730951 * dp.ket([2, 'g'])
        total = LHS == RHS_1 == RHS_2
        if total:
            return True, "Pass"
        else:
            return False, "Multiplication of Same Operators: Fail"

class unittest_dynamics(unittest_base):
    def __init__(self):
        super().__init__()

        self.cav = dp.fock_subspace(index = 0)
        self.atom = dp.two_level_subspace(index = 1)
        
        kappa =123179248988002.58
        ham = self.generate_ham()
        initialbasis = [(1,'g'),(0,'e')]
        self.sys1 = dp.qsys(ham,initialstates = initialbasis, n_int = 0)
        self.sys2 = dp.qsys(ham,initialstates = initialbasis, n_int = 0, 
                            jump_ops = [np.sqrt(kappa) * self.cav.a])
        
        gs = dp.ket((0,'e'))
        rho_0_op = gs * gs.conj()
        self.rho_0_1 = self.sys1.matrix(rho_0_op)
        self.rho_0_2 = self.sys2.matrix(rho_0_op)
        
        T = T = 1e-12
        self.times = np.linspace(0, T, 500)
        
    
    def unit_test_dyn(self, end_test = False):
        evo_r,evo_m = self.unit_evo_test()
        linb_r,linb_m = self.lindblad_test()
        vn_r, vn_m  = self.VonNeumann_test()
        liou_r, liou_m  = self.liouville_test()
        sint_r, sint_m  = self.schrodint_test()
        qj_r, qj_m  = self.quantumjumps_test()
        
        total_r = [evo_r, linb_r, vn_r, liou_r, sint_r, qj_r]
        total_m = np.array([evo_m, linb_m, vn_m, liou_m, sint_m, qj_m])
        if np.prod(total_r):
            return True,  "All QOP Tests Passed"
        else:
            if not end_test:
                return False, total_m[np.invert(total_r)]
            else:
                output = total_m[np.invert(total_r)].tolist()
                print(*output, sep="\n")
                
    def generate_ham(self):
        delta_c = 2127293674469188.0
        delta_a = 2127293674469188.0
        g = 60112183589120.08
        ham = delta_c * self.cav.n +  delta_a * self.atom.sigma_z
        ham_int = g * (self.cav.a * self.atom.sigma_plus + 
                       self.cav.adag * self.atom.sigma_minus)
        return ham + ham_int
    
    def error_tol(self,f1, f2, tol):
        return abs(f1 - f2) < tol
    
    def unit_evo_test(self, tol = 1e-8):
        soln = dp.unitaryevolution(self.rho_0_1, self.times, self.sys1)
        soln.solve()
        cav_pop = np.real(soln.soln[:,0,0])[-1]
        test = self.error_tol(cav_pop, 0.16770312820379182,tol)
        if test:
            return True, "Pass"
        else:
            return False, "Unitary Evolution Test: Fail"
    
    def lindblad_test(self, tol = 1e-8):
        soln = dp.lindbladint(self.rho_0_2, self.times, self.sys2)
        soln.solve()
        cav_pop = np.real(soln.soln[:,0,0])[-1]
        test = self.error_tol(cav_pop, 1.0402581147105376e-13,tol)
        if test:
            return True, "Pass"
        else:
            return False, "Linbland Test: Fail"
    
    def VonNeumann_test(self, tol = 1e-8):
        soln = dp.vonneumannint(self.rho_0_2, self.times, self.sys2)
        cav_pop = np.real(soln.soln[:,0,0])[-1]
        test = self.error_tol(cav_pop, 0.167700032724725,tol)
        if test:
            return True, "Pass"
        else:
            return False, "VonNeumann Test: Fail"
        
    def liouville_test(self, tol = 1e-8):
        soln = dp.liouville(self.rho_0_2, self.times, self.sys2)
        soln.solve()
        cav_pop = np.real(soln.soln[:,0,0])[-1]
        #The answer is far smaller than default error
        test = self.error_tol(cav_pop, 2.315584825739848e-27,tol)
        if test:
            return True, "Pass"
        else:
            return False, "Liouville Test: Fail"
    
    def schrodint_test(self, tol=1e-8):
        rho_0_sint = np.array([0,1]) 
        soln = dp.schrodint(rho_0_sint, self.times, self.sys1)
        soln.solve()
        cav_pop = np.abs(soln.soln[:,0])[-1]**2
        test = self.error_tol(cav_pop, 0.1677069988800502,tol)
        if test:
            return True, "Pass"
        else:
            return False, "Schrodint Test: Fail"
        
    def quantumjumps_test(self, tol=1e-8):
        psi_0_qj = np.array([0,1,0])
        times = np.linspace(0,5e-14,100)
        jpsolve = dp.quantumjumps(psi_0_qj, times, self.sys2, test=True)
        # Deterministic evolution of basis states
        jpsolve.gen_bstate_evolution()
        # Monte Carlo simulation
        mean_rho = jpsolve.calc_rho(100)
        cav_pop = np.real(mean_rho[-1,0,0])
        test = self.error_tol(cav_pop, 0.017304829099494597,tol)
        if test:
            return True, "Pass"
        else:
            return False, "Quantum Jumps Test: Fail"
    
class unit_test(unittest_qop, unittest_innerproduct, unittest_vector, 
                unittest_dynamics):
    def __init__(self):
        super().__init__()
        unittest_dynamics.__init__(self)
    def unit_test(self):
        test_vec_b, test_vec_m = self.unit_test_vec()
        test_ip_b, test_ip_m = self.unit_test_ip()
        test_qop_b, test_qop_m = self.unit_test_qop()
        test_dyn_b, test_dyn_m = self.unit_test_dyn()
        tot_b = test_vec_b * test_ip_b * test_qop_b * test_dyn_b
         
        tot_m = np.array([])
        
        if not test_vec_b:
            tot_m = np.append(tot_m, test_vec_m)
        
        if not test_ip_b:
            tot_m = np.append(tot_m, test_ip_m)
        if not test_qop_b:
            tot_m = np.append(tot_m, test_qop_m)
        if not test_dyn_b:
            tot_m = np.append(tot_m, test_dyn_m)
                
        if np.prod(tot_b):
            print("All Tests Passed")
            return True, "All Tests Passed"
        else:
            output = tot_m.tolist()
            print(*output, sep="\n")
            
if __name__ == "__main__":
    unit_test().unit_test()
    