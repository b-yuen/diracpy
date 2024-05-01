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
    
    def unit_test(self, end_test = False):
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
            return False, "Commutivity Fail:both"
        elif not bra:
            return False, "Commutivity Fail:bra"
        elif not ket:
            return False, "Commutivity Fail:ket"
    
    def asoc_test(self):
        bra = self.assoc_test_bra()
        ket = self.assoc_test_ket()
        total = bra * ket
        if total:
            return True, "Pass"
        elif bra == ket:
            return False, "Associativity Fail:both"
        elif not bra:
            return False, "Associativity Fail:bra"
        elif not ket:
            return False, "Associativity Fail:ket"
    
    def id_test(self):
        bra = self.id_test_bra()
        ket = self.id_test_ket()
        total = bra * ket
        if total:
            return True, "Pass"
        elif bra == ket:
            return False, "Identity Fail:both"
        elif not bra:
            return False, "Identity Fail:bra"
        elif not ket:
            return False, "Identity Fail:ket"
        
    def inv_test(self):
        bra = self.inv_test_bra()
        ket = self.inv_test_ket()
        total = bra * ket
        if total:
            return True, "Pass"
        elif bra == ket:
            return False, "Inverse Fail:both"
        elif not bra:
            return False, "Inverse Fail:bra"
        elif not ket:
            return False, "Inverse Fail:ket"
    
    def sca_mul_test(self):
        bra = self.sca_mul_test_bra()
        ket = self.sca_mul_test_ket()
        total = bra * ket
        if total:
            return True, "Pass"
        elif bra == ket:
            return False, "Scalar Multiplication Fail:both"
        elif not bra:
            return False, "Scalar Multiplication Fail:bra"
        elif not ket:
            return False, "Scalar Multiplication Fail:ket"
    
    def id_scal_mul_test(self):
        bra = self.id_sca_mul_test_bra()
        ket = self.id_sca_mul_test_ket()
        total = bra * ket
        if total:
            return True, "Pass"
        elif bra == ket:
            return False, "Identity Scalar Multiplication Fail:both"
        elif not bra:
            return False, "Identiy Scalar Multiplication Fail:bra"
        elif not ket:
            return False, "Identity Scalar Multiplication Fail:ket"
        
    def dist_sca_add_test(self):
        bra = self.dist_sca_add_test_bra()
        ket = self.dist_sca_add_test_ket()
        total = bra * ket
        if total:
            return True, "Pass"
        elif bra == ket:
            return False, "Distributivity of Scalar Multiplication wrt Vector Addition Fail:both"
        elif not bra:
            return False, "Distributivity of Scalar Multiplication wrt Vector Addition Fail:bra"
        elif not ket:
            return False, "Distributivity of Scalar Multiplication wrt Vector Addition Fail:ket"
    
    def dist_vec_add_test(self):
        bra = self.dist_vec_add_test_bra()
        ket = self.dist_vec_add_test_ket()
        total = bra * ket
        if total:
            return True, "Pass"
        elif bra == ket:
            return False, "Distributivity of Scalar Multiplication wrt Scalar Addition Fail:both"
        elif not bra:
            return False, "Distributivity of Scalar Multiplication wrt Scalar Addition Fail:bra"
        elif not ket:
            return False, "Distributivity of Scalar Multiplication wrt Scalar Addition Fail:ket"

class unittest_innerproduct(unittest_base):
    def __init__(self):
        unittest_base.__init__(self)
        self.ket_u = 1/(2**(0.5)) * dp.ket([2]) + 1/(8**(0.5)) * dp.ket([4])
        self.bra_u = 1/(2**(0.5)) * dp.bra([2]) + 1/(8**(0.5)) * dp.bra([4])
        self.bra_v = (1 + 2**(0.5)) * dp.bra([6]) - 4*1j * dp.bra([7])
        self.ket_v = (1 + 2**(0.5)) * dp.ket([6]) + 4*1j * dp.ket([7])
    
    def unit_test(self, end_test = False):
        zero_r,zero_m = self.zero_prod()
        norm_r,norm_m = self.norm_test()
        ortho_r, ortho_m = self.ortho_test()
        sym_r,sym_m = self.sym_test()
        add_ax_r,add_ax_m = self.add_axiom_test()
        add_form_r,add_form_m = self.add_form_test()
        sca_mul_lef_r,sca_mul_lef_m = self.sca_mul_axiom_left_test()
        sca_mul_rig_r,sca_mul_rig_m = self.sca_mul_axiom_right_test()
        mult_sp_r,mult_sp_m = self.mult_subspaces_test()
        
        total_r = [ zero_r, norm_r, ortho_r, sym_r, add_ax_r, add_form_r,sca_mul_lef_r,
                   sca_mul_rig_r, mult_sp_r]
        total_m = np.array([ zero_m, norm_m, ortho_m, sym_m, add_ax_m, add_form_m,sca_mul_lef_m,
                   sca_mul_rig_m, mult_sp_m])
        if np.prod(total_r):
            return True, "All Vector Space Tests Passed"
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
        cav = dp.fock_subspace(index = 0)
        atom = dp.two_level_subspace(index = 1)
        self.example_ground_ket = dp.ket([1,'g'])
        self.example_excited_ket = dp.ket([0,'e'])
        self.a = cav.a
        self.adag = cav.adag
        self.sig_p = atom.sigma_plus
        self.sig_m = atom.sigma_minus
    
    def unit_test(self, end_test = False):
        Focka_r, Focka_m = self.Fock_a_op_test()
        Fockadag_r, Fockadag_m = self.Fock_adag_op_test()
        TwoLevelsp_r, TwoLevelsp_m = self.Two_Level_sp_op_test()
        TwoLevelsm_r, TwoLevelsm_m = self.Two_Level_sm_op_test()
        #Floqet1_r, Floqet1_m = self.Floquet_1_op_test()
        
        
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
        elif not ket:
            return False, "Two Level 'sigma plus' Operators Acting on Ket: Fail"
        elif not bra:
            return False, "Two Level 'sigma plus' Operators Acting on Bra: Fail"
    
    def Floquet_1_op_test(self):
        #NOT IMPLEMENTED
        pass
    