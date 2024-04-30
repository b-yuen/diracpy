import numpy as np
import diracpy as dp


class unittest_base:
    def __init__(self):
        self.cav = dp.fock_subspace(index = 0)
        self.atom = dp.two_level_subspace(index = 1)
        self.example_ground_ket = dp.ket([1,'g'])
        self.example_excited_ket = dp.ket([0,'e'])
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
                return False, total_message[total_results]
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
    pass

class unittest_qop(unittest_base):
    pass