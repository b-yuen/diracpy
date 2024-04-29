import numpy as np
import diracpy as dp


class unittest_base:
    def __init__(self):
        self.pump = dp.floquet_subspace(index = 0)
        self.cav = dp.fock_subspace(index = 1)
        self.atom = dp.two_level_subspace(index = 2)
        self.example_ground_ket = dp.ket([0,1,'g'])
        self.example_excited_ket = dp.ket([0,0,'e'])
    
    def string_test(self):
        psi0 = self.example_ground_ket
        psi0_str = str(psi0)
        print(psi0_str)
    
    
class unittest_states(unittest_base):
    def string_test(self):
        psi0 = self.example_ground_ket
        psi0_str = str(psi0)
        print(psi0_str)
    
    def empty_action_ket(self):
        atom = self.atom
        psi0 = dp.ket()
        psi_m = atom.sigma_minus * psi0
        psi_p = atom.sigma_plus * psi0
        psi_mp = atom.sigma_minus * atom.sigma_plus * psi0
        psi_pm = atom.sigma_plus * atom.sigma_minus * psi0
        ##psi1 = cav.a * psi0
        #print('psi0 = ', end='')
        #psi0.print()
        #print('psi1 = ', end='')
        #psi1.print()
    
    def state_scalar_multiplacation_test(self):
        test_state = self.example_ground_ket
        new_state = np.sqrt(2) * self.cav.a  
        
        all_true = bool1 + bool2
        if all_true == True:
            
            print("[TEST ITEM] is working for all")
            pass
        else:
            return testarray



class unittest_qop(unittest_base):
    

def blank_test():
    
    bool1 = True
    bool2 = False    
    testarray = [["test_item1", bool1],["test_item2", bool2]] #This can be tidied
    
    
    all_true = bool1 + bool2
    if all_true == True:
        
        print("[TEST ITEM] is working for all")
        pass
    else:
        return testarray




