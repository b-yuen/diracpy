import numpy as np
import diracpy as dp


class unittests:
    def __init__(self):
        self.pump = dp.floquet_subspace(index = 0)
        self.cav = dp.fock_subspace(index = 1)
        self.atom = dp.two_level_subspace(index = 2)
        self.example_ground_ket = dp.ket([0,1,'g'])
        self.example_excited_ket = dp.ket([0,0,'e'])
    
    def state_scalar_multiplacation_test(self):
        test_op = self.example_ground_ket
        new_op = np.sqrt(2) * self.cav.a  
        
        all_true = bool1 + bool2
        if all_true == True:
            
            print("[TEST ITEM] is working for all")
            pass
        else:
            return testarray





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




