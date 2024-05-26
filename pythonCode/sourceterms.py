"""
Student ID: 31827379
"""

import numpy as np

class PsiParameters:
    """ """
    
    def __init__(self):
        """ """
        
        self.psi0  = 10
        self.N     = 0.02
        self.omega = 1.25e-4
        self.lx    = 160e3

# Parent class?

class Psi:
    """ """
    def __init__(self):
        """ """
        self.params = PsiParameters()
        
    def __call__(self, x, t):
        """ """
        
        b = np.pi/self.params.lx
        return self.params.psi0*(b*x)*np.sin(self.params.omega*t)*np.exp(-(b*x)**2)
        

class DPsiDx:
    """ """
    
    def __init__(self):
        """ """
        self.params = PsiParameters()
    
    def __call__(self, x, t):
        """ """
        
        b = np.pi/self.params.lx
        return (1j*self.params.N*self.params.psi0*b*
                (1 - 2*(b*x)**2)*np.sin(self.params.omega*t)*np.exp(-(b*x)**2))