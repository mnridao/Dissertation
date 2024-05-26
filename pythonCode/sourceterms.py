"""
Student ID: 31827379
"""

import numpy as np

class DPsiDx:
    """ """
    
    def __init__(self):
        """ """
        
        self.psi0  = 10
        self.N     = 0.02
        self.omega = 1.25e-4
        self.lx    = 160e3
    
    def __call__(self, x, t):
        """ """
        
        b = np.pi/self.lx
        return (1j*self.N*self.psi0*b*
                (1 - 2*(b*x)**2)*np.sin(self.omega*t)*np.exp(-(b*x)**2))