"""
Student ID: 31827379
"""

from abc import ABC, abstractmethod

import numpy as np

class IntegrationStepper(ABC):
    """ 
    """
    
    def __init__(self, f, h, nt, args, store=True):
        """ 
        """
        # Integrand calculation.
        self.f = f 
        self.args = args
        
        # Things for the integration step.
        self.h  = h
        self.nt = nt
        self.y0 = self.f(0, *self.args)
        
        # Accumulated integral.
        self.I1 = np.zeros_like(self.y0)
        self.I = np.zeros(shape=(self.nt, *self.y0.shape), 
                          dtype=self.y0.dtype) if store else None
        self.store = store
            
    def __call__(self, i):
        
        # Just in case ig.
        if i == 0:
            return
        
        # Calculate the area of the subinterval using trapezoidal method.
        self.calculateSubInterval(self.h*i)
        
        # Store the value at current grid point.
        if self.store:
            self.I[i-1, ...] = self.I1
    
    @abstractmethod
    def calculateSubInterval(self, t):
        pass

class Trapezoidal(IntegrationStepper):
    """ 
    """
    
    def __init__(self, f, h, nt, args=(), store=True):
        
        # Initialise parent class.
        super().__init__(f, h, nt, args, store)
                        
    def calculateSubInterval(self, t):
        """ 
        """
        
        # Approximate subinterval as trapezoid and calculate area.
        y1 = self.f(t, *self.args)
        self.I1 += 0.5*self.h*(self.y0 + y1)
        
        # Update previous subinterval function value.
        self.y0 = y1
        
class Simpsons(IntegrationStepper):
    """ 
    """
    
    def __init__(self, f, h, nt, args, store=True):
        
        # Initialise parent class.
        super().__init__(f, h, nt, args, store)
    
    def calculateSubInterval(self, t):
        pass