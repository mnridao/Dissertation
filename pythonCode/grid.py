"""
Student ID: 31827379
"""

import numpy as np

class Grid1D:
    """ """
    
    def __init__(self, xbounds, nx):
        """ """
        
        self.nx = nx
        self.xbounds = xbounds
        
        # Calculate grid spacing from inputs.
        self.dx = (self.xbounds[1] - self.xbounds[0]) / self.nx
        
        # Setup array representing the grid.
        self.X = np.linspace(self.xbounds[0], self.xbounds[1], self.nx+1)
        
        # Initialise the default state variable.
        self.resetFields()
        
    def resetFields(self):
        """ 
        """
        self.phi = np.zeros(shape=self.nx+1)