"""
Student ID: 31827379
"""

import numpy as np
import pythonCode.plotters as plotters

class Solver:
    """ 
    """
    
    def __init__(self, grid, scheme, dt, nt):
        """ """
        
        # Solver parameters.
        self.grid   = grid
        self.scheme = scheme
        self.dt     = dt
        self.nt     = nt
        
        # Storage parameters.
        self.store = False
        self.history = None
        
        # Plotter parameters.
        self.plotResults = True
        self.plotEveryNTimesteps = 1
        self.plotter = plotters.defaultPlotter
    
    def run(self, u0):
        
        # Initialise storage array.
        if self.store:
            self.history = np.zeros(shape=(self.nt+1, *self.grid.phi.shape))
            self.history[0, ...] = self.grid.phi
         
        for i in range(1, self.nt+1):
            
            # Calculate new time step value.
            phi = self.scheme(self.grid, u0, self.dt, i)
            
            # Update the old timestep value.
            self.grid.phi = phi.copy()
            
            # Store if necessary.
            if self.store:
                self.history[i, ...] = self.grid.phi
            
            # Plot the thing.
            if self.plotResults:
                if i % self.plotEveryNTimesteps == 0:
                    self.plotter(self.grid, i*self.dt, u0)