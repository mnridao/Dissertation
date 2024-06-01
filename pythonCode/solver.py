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
        
        # Functions that should be called each iteration stored here.
        self.customEquations = {}
        
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
            
            # Evaluate any functions added by user (e.g. energy)
            for eqn in self.customEquations.values():
                if eqn["data"] is not None:
                    if eqn["data"].ndim > 1:
                        eqn["data"][i+1] = eqn["func"](self, i)
                    else:
                        eqn["data"] = eqn["func"](self, i)
                else:
                    eqn["func"](self, i) # Don't store.
                            
            # Store if necessary.
            if self.store:
                self.history[i, ...] = self.grid.phi
            
            # Plot the thing.
            if self.plotResults:
                if i % self.plotEveryNTimesteps == 0:
                    self.plotter(self)

    def addCustomEquation(self, key, customEqn, nres=1, store=True):
        """ 
        Add an equation that will be evaluated at each timestep.
        
        Inputs
        -------
        key       : string
                    Key for the equation being added. This is for easy
                    accessing later.
        customEqn : callable object that takes a Model object as its argument.
                    Custom equation, e.g. to calculate energy, that will be
                    evaluated every timestep.
        nres      : int
                    Number of results returned from customEqn. Default is 1.
        """        
        
        # Initialise results data for custom function.
        if store:
            data = np.zeros(shape=(self.nt+1, nres))
            data[0] = customEqn(self, 0)
        else:
            data = customEqn(self, 0) # Save last value only.
            
        # Store in dict for easy accessing.
        self.customEquations[key] = {"func": customEqn, "data": data}
        
    def getCustomData(self, key):
        """ 
        Getter for the data obtained from evaluating the custom equation 
        specified by key at every timestep.
        
        Inputs
        ------
        key : string
              Key for obtaining the data from the dictionary storing the 
              custom equation results.
        """
        return self.customEquations[key]["data"]
    
    def getCustomEquation(self, key):
        """ 
        """
        
        return self.customEquations[key]["func"]