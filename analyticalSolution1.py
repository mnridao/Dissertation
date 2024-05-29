"""
Student ID: 31827379
"""

import numpy as np

from pythonCode.timeschemes import ExplicitImplicitRK2
from pythonCode.grid import Grid1D 
from pythonCode.solver import Solver 
import pythonCode.plotters as plotters
import pythonCode.helpers as helpers
import pythonCode.analyticalSolutions as analy

if __name__ == "__main__":
    
    # Setup grid.
    x0 = 0
    xL = 12e6
    dx = 10e3
    nx = int((xL - x0)/dx)
    grid = Grid1D([x0, xL], nx)
    
    # Source terms.
    N = 0.02
    sPhiCoeff = 1j*N
    
    # Set time scheme.
    scheme = ExplicitImplicitRK2(sPhiCoeff, )
    u = 50.
    
    # Solver parameters.
    endtime = 4e6
    dt = 10
    nt = int(np.ceil(endtime/dt))
    solver = Solver(grid, scheme, dt, nt)
    
    # Initial condition.
    solver.grid.phi = analy.analyticalSolution1(helpers.middleWaveIC, sPhiCoeff, 
                                                u, grid.X, 0)
        
    # Add plotter.
    solver.plotEveryNTimesteps = 1
    solver.plotter = plotters.plotWithAnalytical1

    # Run the solver.
    solver.run(u)