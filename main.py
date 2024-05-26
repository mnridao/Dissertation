"""
Student ID: 31827379
"""

import numpy as np

from pythonCode.timeschemes import ExplicitImplicitRK2
from pythonCode.grid import Grid1D 
from pythonCode.solver import Solver
from pythonCode.sourceterms import DPsiDx
import pythonCode.plotters as plotters

if __name__ == "__main__":
    
    # Setup grid.
    x0 = 0
    xL = 12e6
    dx = 10e3
    nx = int((xL - x0)/dx)
    grid = Grid1D([x0, xL], nx)
    
    # Source terms.
    sFunc = DPsiDx()
    sPhiCoeff = 1j*sFunc.params.N
    
    # Set time scheme.
    scheme = ExplicitImplicitRK2(sPhiCoeff, sFunc)
    
    # Solver parameters.
    endtime = 4e6
    dt = 10
    nt = int(np.ceil(endtime/dt))
    solver = Solver(grid, scheme, dt, nt)
    
    # Initial condition.
    solver.grid.phi = np.zeros_like(solver.grid.X, dtype=np.complex128)
    # u = 5.
    u = 0.
    
    # Add plotter.
    solver.plotter = plotters.plotter1
    
    # Run the solver.
    solver.run(u)