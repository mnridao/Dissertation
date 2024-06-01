"""
Student ID: 31827379
"""

import numpy as np
import matplotlib.pyplot as plt

from pythonCode.timeschemes import ExplicitImplicitRK2
from pythonCode.grid import Grid1D 
from pythonCode.sourceterms import PsiParameters, DPsiDx
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
    sPhiCoeff = 1j*0
    sFunc = DPsiDx()
        
    # Set time scheme.
    scheme = ExplicitImplicitRK2(sPhiCoeff, sFunc)
    u = 0.
    
    # Solver parameters.
    endtime = 4e6
    dt = 10
    nt = int(np.ceil(endtime/dt))
    solver = Solver(grid, scheme, dt, nt)
    solver.plotter = plotters.plotWithAnalytical3
    
    # Initial condition.
    solver.grid.phi = np.zeros_like(solver.grid.X, dtype=np.complex128)
    
    solver.run(u)