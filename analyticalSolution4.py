"""
Student ID: 31827379
"""

import numpy as np

from pythonCode.timeschemes import ExplicitImplicitRK2
from pythonCode.grid import Grid1D 
from pythonCode.solver import Solver 
from pythonCode.sourceterms import DPsiDx, PsiParameters
import pythonCode.plotters as plotters
import pythonCode.helpers as helpers
from pythonCode.integrators import Trapezoidal

def integrand(t, X, u, params):
    
    b = np.pi/params.lx 
    a = X + u*t
    return (1j*params.N*
            params.psi0*b*(1 - 2*(a*b)**2)*np.sin(params.omega*t)*np.exp(-(a*b)**2))

if __name__ == "__main__":
    
    # Setup grid.
    x0 = 0
    xL = 12e6
    dx = 10e3
    nx = int((xL - x0)/dx)
    grid = Grid1D([x0, xL], nx)
    
    # Source terms.
    N = 0.02
    sPhiCoeff = 1j*N*0
    sDPsiDx = DPsiDx()
    
    # Set time scheme.
    scheme = ExplicitImplicitRK2(sPhiCoeff, sDPsiDx)
    u = 5
    
    # Solver parameters.
    endtime = 4e6
    dt = 10
    nt = int(np.ceil(endtime/dt))
    solver = Solver(grid, scheme, dt, nt)
    
    # Initial condition.
    solver.grid.phi = (1+1j)*np.zeros_like(grid.X)
        
    # Add plotter.
    solver.plotEveryNTimesteps = 1
    solver.plotter = plotters.plotWithAnalytical4
    
    # Add custom equation (analytical solution).
    params = PsiParameters()
    integrator = Trapezoidal(integrand, dt, nt, args=(grid.X, u, params), store=False)
    
    solver.addCustomEquation("analytic", integrator, store=False)
    
    # Run the solver.
    solver.run(u)