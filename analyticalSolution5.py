"""
Student ID: 31827379
"""

import pythonCode.plotters as plotters
import pythonCode.helpers as helpers
from pythonCode.integrators import Trapezoidal

from pythonCode.analyticalSolutions import integrand5, analyticalSolution5
    
if __name__ == "__main__":
        
    # Setup solver.
    solver = helpers.setupSolver(xbounds=[0, 12e6], dx=10e3, 
                                 endtime=4e6, dt=10, schemeKey=1, 
                                 sPhiCoeffFlag=True, sFuncFlag=True)
    
    # Initial condition.
    solver.model.grid.phi = helpers.complexZerosIC(solver.model.grid.X)
        
    # Add plotter.
    solver.plotEveryNTimesteps = 1
    solver.plotter = plotters.plotWithAnalytical5
    
    # Add custom equation (analytical solution).
    u = 5
    integrator = Trapezoidal(integrand5, solver.dt, solver.nt, 
                             args=(solver.model.grid.X, u, solver.model.params), 
                             store=False)
    
    solver.addCustomEquation("integrator", integrator, store=False)
    solver.addCustomEquation("analytic", analyticalSolution5, args=(solver,), store=False)
    
    # Run the solver.
    solver.run(u)