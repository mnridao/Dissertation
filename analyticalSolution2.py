"""
Student ID: 31827379
"""

import pythonCode.plotters as plotters
import pythonCode.helpers as helpers
import pythonCode.analyticalSolutions as analy
    
if __name__ == "__main__":
        
    # Setup solver.
    solver = helpers.setupSolver(xbounds=[0, 12e6], dx=10e3, 
                                 endtime=4e6, dt=10, schemeKey=1, 
                                 sPhiCoeffFlag=True, sFuncFlag=True)
    
    # Add plotter.
    solver.plotter = plotters.plotWithAnalytical2
    
    # Initial condition.
    solver.model.grid.phi = helpers.complexZerosIC(solver.model.grid.X)
    
    # Add analytical solution that will be evaluated each iteration.
    solver.addCustomEquation("analytic", analy.analyticalSolution2, 
                             args=(helpers.complexZerosIC, solver), 
                             store=False)
    
    u = 0.
    solver.run(u)