"""
Student ID: 31827379
"""

import pythonCode.plotters as plotters
import pythonCode.helpers as helpers
import pythonCode.analyticalSolutions as analy
   
if __name__ == "__main__":
        
    # Setup solver.
    solver = helpers.setupSolver(xbounds=[0, 4e5], dx=1e3, 
                                 endtime=4e6, dt=1, schemeKey=1, 
                                 sPhiCoeffFlag=True, sFuncFlag=True)
    
    # Initial condition.
    solver.model.grid.phi = helpers.complexZerosIC(solver.model.grid.X)
        
    # Add plotter.
    solver.plotResults=True
    solver.plotEveryNTimesteps = 5
    solver.plotter = plotters.plotWithAnalytical5
    
    # Add custom equation (analytical solution).
    u = 50
    solver.addCustomEquation("analytic", analy.analyticalSolution5_v2,
                              args=(helpers.complexZerosIC, u, solver), 
                              nres=solver.model.grid.phi.shape[0],
                              store=False)
    
    # Run the solver.
    solver.store=False
    solver.run(u)