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
                                  sPhiCoeffFlag=True, sFuncFlag=False)
    

    
    # Initial condition.
    u = 5.
    solver.model.grid.phi = analy.analyticalSolution1(0, helpers.middleWaveIC, 
                                                      u, solver)
    
    # Add analytical solution that will be evaluated each iteration.
    solver.addCustomEquation("analytic", analy.analyticalSolution1, 
                              args=(helpers.middleWaveIC, u, solver), 
                              nres=solver.model.grid.phi.shape[0],
                              store=True)
    
    # Add plotter.
    solver.plotResults=True
    solver.plotEveryNTimesteps = 1
    solver.plotter = plotters.plotWithAnalytical1

    # Run the solver.
    solver.store = False
    solver.run(u)