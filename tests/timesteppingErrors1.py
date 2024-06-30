"""
Student ID: 31827379
"""
import numpy as np
import matplotlib.pyplot as plt

import pythonCode.plotters as plotters
import pythonCode.helpers as helpers
import pythonCode.analyticalSolutions as analy
from pythonCode.sourceterms import PsiParameters

def runDifferentTimesteps(u, solver, dts, endtime, icKey, icArgs):
    """ 
    """
    
    errorsIm = np.zeros_like(dts)
    errorsRe = np.zeros_like(dts)

    # Iterate through the different time resolutions.
    for i, dti in enumerate(dts):
                
        # Reset Initial condition.
        solver.model.grid.phi = helpers.setInitialCondition(icKey, *icArgs)
            
        # Change the timestep.
        solver.setNewTimestep(dti, endtime)
        
        # Run the solver.
        solver.run(u)
        
        ## Calculate error ##
        phiA = solver.getCustomData("analytic")
        phiN = solver.model.grid.phi
        
        errorsIm[i] = helpers.l2(phiN.imag, phiA.imag, solver.dt)
        errorsRe[i] = helpers.l2(phiN.real, phiA.real, solver.dt)
    
    return errorsIm, errorsRe    

def runDifferentTimestepsA1(u, solver, dts, endtime):
    icKey = "a1"
    icArgs = (solver, )
    return runDifferentTimesteps(u, solver, dts, endtime, icKey, icArgs)

def plotl2(dts, errorsIm, errorsRe):
    
    plt.figure(figsize=(5, 5))
    plt.plot(dts, errorsIm, label="w")
    plt.plot(dts, errorsRe, label="b")
    plt.yscale("log")
    plt.xscale("log")
    plt.legend(fontsize=8)
    plt.grid()
    plt.show()
    plt.close()

if __name__ == "__main__":
    
    # Setup solver.
    params = PsiParameters()
    params.omega = 2*np.pi/100
    params.lx = 3
    endtime = 400 # 6 periods.
    dt = 0.5
    solver = helpers.setupSolver(xbounds=[0, 10], dx=10e-3, 
                                 endtime=endtime, dt=dt, schemeKey=1, 
                                 sPhiCoeffFlag=True, sFuncFlag=True,
                                 params=params)
        
    # Reset initial condition.
    solver.model.grid.phi = helpers.complexZerosIC(solver.model.grid.X)
            
    # Add custom equation (analytical solution).
    u = 0.01    
    solver.addCustomEquation("analytic", analy.analyticalSolution5_v2,
                              args=(helpers.complexZerosIC, u, solver), 
                              nres=solver.model.grid.phi.shape[0],
                              store=True)
        
    dts = np.array([50, 10, 5, 1, 0.5, 0.1])
    errorsIm, errorsRe = runDifferentTimestepsA1(u, solver, dts, endtime)
    
    # Plot the l2 errors.
    plotl2(dts, errorsIm, errorsRe)
    
    slope1, _ = np.polyfit(np.log(dts), np.log(errorsIm), 1)
    slope2, _ = np.polyfit(np.log(dts), np.log(errorsRe), 1)
    
    #%% Error over time.
        
    # Setup solver.
    params = PsiParameters()
    params.N = 0.005
    T = 2*np.pi/params.N # Period of oscillation
    endtime = 10*T 
    dt = 0.5
    solver = helpers.setupSolver(xbounds=[0, 12], dx=10e-2, 
                                 endtime=endtime, dt=dt, schemeKey=1, 
                                 sPhiCoeffFlag=True, sFuncFlag=False,
                                 params=params)
    
    # Initial condition.
    u = 0.
    solver.model.grid.phi = analy.analyticalSolution1(0, helpers.middleWaveIC, 
                                                      u, solver)
        
    # Add analytical solution that will be evaluated each iteration.
    solver.addCustomEquation("analytic", analy.analyticalSolution1, 
                             args=(helpers.middleWaveIC, u, solver), 
                             nres=solver.model.grid.phi.shape[0],
                             store=False)
    solver.addCustomEquation("l2Re", helpers.l2ReCustomEqn, args=(solver, ), store=True)
    solver.addCustomEquation("l2Im", helpers.l2ImCustomEqn, args=(solver, ), store=True)
    
    # Add plotter.
    solver.plotResults = False
    # solver.plotEveryNTimesteps = 10
    # solver.plotter = plotters.plotWithAnalytical1

    # Run the solver.
    solver.store = False
    solver.run(u)
    
    #%%
    
    l2Re = solver.getCustomData("l2Re")
    l2Im = solver.getCustomData("l2Im")
    
    time = np.arange(0, solver.dt*(solver.nt+1), solver.dt)
    plt.figure(figsize=(10, 10))
    plt.plot(time, l2Re, label="b")
    plt.plot(time, l2Im, label="w")
    plt.legend(fontsize=15)
    # plt.yscale("log")
    plt.grid()
    plt.show()
    plt.close()