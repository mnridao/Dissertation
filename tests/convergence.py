"""
Student ID: 31827379
"""
import numpy as np
import matplotlib.pyplot as plt

import pythonCode.plotters as plotters
import pythonCode.helpers as helpers
import pythonCode.analyticalSolutions as analy
from pythonCode.sourceterms import PsiParameters

from pythonCode.grid import Grid1D

def runDifferentTimesteps(u, solver, dts, endtime, icKey, icArgs):
    """ 
    """
    # Calculate courant number.
    c = u*solver.dt/solver.model.grid.dx
    
    errorsIm = np.zeros_like(dts)
    errorsRe = np.zeros_like(dts)
    
    intdt = dts.min()
    intnt = int(np.ceil(endtime/dt))
    
    # Iterate through the different time resolutions.
    for i, dti in enumerate(dts):
                    
        # Change the timestep.
        solver.setNewTimestep(dti, endtime)
        
        # Adjust dx for fixed c.
        dxi = u*solver.dt/c
        nxi = int((solver.model.grid.xbounds[1] - solver.model.grid.xbounds[0])/dxi)
        solver.model.grid = Grid1D(solver.model.grid.xbounds, nxi)
        # solver.model.grid.setNewGridSpacing(dxi)
        
        # Reset Initial condition.
        solver.model.grid.phi = helpers.setInitialCondition(icKey, *icArgs)
        
        print(f"c={u*dti/solver.model.grid.dx}, \tdx={solver.model.grid.dx}, \tdt={solver.dt}")
        
        # Add analytic solution.
        solver.addCustomEquation("analytic", analy.analyticalSolution5_v2,
                                  args=(helpers.complexZerosIC, u, solver, intdt, intnt), 
                                  nres=solver.model.grid.phi.shape[0],
                                  store=False)
        
        # Run the solver.
        solver.plotResults=False
        solver.plotter = plotters.plotWithAnalytical3
        solver.plotEveryNTimestep = 10
        solver.store=False
        solver.run(u)
        
        ## Calculate error ##
        phiA = solver.getCustomData("analytic")
        phiN = solver.model.grid.phi
        
        errorsIm[i] = helpers.l2(phiN.imag, phiA.imag, solver.dt)
        errorsRe[i] = helpers.l2(phiN.real, phiA.real, solver.dt)
    
    return errorsIm, errorsRe    

def plotl2(dts, errorsIm, errorsRe):
    
    fontsize=20
    
    plt.figure(figsize=(10, 10))
    plt.plot(dts, errorsIm, '-o', label="w")
    plt.plot(dts, errorsRe, '-o', label="b")
    plt.yscale("log")
    plt.xscale("log")
    plt.legend(fontsize=15)
    plt.grid()
    plt.xlabel("$\Delta t$ [s]", fontsize=fontsize)
    plt.ylabel("l2 error", fontsize=fontsize)
    plt.tick_params(axis='both', which='both', labelsize=15)

if __name__ == "__main__":
    
    # Setup solver.
    params = PsiParameters()
    params.omega = 2*np.pi/100
    params.lx = 3
    endtime = 200 # 6 periods.
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
                              args=(helpers.complexZerosIC, u, solver, 
                                    solver.dt, solver.nt), 
                              nres=solver.model.grid.phi.shape[0],
                              store=False)
    
    alphaKey=3
    c = u*solver.model.grid.dx/solver.dt
    if alphaKey == 1:
        solver.scheme.alpha = np.maximum(0.5, 1 - 1/c)
    if alphaKey == 2:
        solver.scheme.alpha = np.maximum(0.5, 1 - 1/(c + c*params.N*solver.dt))
    if alphaKey == 3:
        solver.scheme.alpha = 0.5
    
    # Change parameters.
    dts = np.array([25, 10, 5, 1, 0.5, 0.1, 0.05])
    errorsIm, errorsRe = runDifferentTimesteps(u, solver, dts, endtime, 
                                               "zeros", (solver,))
    #%%
    # Plot the l2 errors.
    plotl2(dts, errorsIm, errorsRe)
    
    slope1, _ = np.polyfit(np.log(dts), np.log(errorsIm), 1)
    slope2, _ = np.polyfit(np.log(dts), np.log(errorsRe), 1)
    
    # Shows spatial discretisation errors
    # need fixed courant number, run for short time (dominated by truncation errors),
    # Ndt = constant