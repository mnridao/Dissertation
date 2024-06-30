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
    # Calculate courant number.
    c = u*solver.dt/solver.model.grid.dx
    c = 0.5
    dxs = u*dts/c
    
    errorsIm = np.zeros_like(dts, dtype=np.float64)
    errorsRe = np.zeros_like(dts, dtype=np.float64)
    
    intdt = dts.min()
    intnt = int(np.ceil(endtime/dt))
    
    # Iterate through the different time resolutions.
    for i, (dti, dxi) in enumerate(zip(dts, dxs)):
                        
        # Change the timestep.
        solver.setNewTimestep(dti, endtime)
        
        # Adjust dx for fixed c.
        solver.model.grid = helpers.createNewGrid(solver.model.grid.xbounds, dxi)
        
        # Reset Initial condition.
        # solver.model.grid.phi = helpers.setInitialCondition(icKey, *icArgs)
        solver.model.grid.phi = helpers.complexZerosIC(solver.model.grid.X)
        # solver.model.grid.phi = analy.analyticalSolution1(0, helpers.piecewiseWaveIC, u, solver)
        # solver.model.grid.phi = analy.advectionSolution(0, helpers.piecewiseWaveIC, u, solver)
        
        # u = 0
        solver.addCustomEquation("analytic", analy.analyticalSolution5_v2,
                                  args=(helpers.complexZerosIC, u, solver, 
                                        intdt, intnt), 
                                  nres=solver.model.grid.phi.shape[0],
                                  store=False)
        
        # solver.addCustomEquation("analytic", aFunc, *aArgs, 
        #                           nres=solver.model.grid.phi.shape[0],
        #                           store=False)
        
        print(f"c={u*dti/solver.model.grid.dx}, \tdx={solver.model.grid.dx}, \tnx={solver.model.grid.nx}, \tdt={solver.dt}, \tnt={solver.nt}, \tnt*dt={solver.nt*solver.dt}, \tnx*dx={solver.model.grid.nx*solver.model.grid.dx}")
                
        # Run the solver.
        solver.plotResults=False
        solver.plotter = plotters.plotWithAnalytical3
        solver.plotEveryNTimestep = 10
        solver.store=False
        solver.run(u)
        
        # plotters.plotWithAnalytical3(solver) # debugging.
        
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
    
    # Plot 2nd order line.
    e0 = errorsRe[0]
    e1 = e0*(dts[-1]/dts[0])**2
    plt.plot([dts[0], dts[-1]], [e0, e1], 'k--', label='2-order')
    
    plt.yscale("log")
    plt.xscale("log")
    plt.legend(fontsize=15)
    plt.grid()
    plt.xlabel("$\Delta x$", fontsize=fontsize)
    plt.ylabel("l2 error", fontsize=fontsize)
    plt.tick_params(axis='both', which='both', labelsize=15)
#%%
if __name__ == "__main__":
    
    # Setup solver.
    params = PsiParameters()
    # params.omega = 0.005
    # params.omega = 2*np.pi/100
    params.lx = 3
    endtime = 50 # 6 periods.
    # endtime = 2000
    dt = 0.5
    # solver = helpers.setupSolver(xbounds=[0, 10], dx=10e-3, 
    #                              endtime=endtime, dt=dt, schemeKey=1, 
    #                              sPhiCoeffFlag=True, sFuncFlag=True,
    #                              params=params)
    solver = helpers.setupSolver(xbounds=[-5, 10], dx=10e-3, 
                                 endtime=endtime, dt=dt, schemeKey=1, 
                                 sPhiCoeffFlag=True, sFuncFlag=True,
                                 params=params)
        
    # Reset initial condition.
    # solver.model.grid.phi = helpers.complexZerosIC(solver.model.grid.X)
            
    # Add custom equation (analytical solution).
    u = 0.005
    # u = 0.
    
    intdt = 0.01
    intnt = int(np.ceil(endtime/dt))
    solver.addCustomEquation("analytic", analy.analyticalSolution5_v2,
                              args=(helpers.complexZerosIC, u, solver, 
                                    intdt, intnt), 
                              nres=solver.model.grid.phi.shape[0],
                              store=False)
    # solver.addCustomEquation("analytic", analy.analyticalSolution4_v2,
    #                           args=(helpers.complexZerosIC, u, solver, 
    #                                 solver.dt, solver.nt), 
    #                           nres=solver.model.grid.phi.shape[0],
    #                           store=False)
    # solver.addCustomEquation("analytic", analy.analyticalSolution1, 
    #                           args=(helpers.piecewiseWaveIC, u, solver), 
    #                           nres=solver.model.grid.phi.shape[0],
    #                           store=False)
    
    alphaKey=3
    c = u*solver.model.grid.dx/solver.dt
    if alphaKey == 1:
        solver.scheme.alpha = np.maximum(0.5, 1 - 1/c)
    if alphaKey == 2:
        solver.scheme.alpha = np.maximum(0.5, 1 - 1/(c + c*params.N*solver.dt))
    if alphaKey == 3:
        solver.scheme.alpha = 0.5
    
    # Change parameters.
    dts = np.array([50, 25, 10, 5, 1, 0.5, 0.1])
    # dts = np.array([50, 25, 10, 5, 1, 0.5])
    errorsIm, errorsRe = runDifferentTimesteps(u, solver, dts, endtime, 
                                                "zeros", (solver,))
    
    #%%
    # Plot the l2 errors.
    plotl2(dts, errorsIm, errorsRe)
    
    dxs = u*dts/(u*solver.dt/solver.model.grid.dx)
    plotl2(dxs, errorsIm, errorsRe)
    
    slope1, _ = np.polyfit(np.log(dts), np.log(errorsIm), 1)
    slope2, _ = np.polyfit(np.log(dts), np.log(errorsRe), 1)
    
    slope3, _ = np.polyfit(np.log(dxs), np.log(errorsIm), 1)
    slope4, _ = np.polyfit(np.log(dxs), np.log(errorsRe), 1)
    
    # Shows spatial discretisation errors
    # need fixed courant number, run for short time (dominated by truncation errors),
    # Ndt = constant