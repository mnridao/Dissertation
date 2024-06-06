"""
Student ID: 31827379
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.animation import PillowWriter
from functools import partial

import pythonCode.plotters as plotters
import pythonCode.helpers as helpers
import pythonCode.analyticalSolutions as analy
from pythonCode.integrators import Trapezoidal
from pythonCode.sourceterms import PsiParameters, DPsiDx, Psi

def animateSolution(solver, phiN, phiA, ylabelIm, ylabelRe, figname):
        
    def updateSolution(frame, yboundsIm, yboundsRe):
        """ 
        """
        t = solver.dt*frame
        
        # Update the data.
        axs[0].clear()
        axs[1].clear()
        
        # Plot the real part on the left subplot
        if type(phiA) != type(None):
            axs[0].plot(solver.model.grid.X, phiA[frame].real, 'k--')
        axs[0].plot(solver.model.grid.X, phiN[frame].real)
        axs[0].set_xlabel('X [m]', fontsize=fontsize)
        axs[0].set_ylabel(ylabelRe, fontsize=fontsize)
        axs[0].grid(True)
        axs[0].set_ylim(yboundsRe)
        
        # Plot the imaginary part on the right subplot   
        if type(phiA) != type(None):
            axs[1].plot(solver.model.grid.X, phiA[frame].imag, 'k--', label="analytic")     
        axs[1].plot(solver.model.grid.X, phiN[frame].imag, label="numerical")
        axs[1].set_xlabel('X [m]', fontsize=fontsize)
        axs[1].set_ylabel(ylabelIm, fontsize=fontsize)
        axs[1].grid(True)
        axs[1].text(0.05, textposy, f'Time: {t:.2f} s', 
                    # transform=ax1.transAxes, fontsize=20, verticalalignment='top'
                    )
        axs[1].legend(fontsize=fontsize)
        axs[1].set_ylim(yboundsIm)
        
        return axs,
             
    imMax = phiN.imag.max()
    imMin = phiN.imag.min()
    reMax = phiN.real.max()
    reMin = phiN.real.min()
    
    imBuffer = 0.1*(imMax - imMin)
    reBuffer = 0.1*(reMax - reMin)
    
    textposy = imMax
    
    yboundsIm = [imMin-imBuffer, imMax+imBuffer]
    yboundsRe = [reMin-reBuffer, reMax+reBuffer]
    
    # Set the fontsize here.
    fontsize=15
        
    # Plot the thing.
    fig, axs = plt.subplots(1, 2, figsize=(15, 7))

    anim = FuncAnimation(fig, partial(updateSolution, yboundsRe=yboundsRe, yboundsIm=yboundsIm), 
                         frames=phiN.shape[0], interval=20)
    
    plt.show()
    
    # Save the animation as a GIF
    writer = PillowWriter(fps=30)
    anim.save(figname, writer=writer)

if __name__ == "__main__":
    
    # Setup solver.
    params = PsiParameters()
    params.N = 0.02
    T = 2*np.pi/params.N # Period of oscillation
    endtime = T 
    dt = 0.5
    solver = helpers.setupSolver(xbounds=[0, 10], dx=10e-2, 
                                  endtime=endtime, dt=dt, schemeKey=1, 
                                  sPhiCoeffFlag=True, sFuncFlag=False,
                                  params=params)
    
    # Initial condition.
    u = 0.01
    # solver.model.grid.phi = helpers.setInitialCondition("a1", solver)
    solver.model.grid.phi = analy.analyticalSolution1(0, helpers.piecewiseWaveIC, u, solver)
            
    # Add analytical solution that will be evaluated each iteration.
    solver.addCustomEquation("analytic", analy.analyticalSolution1, 
                              args=(helpers.piecewiseWaveIC, u, solver), 
                              nres=solver.model.grid.phi.shape[0],
                              store=True)
    
    # Run the solver.
    solver.plotResults = False
    solver.plotter = plotters.plotWithAnalytical1
    solver.store = True
    solver.run(u)
    
    
    
    #%% Do the animation.
    
    # animateSolution(solver, solver.history, solver.getCustomData("analytic"), 
    #                 "Nw", "b", "oscillationsNoAdvection.gif")
    
    # Increase N.
    # solver.model.params.N = 0.001
    # solver.setNewTimestep(10, endtime)
    # solver.run(u)
    animateSolution(solver, solver.history, solver.getCustomData("analytic"), 
                    "Nw", "b", "oscillationsAdvection.gif")
    
    #%% Streamfunction.
            
    # Source terms.
    sPsi = Psi()
    sDPsiDx = DPsiDx()
    
    # Change the parameters.
    sPsi.params.lx = 3
    sDPsiDx.params.lx = 3
    
    sPsi.params.omega = 2*np.pi/100
    sDPsiDx.params.omega = 2*np.pi/100
        
    streamfunction = np.zeros(shape=(solver.nt, solver.model.grid.X.shape[0]))
    grad = np.zeros_like(streamfunction)
    for i in range(solver.nt):
        
        grad[i] = sDPsiDx(solver.model.grid.X, i*solver.dt).imag/sPsi.params.N
        streamfunction[i] = sPsi(solver.model.grid.X, i*solver.dt)
     
    var = streamfunction + 1j*grad
    animateSolution(solver, var, None, "dpsidx", "$\psi$", "streamfunction.gif")
    
    #%% Streamfunction in equation (no advection).
    
    # Setup solver.
    params = PsiParameters()
    params.omega = 2*np.pi/100
    params.lx = 3
    endtime = 200
    dt = 0.5
    solver = helpers.setupSolver(xbounds=[0, 10], dx=10e-3, 
                                 endtime=endtime, dt=dt, schemeKey=1, 
                                 sPhiCoeffFlag=False, sFuncFlag=True, params=params) # i think params arent change
    
    # Initial condition.
    solver.model.grid.phi = helpers.complexZerosIC(solver.model.grid.X)
    
    # Add analytical solution that will be evaluated each iteration.
    solver.addCustomEquation("analytic", analy.analyticalSolution3, 
                             args=(helpers.complexZerosIC, solver), 
                             nres=solver.model.grid.phi.shape[0],
                             store=True)
    
    # Add plotter.
    solver.plotResults = False
    
    # Run without advection.
    u = 0.
    solver.store = True
    solver.run(u)
    
    animateSolution(solver, solver.history, solver.getCustomData("analytic"), 
                    "Nw", "b", "streamfunctionNoAdvection.gif")
        
    #%%
    
    # Setup solver.
    params = PsiParameters()
    params.omega = 2*np.pi/100
    params.lx = 3
    endtime = 600 # 6 periods.
    dt = 0.5
    solver = helpers.setupSolver(xbounds=[0, 10], dx=10e-3, 
                                 endtime=endtime, dt=dt, schemeKey=1, 
                                 sPhiCoeffFlag=False, sFuncFlag=True,
                                 params=params)
    
    # Initial condition.
    solver.model.grid.phi = helpers.complexZerosIC(solver.model.grid.X)
        
    # # Add plotter.
    # solver.plotEveryNTimesteps = 1
    # solver.plotter = plotters.plotWithAnalytical4
    
    # Add custom equation (analytical solution).
    u = 0.01
    integrator = Trapezoidal(analy.integrand4, solver.dt, solver.nt, 
                             args=(solver.model.grid.X, u, solver.model.params), 
                             store=False)
    
    solver.addCustomEquation("integrator", integrator, store=False)
    solver.addCustomEquation("analytic", analy.analyticalSolution4, args=(solver,),
                             nres = solver.model.grid.phi.shape[0],
                             store=True)
    
    # Run the solver.
    solver.plotResults=False
    solver.store = True
    solver.run(u)
    
    animateSolution(solver, solver.history, solver.getCustomData("analytic"), 
                    "Nw", "b", "streamfunctionAdvection.gif")
    
    #%%
    
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
    
    # Add plotter.
    solver.plotter = plotters.plotWithAnalytical2
    
    # Initial condition.
    solver.model.grid.phi = helpers.complexZerosIC(solver.model.grid.X)
    
    # Add analytical solution that will be evaluated each iteration.
    solver.addCustomEquation("analytic", analy.analyticalSolution2, 
                             args=(helpers.complexZerosIC, solver), 
                             nres = solver.model.grid.phi.shape[0],
                             store=True)
    
    u = 0.
    
    # Run the solver.
    solver.plotResults=False
    solver.store = True
    solver.run(u)
    solver.run(u)
    
    animateSolution(solver, solver.history, solver.getCustomData("analytic"), 
                    "Nw", "b", "allNoAdvection.gif")
    
    #%% All with Advection.
    
    # Change the endtime (nt)
    endtime = 600 # 6 periods.
    nt = int(np.ceil(endtime/solver.dt))
    solver.nt = nt
    
    # Reset initial condition.
    solver.model.grid.phi = helpers.complexZerosIC(solver.model.grid.X)
            
    # Add custom equation (analytical solution).
    u = 0.01
    integrator = Trapezoidal(analy.integrand5, solver.dt, solver.nt, 
                             args=(solver.model.grid.X, u, solver.model.params), 
                             store=False)
    
    solver.addCustomEquation("integrator", integrator, store=False)
    solver.addCustomEquation("analytic", analy.analyticalSolution5, args=(solver,), 
                             nres = solver.model.grid.phi.shape[0],
                             store=True)
    
    # Run the solver.
    solver.run(u)
    
    animateSolution(solver, solver.history, solver.getCustomData("analytic"), 
                    "Nw", "b", "allAdvection.gif")
    