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
    solver = helpers.setupSolver(xbounds=[0, 30], dx=10e-2, 
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
                              store=False)
    
    # Run the solver.
    solver.plotResults = True
    solver.plotter = plotters.plotWithAnalytical1
    solver.plotEveryNTimesteps = 10
    solver.store = False
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
    
        plt.plot(solver.model.grid.X, streamfunction[i])    
        plt.ylim([-10, 10])
        plt.show(), plt.close()
        
    var = streamfunction + 1j*grad
    # animateSolution(solver, var, None, "dpsidx", "$\psi$", "streamfunction.gif")
    
    #%% Streamfunction in equation (no advection).
    
    import pythonCode.analyticalSolutions as analy
    
    # Setup solver.
    params = PsiParameters()
    params.omega = 2*np.pi/100
    params.lx = 3
    endtime = 200
    dt = 0.5
    solver = helpers.setupSolver(xbounds=[0, 10], dx=10e-3, 
                                 endtime=endtime, dt=dt, schemeKey=1, 
                                 sPhiCoeffFlag=False, sFuncFlag=True, params=params)
    
    # Initial condition.
    solver.model.grid.phi = helpers.complexZerosIC(solver.model.grid.X)
    
    # # Add analytical solution that will be evaluated each iteration.
    # solver.addCustomEquation("analytic", analy.analyticalSolution3, 
    #                           args=(helpers.complexZerosIC, solver), 
    #                           nres=solver.model.grid.phi.shape[0],
    #                           store=True)
    
    # integrator = Trapezoidal(analy.integrand4, solver.dt, solver.nt, 
    #                           args=(solver.model.grid.X[1:], u, solver.model.params), 
    #                           store=False)
    # solver.addCustomEquation("integrator", integrator, store=False)
    # solver.addCustomEquation("analytic", analy.analyticalSolution4, args=(solver,),
    #                           nres = solver.model.grid.phi.shape[0],
    #                           store=False)
    
    integrator = Trapezoidal(analy.integrand4, solver.dt, solver.nt, 
                              args=(0, solver.model.grid.X, u, solver.model.params), 
                              store=False)
    solver.addCustomEquation("integrator", integrator, store=False)
    solver.addCustomEquation("analytic", analy.analyticalSolution4, 
                             args=(helpers.complexZerosIC, u, solver),
                             nres = solver.model.grid.phi.shape[0],
                             store=False)
    
    # Add plotter.
    solver.plotResults = True
    solver.plotter = plotters.plotWithAnalytical3
    
    # Run without advection.
    u = 0.
    solver.nt=1
    solver.store = False
    solver.run(u)
    
    animateSolution(solver, solver.history, solver.getCustomData("analytic"), 
                    "Nw", "b", "streamfunctionNoAdvection.gif")
    
    #%%
    u = 0.01
    X = solver.model.grid.X
    for t in range(solver.nt):
        
        plt.figure(figsize=(10, 10))
        # x0 = solver.model.grid.X - u*t*solver.dt
        plt.plot(X+u*t*solver.dt, analy.integrand4(t*solver.dt, X, u, solver.model.params).imag)
        plt.ylim([-0.25, 0.25])
        plt.show()
        plt.close()
        
    #%%
    import pythonCode.analyticalSolutions as analy
    
    # Setup solver.
    params = PsiParameters()
    params.omega = 2*np.pi/100
    params.lx = 4
    endtime = 300 # 1 period = 100s.
    dt = 0.5
    solver = helpers.setupSolver(xbounds=[0, 10], dx=10e-3, 
                                 endtime=endtime, dt=dt, schemeKey=1, 
                                 sPhiCoeffFlag=False, sFuncFlag=True,
                                 params=params)
    
    # Initial condition.
    solver.model.grid.phi = helpers.complexZerosIC(solver.model.grid.X)
            
    # Add custom equation (analytical solution).
    u = 0.01
    print(f"c = {u*solver.dt/solver.model.grid.dx}")
    
    solver.addCustomEquation("analytic", analy.analyticalSolution4_v2,
                              args=(helpers.complexZerosIC, u, solver), 
                              nres=solver.model.grid.phi.shape[0],
                              store=True)
        
    # Run the solver.
    solver.plotResults=False
    solver.plotter=plotters.plotWithAnalytical3
    solver.plotEveryNTimesteps = 10
    solver.store = True
    solver.run(u)
    
    animateSolution(solver, solver.history, solver.getCustomData("analytic"), 
                    "Nw", "b", "streamfunctionAdvection22.gif")
    
    #%%
    
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
    solver.plotResults = False
    solver.plotter = plotters.plotWithAnalytical2
    solver.store = True
    solver.run(u)
    
    animateSolution(solver, solver.history, solver.getCustomData("analytic"), 
                    "Nw", "b", "allNoAdvection.gif")
    
    #%% All with Advection.
    
    import pythonCode.analyticalSolutions as analy
    
    # Setup solver.
    params = PsiParameters()
    params.omega = 2*np.pi/100
    # params.omega = 0.005
    params.lx = 3
    endtime = 400 # 6 periods.
    dt = 0.5
    solver = helpers.setupSolver(xbounds=[-5, 10], dx=10e-3, 
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
                              store=False)
    
    # Run the solver.
    solver.plotResults=True
    solver.plotter=plotters.plotWithAnalytical3
    solver.plotEveryNTimesteps = 5
    solver.store = False
    solver.run(u)
    # animateSolution(solver, solver.history, solver.getCustomData("analytic"), 
    #                 "Nw", "b", "allAdvection.gif")
    