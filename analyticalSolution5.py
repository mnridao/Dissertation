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
            params.psi0*b*(1 - 2*(a*b)**2)*np.sin(params.omega*t)*
            np.exp(-(a*b)**2)*np.exp(-1j*params.N*t))

def analyticEqn(solver, i):
    
    # Wrapper for integrator result.
    integrator = solver.getCustomEquation("integrator")
    params = PsiParameters()
    
    return integrator.I1*np.exp(1j*params.N*i*solver.dt)
    
if __name__ == "__main__":
    
    # Setup grid.
    x0 = 0
    xL = 12e6
    dx = 10e3
    nx = int((xL - x0)/dx)
    grid = Grid1D([x0, xL], nx)
    
    # Source terms.
    N = 0.02
    sPhiCoeff = 1j*N
    sDPsiDx = DPsiDx()
    
    # Set time scheme.
    scheme = ExplicitImplicitRK2(sPhiCoeff, sDPsiDx)
    u = 5
    
    # Solver parameters.
    endtime = 4e6
    dt = 5
    nt = int(np.ceil(endtime/dt))
    solver = Solver(grid, scheme, dt, nt)
    
    # Initial condition.
    solver.grid.phi = (1+1j)*np.zeros_like(grid.X)
        
    # Add plotter.
    solver.plotEveryNTimesteps = 5
    solver.plotter = plotters.plotWithAnalytical5
    # solver.plotter = plotters.plotter1
    
    # Add custom equation (analytical solution).
    params = PsiParameters()
    integrator = Trapezoidal(integrand, dt, nt, args=(grid.X, u, params), store=False)
    
    solver.addCustomEquation("integrator", integrator, store=False)
    solver.addCustomEquation("analytic", analyticEqn, store=False)
    
    # Run the solver.
    solver.run(u)
    
    #%%
    import matplotlib.pyplot as plt
    
    integrator = Trapezoidal(integrand, dt, nt, args=(grid.X, u, params), store=False)

    for i in range(1, nt+1):
        
        # Update the integral.
        integrator(None, i)
        
        phiAnaly = integrator.I1*np.exp(1j*params.N*i*dt)
        
        # Plot?
        # plt.figure(figsize=(10, 10))
        # plt.plot(grid.X, integrator.I1.imag)
        # plt.ylim([-1e-6, 1e-6])
        # plt.grid()
        # plt.show()
        # plt.close()
        
        fontsize=15
        
        # Create a figure with 1 row and 2 columns
        fig, axs = plt.subplots(1, 2, figsize=(15, 7))
        
        # Plot the real part on the left subplot
        axs[0].plot(grid.X, phiAnaly.real)
        axs[0].set_xlabel('X', fontsize=fontsize)
        axs[0].set_ylabel('b', fontsize=fontsize)
        axs[0].grid(True)
        # axs[0].set_ylim([-1e-4, 1e-4])
        axs[0].set_ylim([-1e-6, 1e-6])
        # axs[0].set_xlim(grid.xbounds)
        
        # Plot the imaginary part on the right subplot        
        axs[1].plot(grid.X, phiAnaly.imag)
        axs[1].set_xlabel('X', fontsize=fontsize)
        axs[1].set_ylabel('Nw', fontsize=fontsize)
        axs[1].grid(True)
        axs[1].set_ylim([-0.3e-5, 0.3e-5])
        # axs[1].set_ylim([-1e-3, 1e-3])
        # axs[1].set_ylim([-1e-6, 1e-6])
        # axs[1].set_xlim(grid.xbounds)
        
        # Adjust layout to prevent overlap
        plt.tight_layout()
        
        # Display the plot
        plt.show()

        # Close the plot
        plt.close()