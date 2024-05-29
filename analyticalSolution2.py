"""
Student ID: 31827379
"""

import numpy as np
import matplotlib.pyplot as plt

from pythonCode.timeschemes import ExplicitImplicitRK2
from pythonCode.grid import Grid1D 
from pythonCode.sourceterms import PsiParameters, DPsiDx
from pythonCode.solver import Solver 
import pythonCode.plotters as plotters
import pythonCode.helpers as helpers
import pythonCode.analyticalSolutions as analy

def analyticalSolution(icFunc, sPhiCoeff, X, t):
    
    # Get the default parameters for the problem.
    params = PsiParameters()
    
    b = np.pi/params.lx
    term2 = (-1j*params.N*params.psi0*b*np.exp(-(b*X)**2)*(1 - 2*(b*X)**2)*
             (1j*params.N*np.sin(params.omega*t) + params.omega*np.cos(params.omega*t) -
             params.omega*np.exp(sPhiCoeff*t)))
    term2 /= (params.omega**2 - params.N**2)
    
    return icFunc(X)*np.exp(sPhiCoeff*t) + term2
    
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
    # sPhiCoeff = 1j*0
    sFunc = DPsiDx()
        
    # Set time scheme.
    scheme = ExplicitImplicitRK2(sPhiCoeff, sFunc)
    u = 0.
    
    # Solver parameters.
    endtime = 4e6
    dt = 10
    nt = int(np.ceil(endtime/dt))
    solver = Solver(grid, scheme, dt, nt)
    solver.plotter = plotters.plotWithAnalytical2
    
    # Initial condition.
    solver.grid.phi = np.zeros_like(solver.grid.X, dtype=np.complex128)
    
    solver.run(u)
    
    # #%%
    
    # ic = lambda x: (1 + 1j)*0
    
    # # test = []
    # for i in range(nt):
        
    #     # Calculate analytical solution for the current time step.
    #     phiAnaly = analyticalSolution(ic, sPhiCoeff, grid.X, i*dt)
        
    #     # test.append(phiAnaly.imag)
                                
    #     # Create a figure and subplots
    #     fig, axs = plt.subplots(1, 2, figsize=(20, 10))
                  
    #     # Plot the real.
    #     axs[0].plot(grid.X / 1e3, phiAnaly.real)
    #     axs[0].set_xlabel("X [km]", fontsize=15)
    #     axs[0].set_ylabel("b", fontsize=15)
    #     axs[0].set_ylim([-1e-4, 1e-4])
    #     axs[0].grid()
                                
    #     # Plot the imag.
    #     axs[1].plot(grid.X / 1e3, phiAnaly.imag)
    #     axs[1].set_xlabel("X [km]", fontsize=15)
    #     axs[1].set_ylabel("wN", fontsize=15)
    #     axs[1].set_ylim([-0.3e-5, 0.3e-5])
    #     axs[1].grid()
        
    #     # Add text to the second subplot showing the current time in hours
    #     hrs = i*dt / (60**2)
    #     axs[1].text(0.05, 0.95, f'Time: {hrs:.2f} hrs', fontsize=15,
    #                   horizontalalignment='left', verticalalignment='top', transform=axs[1].transAxes)
        
    #     # Adjust layout to prevent overlap
    #     plt.tight_layout()
        
    #     # Display the plot
    #     plt.show()