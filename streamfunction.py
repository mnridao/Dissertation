"""
Student ID: 31827379
"""

import numpy as np
import matplotlib.pyplot as plt

from pythonCode.sourceterms import DPsiDx, Psi
from pythonCode.grid import Grid1D 

if __name__ == "__main__":
    
    # Setup grid.
    x0 = 0
    xL = 12e6
    dx = 10e3
    nx = int((xL - x0)/dx)
    grid = Grid1D([x0, xL], nx)
        
    # Source terms.
    sPsi = Psi()
    sDPsiDx = DPsiDx()
    
    # Time loop parameters.
    endtime = 4e6
    dt = 10
    nt = int(np.ceil(endtime/dt))
    
    plotEveryN = 10
    for i in range(nt):
        
        # Calculate the gradient of the streamfunction for current timestep.
        dpsidxi = sDPsiDx(grid.X, i*dt)
        psi = sPsi(grid.X, i*dt)
        
        if i % plotEveryN == 0:
                
            # Create a figure and subplots
            fig, axs = plt.subplots(1, 2, figsize=(20, 10))
                      
            # Plot the streamfunction.
            axs[0].plot(grid.X / 1e3, psi)
            axs[0].set_xlabel("X [km]", fontsize=15)
            axs[0].set_ylabel(r'$\psi$', fontsize=15)
            axs[0].set_ylim([-4.5, 4.5])
            axs[0].grid()
                                    
            # Plot the gradient of the streamfunction.
            axs[1].plot(grid.X / 1e3, dpsidxi.imag/sPsi.params.N)
            axs[1].set_xlabel("X [km]", fontsize=15)
            axs[1].set_ylabel(r'$\frac{d\psi}{dx}$', fontsize=15)
            axs[1].set_ylim([-2e-4, 2e-4])
            axs[1].grid()
            
            # Add text to the second subplot showing the current time in hours
            hrs = i*dt / (60**2)
            axs[1].text(0.05, 0.95, f'Time: {hrs:.2f} hrs', fontsize=15,
                         horizontalalignment='left', verticalalignment='top', transform=axs[1].transAxes)
            
            # Adjust layout to prevent overlap
            plt.tight_layout()
            
            # Display the plot
            plt.show()