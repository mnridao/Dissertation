"""
Student ID: 31827379
"""
import numpy as np
import matplotlib.pyplot as plt

import pythonCode.helpers as helpers
import pythonCode.analyticalSolutions as analy

def defaultPlotter(grid, t=None, u=None):
    """ 
    """
    plt.figure(figsize=(10, 10))
    plt.plot(grid.X, grid.phi.real)
    plt.grid()
    plt.show(), plt.close()
    
def plotter1(grid, t=None, u=None):
    """ 
    """
    fontsize=15
    
    # Create a figure with 1 row and 2 columns
    fig, axs = plt.subplots(1, 2, figsize=(15, 7))
    
    # Plot the real part on the left subplot
    axs[0].plot(grid.X, grid.phi.real)
    axs[0].set_xlabel('X', fontsize=fontsize)
    axs[0].set_ylabel('b', fontsize=fontsize)
    axs[0].grid(True)
    axs[0].set_ylim([-1e-4, 1e-4])
    # axs[0].set_ylim([-1e-6, 1e-6])
    # axs[0].set_xlim(grid.xbounds)
    
    # Plot the imaginary part on the right subplot        
    axs[1].plot(grid.X, grid.phi.imag)
    axs[1].set_xlabel('X', fontsize=fontsize)
    axs[1].set_ylabel('Nw', fontsize=fontsize)
    axs[1].grid(True)
    # axs[1].set_ylim([-0.3e-5, 0.3e-5])
    axs[1].set_ylim([-1e-3, 1e-3])
    # axs[1].set_ylim([-1e-6, 1e-6])
    # axs[1].set_xlim(grid.xbounds)
    
    # Adjust layout to prevent overlap
    plt.tight_layout()
    
    # Display the plot
    plt.show()

    # Close the plot
    plt.close()

def plotWithAnalytical1(grid, t, u):
    
    # Source term.
    N = 0.02
    sPhiCoeff = 1j*N
    
    # Calculate analytical solution.
    phiAnaly = analy.analyticalSolution1(helpers.middleWaveIC, sPhiCoeff, u, grid.X, t)
    
    # Set the fontsize here.
    fontsize=15
    
    # Plot the thing.
    fig, axs = plt.subplots(1, 2, figsize=(15, 7))
    
    # Plot the real part on the left subplot
    axs[0].plot(grid.X, phiAnaly.real, 'k--')
    axs[0].plot(grid.X, grid.phi.real)
    axs[0].set_xlabel('X', fontsize=fontsize)
    axs[0].set_ylabel('b', fontsize=fontsize)
    axs[0].grid(True)
    axs[0].set_ylim([-1e-6, 1e-6])
    
    # Plot the imaginary part on the right subplot   
    axs[1].plot(grid.X, phiAnaly.imag, 'k--')     
    axs[1].plot(grid.X, grid.phi.imag)
    axs[1].set_xlabel('X', fontsize=fontsize)
    axs[1].set_ylabel('Nw', fontsize=fontsize)
    axs[1].grid(True)
    axs[1].set_ylim([-1e-6, 1e-6])
    
    # Adjust layout to prevent overlap
    plt.tight_layout()
    
    # Display the plot
    plt.show()

    # Close the plot
    plt.close()
    
def plotWithAnalytical2(grid, t, u=None):
    
    # Source term.
    N = 0.02
    sPhiCoeff = 1j*N
    
    # Calculate the analytic solution.
    ic = lambda x: (1+1j)*np.zeros_like(x, dtype=np.complex128)
    phiAnaly = analy.analyticalSolution2(ic, sPhiCoeff, grid.X, t)
    
    # Set the fontsize here.
    fontsize=15
    
    # Plot the thing.
    fig, axs = plt.subplots(1, 2, figsize=(15, 7))
    
    # Plot the real part on the left subplot
    axs[0].plot(grid.X/1e3, phiAnaly.real, 'k--')
    axs[0].plot(grid.X/1e3, grid.phi.real)
    axs[0].set_xlabel('X [km]', fontsize=fontsize)
    axs[0].set_ylabel('b', fontsize=fontsize)
    axs[0].grid(True)
    axs[0].set_ylim([-1e-4, 1e-4])
    
    # Plot the imaginary part on the right subplot   
    axs[1].plot(grid.X/1e3, phiAnaly.imag, 'k--')     
    axs[1].plot(grid.X/1e3, grid.phi.imag)
    axs[1].set_xlabel('X [km]', fontsize=fontsize)
    axs[1].set_ylabel('Nw', fontsize=fontsize)
    axs[1].grid(True)
    axs[1].set_ylim([-0.3e-5, 0.3e-5])
    
    # Adjust layout to prevent overlap
    plt.tight_layout()
    
    # Display the plot
    plt.show()

    # Close the plot
    plt.close()