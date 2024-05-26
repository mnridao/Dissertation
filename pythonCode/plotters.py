"""
Student ID: 31827379
"""
import matplotlib.pyplot as plt

def defaultPlotter(grid):
    """ 
    """
    plt.figure(figsize=(10, 10))
    plt.plot(grid.X, grid.phi.real)
    plt.grid()
    plt.show(), plt.close()
    
def plotter1(grid, N=None):
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
    # axs[0].set_xlim(grid.xbounds)
    
    # Plot the imaginary part on the right subplot        
    axs[1].plot(grid.X, grid.phi.imag)
    axs[1].set_xlabel('X', fontsize=fontsize)
    axs[1].set_ylabel('Nw', fontsize=fontsize)
    axs[1].grid(True)
    axs[1].set_ylim([-0.3e-5, 0.3e-5])
    # axs[1].set_xlim(grid.xbounds)
    
    # Adjust layout to prevent overlap
    plt.tight_layout()
    
    # Display the plot
    plt.show()

    # Close the plot
    plt.close()
    