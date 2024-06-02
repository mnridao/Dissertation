"""
Student ID: 31827379
"""
import matplotlib.pyplot as plt

def defaultPlotter(solver):
    """ 
    """
    
    plt.figure(figsize=(10, 10))
    plt.plot(solver.grid.X, solver.grid.phi.real)
    plt.grid()
    plt.show(), plt.close()
    
def plotter1(solver):
    """ 
    """
    fontsize=15
    
    # Create a figure with 1 row and 2 columns
    fig, axs = plt.subplots(1, 2, figsize=(15, 7))
    
    # Plot the real part on the left subplot
    axs[0].plot(solver.grid.X, solver.grid.phi.real)
    axs[0].set_xlabel('X', fontsize=fontsize)
    axs[0].set_ylabel('b', fontsize=fontsize)
    axs[0].grid(True)
    # axs[0].set_ylim([-1e-4, 1e-4])
    axs[0].set_ylim([-1e-6, 1e-6])
    # axs[0].set_xlim(grid.xbounds)
    
    # Plot the imaginary part on the right subplot        
    axs[1].plot(solver.grid.X, solver.grid.phi.imag)
    axs[1].set_xlabel('X', fontsize=fontsize)
    axs[1].set_ylabel('Nw', fontsize=fontsize)
    axs[1].grid(True)
    axs[1].set_ylim([-0.3e-5, 0.3e-5])
    # axs[1].set_ylim([-1e-3, 1e-3])
    # axs[1].set_ylim([-1e-6, 1e-6])
    # axs[1].set_xlim(grid.xbounds)
    
    plt.tight_layout()
    plt.show()
    plt.close()

def plotWithAnalytical1(solver):
    
    yboundsRe = [-1e-6, 1e-6]
    yboundsIm = [-1e-6, 1e-6]
    
    plotWithAnalytical(solver, yboundsRe, yboundsIm)
    
def plotWithAnalytical2(solver):
    
    yboundsRe = [-1e-4, 1e-4]
    yboundsIm = [-0.3e-5, 0.3e-5]
    
    plotWithAnalytical(solver, yboundsRe, yboundsIm)
    
def plotWithAnalytical3(solver):
    
    yboundsRe = [-1e-4, 1e-4]
    yboundsIm = [-0.08, 0.08]
    
    plotWithAnalytical(solver, yboundsRe, yboundsIm)
        
def plotWithAnalytical4(solver):
    
    # Same as plotWithAnalytical1 for now.
    yboundsRe = [-1e-6, 1e-6] 
    yboundsIm = [-1e-6, 1e-6]
    
    plotWithAnalytical(solver, yboundsRe, yboundsIm)
    
def plotWithAnalytical5(solver):
    
    # Same as plotWithAnalytical1 for now.
    yboundsRe = [-1e-6, 1e-6] 
    yboundsIm = [-1e-6, 1e-6]
    
    plotWithAnalytical(solver, yboundsRe, yboundsIm)
        
def plotWithAnalytical(solver, yboundsRe, yboundsIm):
    """ 
    """
    # Analytic solution (from numerical integrator).
    phiAnaly = solver.getCustomData("analytic")
        
    # Set the fontsize here.
    fontsize=15
        
    # Plot the thing.
    fig, axs = plt.subplots(1, 2, figsize=(15, 7))
    
    # Plot the real part on the left subplot
    axs[0].plot(solver.model.grid.X/1e3, phiAnaly.real, 'k--')
    axs[0].plot(solver.model.grid.X/1e3, solver.model.grid.phi.real)
    axs[0].set_xlabel('X [km]', fontsize=fontsize)
    axs[0].set_ylabel('b', fontsize=fontsize)
    axs[0].grid(True)
    axs[0].set_ylim(yboundsRe)
    
    # Plot the imaginary part on the right subplot   
    axs[1].plot(solver.model.grid.X/1e3, phiAnaly.imag, 'k--')     
    axs[1].plot(solver.model.grid.X/1e3, solver.model.grid.phi.imag)
    axs[1].set_xlabel('X [km]', fontsize=fontsize)
    axs[1].set_ylabel('Nw', fontsize=fontsize)
    axs[1].grid(True)
    axs[1].set_ylim(yboundsIm)

    plt.tight_layout()
    plt.show()
    plt.close()