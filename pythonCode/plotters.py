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
    X = solver.model.grid.X
    axs[0].plot(X, solver.model.grid.phi.real)
    axs[0].set_xlabel('X', fontsize=fontsize)
    axs[0].set_ylabel('b', fontsize=fontsize)
    axs[0].grid(True)
    # axs[0].set_ylim([-1e-6, 1e-6])
    
    # Plot the imaginary part on the right subplot        
    axs[1].plot(X, solver.model.grid.phi.imag)
    axs[1].set_xlabel('X', fontsize=fontsize)
    axs[1].set_ylabel('Nw', fontsize=fontsize)
    axs[1].grid(True)
    # axs[1].set_ylim([-0.3e-5, 0.3e-5])
    
    plt.tight_layout()
    plt.show()
    plt.close()

def plotWithAnalytical1(solver, figname='', save=False):
    
    # yboundsRe = [-1e-6, 1e-6]
    # yboundsIm = [-1e-6, 1e-6]
    
    yboundsRe = [-1, 1]
    yboundsIm = [-1, 1]
    
    plotWithAnalytical(solver, yboundsRe, yboundsIm, figname, save)
    
def plotWithAnalytical2(solver, figname='', save=False):
    
    # yboundsRe = [-1e-4, 1e-4]
    # yboundsIm = [-0.3e-5, 0.3e-5]
    
    yboundsRe = [-10, 10]
    yboundsIm = [-10, 10]
    
    plotWithAnalytical(solver, yboundsRe, yboundsIm, figname, save)
    
def plotWithAnalytical3(solver, figname='', save=False):
    
    yboundsRe = [-1e-4, 1e-4]
    yboundsRe = [-10, 10]
    
    # yboundsIm = [-0.08, 0.08]
    yboundsIm = [-10, 10]
    
    plotWithAnalytical(solver, yboundsRe, yboundsIm, figname, save)
        
def plotWithAnalytical4(solver, figname='', save=False):
    
    # Same as plotWithAnalytical1 for now.
    yboundsRe = [-1e-6, 1e-6] 
    # yboundsIm = [-1e-6, 1e-6]
    yboundsIm = [-1e-2, 1e-2]
    
    plotWithAnalytical(solver, yboundsRe, yboundsIm, figname, save)
    
def plotWithAnalytical5(solver, figname='', save=False):
    
    # Same as plotWithAnalytical1 for now.
    yboundsRe = [-1e-6, 1e-6] 
    yboundsIm = [-1e-6, 1e-6]
    
    plotWithAnalytical(solver, yboundsRe, yboundsIm, figname, save)
        
def plotWithAnalytical(solver, yboundsRe, yboundsIm, figname='', save=False):
    """ 
    """
    # Analytic solution (from numerical integrator).
    phiAnaly = solver.getCustomData("analytic")
        
    # Set the fontsize here.
    fontsize=15
        
    # Plot the thing.
    fig, axs = plt.subplots(1, 2, figsize=(15, 7))
    
    # Plot the real part on the left subplot
    X = solver.model.grid.X
    
    # axs[0].plot(solver.model.grid.X/1e3, phiAnaly.real, 'k--')
    axs[0].plot(X/1e3, phiAnaly.real, 'k--')
    axs[0].plot(solver.model.grid.X/1e3, solver.model.grid.phi.real)
    axs[0].set_xlabel('X [km]', fontsize=fontsize)
    axs[0].set_ylabel('b', fontsize=fontsize)
    axs[0].grid(True)
    axs[0].set_ylim(yboundsRe)
    # axs[0].set_xlim(solver.model.grid.xbounds)
    
    # Plot the imaginary part on the right subplot   
    axs[1].plot(X/1e3, phiAnaly.imag, 'k--')     
    axs[1].plot(solver.model.grid.X/1e3, solver.model.grid.phi.imag)
    axs[1].set_xlabel('X [km]', fontsize=fontsize)
    axs[1].set_ylabel('Nw', fontsize=fontsize)
    axs[1].grid(True)
    axs[1].set_ylim(yboundsIm)
    # axs[1].set_xlim(solver.model.grid.xbounds)

    plt.tight_layout()
    if save:
        plt.savefig(figname)
    
    plt.show()
    plt.close()
    
def plotWithoutAnalytical(solver, figname='', save=False):
    """ 
    """            
    plot(solver.model.grid.X, solver.model.grid.phi, 
         figname=figname, save=save)
    
def plot(X, phi, phiA=None, figname='', save=False, 
         yboundsIm=None, yboundsRe=None):
    """ 
    """
    
    # Set the fontsize here.
    fontsize=15
        
    # Plot the thing.
    fig, axs = plt.subplots(1, 2, figsize=(15, 7))
    
    # Plot the real part on the left subplot    
    axs[0].plot(X, phi.real)
    if phiA:
        axs[0].plot(X, phiA.real)
    
    axs[0].set_xlabel('X', fontsize=fontsize)
    axs[0].set_ylabel('b', fontsize=fontsize)
    axs[0].grid(True)
    if yboundsRe:
        axs[0].set_ylim(yboundsRe)
    # axs[0].set_xlim(solver.model.grid.xbounds)
    
    # Plot the imaginary part on the right subplot    
    axs[1].plot(X, phi.imag)
    if phiA:
        axs[1].plot(X, phiA.imag)
    axs[1].set_xlabel('X', fontsize=fontsize)
    axs[1].set_ylabel('Nw', fontsize=fontsize)
    axs[1].grid(True)
    if yboundsIm:
        axs[1].set_ylim(yboundsIm)
    # axs[1].set_xlim(solver.model.grid.xbounds)

    plt.tight_layout()
    if save:
        plt.savefig(figname)
    
    plt.show()
    plt.close()