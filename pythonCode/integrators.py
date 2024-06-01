"""
Student ID: 31827379
"""

from abc import ABC, abstractmethod

import numpy as np

class IntegrationStepper(ABC):
    """ 
    """
    
    def __init__(self, f, h, nt, args, store=True):
        """ 
        """
        # Integrand calculation.
        self.f = f 
        self.args = args
        
        # Things for the integration step.
        self.h  = h
        self.nt = nt
        self.y0 = self.f(0, *self.args)
        
        # Accumulated integral.
        self.I1 = np.zeros_like(self.y0)
        self.I = np.zeros(shape=(self.nt, *self.y0.shape), 
                          dtype=self.y0.dtype) if store else None
        self.store = store
            
    def __call__(self, solver, i):
        
        # Just in case ig.
        if i == 0:
            return
        
        # Calculate the area of the subinterval using trapezoidal method.
        self.calculateSubInterval(self.h*i)
        
        # Store the value at current grid point.
        if self.store:
            self.I[i-1, ...] = self.I1
    
    @abstractmethod
    def calculateSubInterval(self, t):
        pass

class Trapezoidal(IntegrationStepper):
    """ 
    """
    
    def __init__(self, f, h, nt, args=(), store=True):
        
        # Initialise parent class.
        super().__init__(f, h, nt, args, store)
                        
    def calculateSubInterval(self, t):
        """ 
        """
        
        # Approximate subinterval as trapezoid and calculate area.
        y1 = self.f(t, *self.args)
        self.I1 += 0.5*self.h*(self.y0 + y1)
        
        # Update previous subinterval.
        self.y0 = y1
        
class Simpsons(IntegrationStepper):
    """ 
    """
    
    def __init__(self, f, h, nt, args, store=True):
        
        # Initialise parent class.
        super().__init__(f, h, nt, args, store)
    
    def calculateSubInterval(self, t):
        pass
    
#%%
if __name__ == "__main__":
    
    import numpy as np
    import matplotlib.pyplot as plt
    
    integrand = lambda t : np.sin(t)
    
    a = np.ones(shape=10)
    nt = 1000
    
    h = 0.001
    integrator = Trapezoidal(integrand, h, nt)

    for i in range(1, nt+1):
        
        # Update the integral.
        integrator(i)
    
    x = np.linspace(0, h*(nt+1), nt+1)
    plt.figure(figsize=(10, 10))
    plt.plot(integrator.I)
    plt.plot(-np.cos(x)+1, 'k--')
    plt.grid()
    plt.show()
    plt.close()
    
#%% Try with one of the problem integrals.
    
    from pythonCode.sourceterms import PsiParameters
    from pythonCode.grid import Grid1D 

    def integrand1(t, u, X0, params):
        
        b = np.pi/params.lx
        X = u*t + X0
        
        return (1j*params.N*params.psi0*b*(1 - 2*(b*X)**2)*np.exp(-(b*X)**2)*
                np.sin(params.omega*t))

    # Setup grid.
    x0 = 0
    xL = 12e6
    dx = 10e3
    nx = int((xL - x0)/dx)
    grid = Grid1D([x0, xL], nx)

    # Other parameters.
    u  = 100
    params = PsiParameters()
    args = (u, grid.X, params)
    
    # Time loop parameters.
    endtime = 4e6
    dt = 10
    nt = int(np.ceil(endtime/dt))
    
    integrator = Trapezoidal(integrand1, dt, nt, args=args, store=False)
    
    for i in range(1, nt+1):
        
        # Update the integral.
        integrator(None, i)
        
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
        axs[0].plot(grid.X, integrator.I1.real)
        axs[0].set_xlabel('X', fontsize=fontsize)
        axs[0].set_ylabel('b', fontsize=fontsize)
        axs[0].grid(True)
        # axs[0].set_ylim([-1e-4, 1e-4])
        axs[0].set_ylim([-1e-6, 1e-6])
        # axs[0].set_xlim(grid.xbounds)
        
        # Plot the imaginary part on the right subplot        
        axs[1].plot(grid.X, integrator.I1.imag)
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