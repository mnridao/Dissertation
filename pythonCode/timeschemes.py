"""
Student ID: 31827379
"""
import numpy as np
    
class ExplicitImplicitRK2:
    """ 
    """
    
    def __init__(self, sPhiCoeff, sFunc=None, alpha=0.5, beta=1., gamma=1.):
        
        self.alpha = alpha  # Temporal off-centering (0.5 for 2nd order)
        self.beta  = beta   # Implicit (1) / explicit (0) flag.
        self.gamma = gamma  # Flux limiter (for high-order spatial correction).
        
        # Source term coefficients and functions.
        self.sPhiCoeff = sPhiCoeff
        self.sFunc = sFunc if sFunc else lambda x, t: 0.
        
        # Boundary conditions - remove user option for now.
        # ...
        
    def __call__(self, grid, u0, dt, nti):
        
        # Calculate courant number.
        c = u0*dt/grid.dx
        # c = 0.5
        
        # Runge-Kutte iterations (2nd-order).
        phikMinus1 = grid.phi.copy()
        for it in range(0, 2):
            
            # Update the RK prediction.
            phik = self.rkNIteration(grid.phi, phikMinus1, c, nti, dt, grid.dx)
                        
            # To avoid unnecessary copy of array - idk if needed.
            if it == 1:
                break
            
            # Update previous iteration prediction.
            phikMinus1 = phik.copy()
        
        return phik
                
    def rkNIteration(self, phi0, phikMinus1, c, nti, dt, dx):
        
        # Domain length.
        nx = phi0.shape[0]
        
        # Initialise updated phi array.
        phik = np.zeros_like(phi0, dtype=phi0.dtype)
                
        for i in range(1, nx):
            
            phik[i] = (phi0[i] - 
                        (1 - self.alpha)*c*(phi0[i] - phi0[i-1]) - 
                        self.alpha*c*(1 - self.beta)*(phikMinus1[i] - phikMinus1[i-1]) + 
                        self.alpha*c*self.beta*phik[i-1] - 
                       
                        # Cubic correction terms at cell interface.
                        self.gamma*(1 - self.alpha)*c*
                        (self.highOrderFlux(phi0, i) - self.highOrderFlux(phi0, i-1)) - 
                       
                        self.gamma*self.alpha*c*
                        (self.highOrderFlux(phikMinus1, i) - self.highOrderFlux(phikMinus1, i-1)) +
                       
                        # Source terms.
                        self.sPhiCoeff*dt*(1 - self.alpha)*phi0[i] + 
                       
                        self.alpha*dt*self.sFunc(dx*i, dt*(nti-1)) + 
                        (1 - self.alpha)*dt*self.sFunc(dx*i, dt*nti)
                        )
            
            # TODO: Check for division by zero?
            phik[i] /= (1 + self.alpha*self.beta*c - dt*self.alpha*self.sPhiCoeff)
                
        return phik
    
    def highOrderFlux(self, phi, i):
                
        # Check for boundary condition (assumes 0 upper gradient for now).
        if i == phi.shape[0]-1:
            phiPlus1 = phi[i]
        else:
            phiPlus1 = phi[i+1]
            
        return (2*phiPlus1 - phi[i] - phi[i-1])/6