"""
Student ID: 31827379
"""
import numpy as np

import matplotlib.pyplot as plt

#%%
if __name__ == "__main__":
    
    # Grid setup.
    x0 = 0
    xL = 160e3
    dx = 10e3
    # dx = 0.01
    nx = int((xL - x0)/dx)
    X = np.linspace(x0, xL, nx)
    
    # Initial condition.
    c = 0.8
    u = 5.
    
    phi0 = np.zeros_like(X, dtype=np.complex128)
        
    # Parameters.
    alpha = 0.5
    beta  = 1.
    gamma = 1
    N     = 2e-2
    
    # Source term.
    psi0 = 10 
    omega = 1.25e-4    
    dpsidx = lambda x, t: ((1/x - 2*x*(np.pi/xL)**2)*
                           psi0*(np.pi*x/xL)*np.sin(omega*t)*np.exp(-(np.pi*x/xL)**2))
    # s = lambda x, t: 1j*N*(1 + dpsidx(x, t))
    s = lambda x, t: 1j*N* + dpsidx(x, t)
    
    # Loop over time.
    endtime = 4e6
    dt = c*dx/u
    nt = int(np.ceil(endtime/dt))
    for t in range(nt):
        
        # phi2 = np.zeros_like(phi0, dtype=np.complex128)
        # for i in range(1, nx-1):
        #     si = s(i*dx, t*dt)
        #     phi2[i] = phi0[i] - c*(phi0[i] - phi0[i-1]) + si*phi0[i]
        
        # First iteration (k1 step).
        phi1 = np.zeros_like(phi0, dtype=np.complex128)
        
        # Update inlet bc (for now just zero?).
        # ...
        
        for i in range(1, nx-1):   
            si = s(i*dx, t*dt)
            phi1[i] = (phi0[i] - (1 - alpha*beta)*c*(phi0[i] - phi0[i-1]) + 
                        alpha*c*beta*phi1[i-1] - 
                        gamma*c*((2*phi0[i+1] - phi0[i] - phi0[i-1])/6 -
                                 (2*phi0[i] - phi0[i-1] - phi0[i-2])/6) + 
                        si*dt*(1 - alpha))
            phi1[i] /= (1 + alpha*beta*c - si*dt*alpha)
        
        # Update outlet bc (later).
        # ...
        
        # Second iteration (k2 step).
        phi2 = np.zeros_like(phi0, dtype=np.complex128)
        
        # Update inlet bc (for now just zero?).
        # ...
        
        for i in range(1, nx-1):
            si = s(i*dx, t*dt)
            phi2[i] = (phi0[i] - (1 - alpha)*c*(phi0[i] - phi0[i-1]) - 
                        alpha*c*(1 - beta)*(phi1[i] - phi1[i-1]) + 
                        alpha*c*beta*phi2[i-1] - gamma*(1 - alpha)*c*
                        ((2*phi0[i+1] - phi0[i] - phi0[i-1])/6 -
                        (2*phi0[i] - phi0[i-1] - phi0[i-2])/6) - gamma*alpha*c*
                        ((2*phi1[i+1] - phi1[i] - phi1[i-1])/6 -
                        (2*phi1[i] - phi1[i-1] - phi1[i-2])/6) + 
                        si*dt*(1 - alpha)*phi0[i])
            phi2[i] /= (1 + alpha*beta*c - si*dt*alpha)
        
        # Update outlet bc (later).
        # ...
        
        # Update phi0.
        phi0 = phi2.copy()
        # phi0 = phi1.copy()
        
        plt.plot(X, phi2.real)
        plt.ylim([-0.003, 0.003])
        plt.show()