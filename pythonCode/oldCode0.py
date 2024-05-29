"""
Student ID: 31827379
"""
import numpy as np

import matplotlib.pyplot as plt

#%%
if __name__ == "__main__":
    
    # Grid setup.
    x0 = 0
    xL = 12e6
    dx = 10e3
    # dx = 0.01
    nx = int((xL - x0)/dx)
    X = np.linspace(x0, xL, nx)
    
    # Initial condition.
    # c = 0.5
    u = 5.
    
    phi0 = np.zeros_like(X, dtype=np.complex128)
        
    # Parameters.
    alpha = 0.5
    beta  = 0.
    gamma = 1.
    
    # Source term.
    lx = 160e3
    N     = 0.02
    # N = 1e-3
    psi0 = 10 
    omega = 1.25e-4    
    dpsidx = lambda x, t: ((1/x - 2*x*(np.pi/lx)**2)*
                           psi0*(np.pi*x/lx)*np.sin(omega*t)*np.exp(-(np.pi*x/lx)**2))

    # Loop over time.
    endtime = 4e6
    # dt = c*dx/u
    dt = 10
    c = u*dt/dx
    nt = int(np.ceil(endtime/dt))
    for t in range(nt):
        
        # First iteration (k1 step).
        phi1 = np.zeros_like(phi0, dtype=np.complex128)
        
        # Update inlet bc (for now just zero).
        # ...
        
        for i in range(1, nx):
            if i == nx-1:
                phi0Plus1 = phi0[i] 
            else:
                phi0Plus1 = phi0[i+1]
            
            phi1[i] = (phi0[i] - (1 - alpha*beta)*c*(phi0[i] - phi0[i-1]) + 
                        alpha*c*beta*phi1[i-1] - 
                        gamma*c*((2*phi0Plus1 - phi0[i] - phi0[i-1])/6 -
                                 (2*phi0[i] - phi0[i-1] - phi0[i-2])/6) + 
                        
                        # Something wrong with source term i think :(
                        dt*(1 - alpha)*1j*N*phi0[i] + 
                        # 1j*N*dt*dpsidx(i*dx, t*dt)
                        1j*N*dt*((1 - alpha)*dpsidx(i*dx, (t-1)*dt) + 
                                  alpha*dpsidx(i*dx, t*dt))
                        )
            phi1[i] /= (1 + alpha*beta*c - dt*alpha*1j*N)
                
        # Second iteration (k2 step).
        phi2 = np.zeros_like(phi0, dtype=np.complex128)
        
        # Update inlet bc (for now just zero).
        # ...
        
        for i in range(1, nx):
            if i == nx-1:
                phi0Plus1 = phi0[i] 
                phi1Plus1 = phi1[i] 
            else:
                phi0Plus1 = phi0[i+1]
                phi1Plus1 = phi1[i+1]
            
            phi2[i] = (phi0[i] - (1 - alpha)*c*(phi0[i] - phi0[i-1]) - 
                        alpha*c*(1 - beta)*(phi1[i] - phi1[i-1]) + 
                        alpha*c*beta*phi2[i-1] - gamma*(1 - alpha)*c*
                        ((2*phi0Plus1 - phi0[i] - phi0[i-1])/6 -
                        (2*phi0[i] - phi0[i-1] - phi0[i-2])/6) - gamma*alpha*c*
                        ((2*phi1Plus1 - phi1[i] - phi1[i-1])/6 -
                        (2*phi1[i] - phi1[i-1] - phi1[i-2])/6) + 
                        
                        # Something wrong with source term i think :(
                        dt*(1 - alpha)*1j*N*phi0[i] + 
                        1j*N*dt*((1 - alpha)*dpsidx(i*dx, (t-1)*dt) + 
                                  alpha*dpsidx(i*dx, t*dt))
                        )
            phi2[i] /= (1 + alpha*beta*c - 1j*N*dt*alpha)
                
        # Update phi0.
        phi0 = phi2.copy()
        # phi0 = phi1.copy()
        
        # imag -> w 
        # real -> b
        
        plt.plot(X, phi2.real)
        # plt.ylim([-1e-6, 1e-6])
        plt.show()