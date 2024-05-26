"""
Student ID: 31827379
"""
import numpy as np
import matplotlib.pyplot as plt

def dpsidx2(x, t):
    """ 
    """
    # TODO: Think about how to avoid global parameters.
    
    b = np.pi/lx
    return 1j*N*psi0*b*(1 - 2*(b*x)**2)*np.sin(omega*t)*np.exp(-(b*x)**2)


def highOrderFlux(phi, i):
    
    # Check for boundary condition (assumes 0 upper gradient for now).
    if i == phi.shape[0]-1:
        phiPlus1 = phi[i]
    else:
        phiPlus1 = phi[i+1]
        
    return (2*phiPlus1 - phi[i] - phi[i-1])/6

def rkNIteration(phi0, phi1, c, nti, dx):
    
    # Second iteration (k2 step).
    phi2 = np.zeros_like(phi0, dtype=np.complex128)
            
    for i in range(1, nx):
        
        phi2[i] = (phi0[i] - 
                   (1 - alpha)*c*(phi0[i] - phi0[i-1]) - 
                   alpha*c*(1 - beta)*(phi1[i] - phi1[i-1]) + alpha*c*beta*phi2[i-1] - 
                   
                   gamma*(1 - alpha)*c*
                   (highOrderFlux(phi0, i) - highOrderFlux(phi0, i-1)) -
                    
                   gamma*alpha*c*
                   (highOrderFlux(phi1, i) - highOrderFlux(phi1, i-1)) +
                    
                   dt*(1 - alpha)*testS*phi0[i] + 
                   # 1j*N*dt*((1 - alpha)*testSFunc(i*dx, (t-1)*dt) + 
                   #          alpha*testSFunc(i*dx, t*dt))
                   dt*((1 - alpha)*testSFunc(i*dx, (t-1)*dt) + 
                            alpha*testSFunc(i*dx, t*dt))
                   )
        phi2[i] /= (1 + alpha*beta*c - testS*dt*alpha)
        
    return phi2

if __name__ == "__main__":
    
    # Grid setup.
    x0 = 0
    xL = 12e6
    dx = 10e3
    nx = int((xL - x0)/dx)
    X = np.linspace(x0, xL, nx)
    
    # Initial condition.
    c = 0.5
    u = 5.
    
    phi0 = np.zeros_like(X, dtype=np.complex128)
        
    # Parameters.
    alpha = 0.5
    beta  = 1.
    gamma = 1.
    
    # Source term.
    lx = 160e3
    N     = 0.02
    # N = 1e-3
    psi0 = 10 
    omega = 1.25e-4    
    dpsidx = lambda x, t: ((1/x - 2*x*(np.pi/lx)**2)*
                           psi0*(np.pi*x/lx)*np.sin(omega*t)*np.exp(-(np.pi*x/lx)**2))*1j*N

    # Loop over time.
    endtime = 4e6
    dt = c*dx/u
    dt = 10
    c = dt*c/dx
    nt = int(np.ceil(endtime/dt))

    testS = 1j*N
    testSFunc = dpsidx2
    for t in range(nt):
        
        phikMinus1 = phi0.copy()        
        for _ in range(0, 2):
        
            phik = rkNIteration(phi0, phikMinus1, c, t, dx)
            
            if _ == 1:
                break
            phikMinus1 = phik.copy()
        
        # Update phi0.
        phi0 = phik.copy()
        
        # imag -> w 
        # real -> b
        
        plt.figure(figsize=(10, 10))
        plt.plot(X, phik.real)
        plt.grid()
        plt.show(), plt.close()