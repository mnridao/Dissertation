"""
Student ID: 31827379
"""

import numpy as np

def smoothWave(X, u, t):
    """ 
    Student ID: 31827379
    
    Solution for the square wave.
        
    Inputs
    -------
    X  : float or np.array
         Domain locations where the solution is evaluated.
    a  : float
         Domain lower bound of the wave at t = 0.
    b  : float
         Domain upper bound of the wave at t = 0.
    u  : float
         Advection speed.
    t  : float
         Time.
    
    Returns
    -------
    np.array
        The analytical solution at time t
    """
    
    a = 0.1 
    b = 0.5
    
    smoothWave = lambda x, a, b : 0.5 * (1 - np.cos(2*np.pi * ((x - a) / (b - a))))
    return np.array([0. if x <= a or x > b else smoothWave(x, a, b) for x in (X - u*t)])

def advectionWithDecay1D(X, u0, t, s, phi0):
    """ 
    """
    return phi0(X - u0*t, u0, 0)*np.exp(s*t)

if __name__ == "__main__":
    
    import matplotlib.pyplot as plt
    
    # Grid setup.
    x0 = 0
    xL = 12e6
    dx = 10e3
    # dx = 0.01
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
    # N = 1e-6
    psi0 = 10 
    omega = 1.25e-4    
    dpsidx = lambda x, t: ((1/x - 2*x*(np.pi/lx)**2)*
                           psi0*(np.pi*x/lx)*np.sin(omega*t)*np.exp(-(np.pi*x/lx)**2))
        
    A = lambda x: (psi0*(np.pi/lx)*np.exp(-(np.pi*x/lx)**2) - 
                   2*x*(np.pi/lx)**2*psi0*(np.pi*x/lx)*np.exp(-(np.pi*x/lx)**2))
    
    
    # Loop over time.
    endtime = 4e6
    dt = 10
    nt = int(np.ceil(endtime/dt))
    for i in range(1, nt):
        
        # Current time.
        t = i*dt
        
        # Analytical solution maybe.
        phi1 = np.zeros_like(phi0, dtype=np.complex128)
        
        phiNew = ((phi0 + 1j*N*omega*A(X - c*t)/(omega**2 - N**2))*np.exp(1j*N*t) -
                  1j*N*A(X)*(omega*np.cos(omega*t) + 1j*N*np.sin(omega*t))/
                  (omega**2 - N**2)).real
        
        plt.plot(X, phiNew)
        # plt.plot(X, A(X - c*t))
        plt.show()