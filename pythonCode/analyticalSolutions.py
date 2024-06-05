"""
Student ID: 31827379
"""

import numpy as np
    
def analyticalSolution1(i, icFunc, u, solver):
    """ 
    Analytical solution to:
        dphi\dt + u dphi\dx = iN phi

    """
    params = solver.model.params 
    X = solver.model.grid.X
    t = i*solver.dt
    
    return icFunc(X -  u*t)*np.exp(1j*params.N*t)

def analyticalSolution2(i, icFunc, solver):
    """ 
    Analytical solution to:
        dphi\dt = iN phi + iN dpsi/dx

    """
    
    # Get the default parameters for the problem.
    params = solver.model.params 
    X = solver.model.grid.X
    t = i*solver.dt
        
    b = np.pi/params.lx
    Q0 = params.psi0*b*(1 - 2*(b*X)**2)*np.exp(-(b*X)**2)
    
    term2 = -1j*params.N*Q0*(1j*params.N*np.sin(params.omega*t) + params.omega*
                           (np.cos(params.omega*t) - np.exp(1j*params.N*t)))
    term2 /= (params.omega**2 - params.N**2)
    
    return icFunc(X)*np.exp(1j*params.N*t) + term2
    
def analyticalSolution3(i, icFunc, solver):
    """ 
    Analytical solution to:
        dphi\dt = iN dpsi/dx
        
    """
    
    # Get the default parameters for the problem.
    params = solver.model.params
    X = solver.model.grid.X
        
    b = np.pi/params.lx 
    Q0 = params.psi0*b*(1 - 2*(b*X)**2)*np.exp(-(b*X)**2)
    
    term2 = -1j*params.N*Q0*(np.cos(params.omega*i*solver.dt) - 1)/params.omega
    
    return icFunc(X) + term2

def integrand4(t, X, u, params):
    
    b = np.pi/params.lx 
    a = X + u*t
    return (1j*params.N*
            params.psi0*b*(1 - 2*(a*b)**2)*np.sin(params.omega*t)*np.exp(-(a*b)**2))

def analyticalSolution4(i, solver):
    """ 
    Numerical integration solution to:
        dphi\dt + u dphi\dx = iN dpsi\dx
    """    
    
    # Wrapper for integrator result.
    return solver.getCustomEquation("integrator").I1

def integrand5(t, X, u, params):
    
    b = np.pi/params.lx 
    a = X + u*t
    return (1j*params.N*
            params.psi0*b*(1 - 2*(a*b)**2)*np.sin(params.omega*t)*
            np.exp(-(a*b)**2)*np.exp(-1j*params.N*t))
    
def analyticalSolution5(i, solver):
    """ 
    Numerical integration solution to:
        dphi\dt + u dphi\dx = iN phi + iN dpsi\dx
    """
    
    # Wrapper for integrator result.
    integrator = solver.getCustomEquation("integrator")
    params = solver.model.params
    
    return integrator.I1*np.exp(1j*params.N*i*solver.dt)