"""
Student ID: 31827379
"""

import numpy as np

from pythonCode.integrators import Trapezoidal
# from pythonCode.helpers import helpers
    
def advectionSolution(i, icFunc, u, solver):
    """ 
    """
    X = solver.model.grid.X
    t = i*solver.dt
    return icFunc(X - u*t)

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

def integrand3(t, X, u, params):
    """ 
    """
    b = np.pi/params.lx 
    return 1j*params.N*params.psi0*b*np.exp(-(b*X)**2)*(1 - 2*(b*X)**2)*np.sin(params.omega*t)

def analyticalSolution3_v2(i, solver):
    """ 
    """
    return solver.getCustomEquation("integrator").I1

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

def integrand4(t, endtime, X, u, params):
        
    b = np.pi/params.lx 
    a = X + u*t - u*endtime
    
    sol = (1j*params.N*
            params.psi0*b*(1 - 2*(a*b)**2)*np.sin(params.omega*t)*np.exp(-(a*b)**2))
    
    # Impose dirichlet boundary condition (phi=0) at lower boundary.
    sol[X+u*t <= u*endtime + X[0]] = 0
    
    return sol
    
def analyticalSolution4(i, icFunc, u, solver):
    """ 
    Numerical integration solution to:
        dphi\dt + u dphi\dx = iN dpsi\dx
    """     
    
    # Wrapper for integrator result.
    X = solver.model.grid.X
    t = i*solver.dt
    return icFunc(X - u*t) + solver.getCustomEquation("integrator").I1

def analyticalSolution4_v2(i, icFunc, u, solver, dt=None, nt=None):
    """ 
    """
    dt = dt if dt else solver.dt
    nt = nt if nt else solver.nt
    
    endtime = solver.dt*i
    integrator = Trapezoidal(integrand4, dt, nt, args=(
        endtime, solver.model.grid.X, u, solver.model.params), store=False)
    
    # Evaluate the integrator up until current endtime.
    integrator.integrateToLimit(endtime)
    
    X = solver.model.grid.X
    return icFunc(X - u*endtime) + integrator.I1
        
def integrand5(t, endtime, X, u, params):
    b = np.pi/params.lx 
    a = X + u*t - u*endtime
    
    sol = (1j*params.N*params.psi0*b*(1 - 2*(a*b)**2)*np.sin(params.omega*t)*
           np.exp(-(a*b)**2)*np.exp(-1j*params.N*t)*
           np.exp(1j*params.N*endtime))
    
    # Impose dirichlet boundary condition (phi=0) at lower boundary.
    sol[X+u*t <= u*endtime + X[0]] = 0
    
    return sol
    
def analyticalSolution5(i, icFunc, u, solver):
    """ 
    Numerical integration solution to:
        dphi\dt + u dphi\dx = iN phi + iN dpsi\dx
    """
    
    # Wrapper for integrator result.
    integrator = solver.getCustomEquation("integrator")
    params = solver.model.params
    
    X = solver.model.grid.X
    t = i*solver.dt
    return icFunc(X - u*t) - integrator.I1*np.exp(1j*params.N*t)

def analyticalSolution5_v2(i, icFunc, u, solver, dt=None, nt=None):
    """ 
    """
    dt = dt if dt else solver.dt
    nt = nt if nt else solver.nt
    
    endtime = solver.dt*i
    integrator = Trapezoidal(integrand5, dt, nt, args=(
        endtime, solver.model.grid.X, u, solver.model.params), store=False)
    
    # Evaluate the integrator up until current endtime.
    integrator.integrateToLimit(endtime)
    
    X = solver.model.grid.X
    return icFunc(X - u*endtime) + integrator.I1