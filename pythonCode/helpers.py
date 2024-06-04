"""
Student ID: 31827379
"""

from scipy.stats import norm
import numpy as np

import pythonCode.timeschemes as ts
from pythonCode.grid import Grid1D 
from pythonCode.solver import Solver 
from pythonCode.model import Model
from pythonCode.sourceterms import PsiParameters, DPsiDx

## MISC ##
def setScheme(key, args=()):
    
    if key == 1 or "ImpExRK2":
        return ts.ExplicitImplicitRK2(*args)
    
    # Free to add more.

def setupSolver(xbounds, dx, endtime, dt, schemeKey, sPhiCoeffFlag, sFuncFlag, 
                params=None):
    """ 
    """
    
    # Setup grid.
    x0, xL = xbounds
    nx = int((xL - x0)/dx)
    grid = Grid1D(xbounds, nx)
    
    # Setup model.
    params = params if params else PsiParameters()
    model = Model(grid, params)
    
    # Set time scheme.
    schemeArgs = [1j*params.N*sPhiCoeffFlag]
    if sFuncFlag:
        schemeArgs.append(DPsiDx())
    scheme = setScheme("ImpExRK2", *(schemeArgs, ))
    
    # Solver parameters.
    nt = int(np.ceil(endtime/dt))
    return Solver(model, scheme, dt, nt)

## ERROR METRICS ##
def l2(phiN, phiA, dx):
    """ 
    """
    return (np.sqrt(np.sum(dx*np.power(phiN - phiA, 2)))/
            np.sqrt(np.sum(dx*np.power(phiA, 2))))

def l2ReCustomEqn(i, solver):
    """ 
    """
    phiN = solver.model.grid.phi.real
    phiA = solver.getCustomData("analytic").real
    
    return l2(phiN, phiA, solver.model.grid.dx)

def l2ImCustomEqn(i, solver):
    """ 
    """
    phiN = solver.model.grid.phi.imag
    phiA = solver.getCustomData("analytic").imag
    
    return l2(phiN, phiA, solver.model.grid.dx)

def rmse(phiN, phiA):
    """ 
    """
    N = phiN.shape[0]
    return np.sqrt(np.sum((phiA - phiN)**2)/N)

def rms(phi):
    """ 
    """
    return np.sqrt(np.mean(phi**2))

## INITIAL CONDITIONS ##
def waveIC(mu, sigma, X):
    """ 
    """
    
    # Generate 1D Gaussian distribution with given mean and standard deviation.
    gaussian = norm(loc=mu, scale=sigma)
    
    return gaussian.pdf(X)

def middleWaveIC(X):
    """ 
    """
    
    mu = (X[-1] - X[0])/2
    sigma = 0.05*X[-1]
    
    return waveIC(mu, sigma, X)

def complexZerosIC(X):
    """ 
    """
    
    return np.zeros_like(X, dtype=np.complex128)