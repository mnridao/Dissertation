"""
Student ID: 31827379
"""

import numpy as np

from pythonCode.sourceterms import PsiParameters

def analyticalSolution1(icFunc, sPhiCoeff, u, X, t):
    """ 
    """
    
    return icFunc(X -  u*t)*np.exp(sPhiCoeff*t)

def analyticalSolution2(icFunc, sPhiCoeff, X, t):
    """ 
    """
    
    # Get the default parameters for the problem.
    params = PsiParameters()
    
    b = np.pi/params.lx
    term2 = (-sPhiCoeff*params.psi0*b*np.exp(-(b*X)**2)*(1 - 2*(b*X)**2)*
             (sPhiCoeff*np.sin(params.omega*t) + params.omega*np.cos(params.omega*t) -
             params.omega*np.exp(sPhiCoeff*t)))
    term2 /= (params.omega**2 - params.N**2)
    
    return icFunc(X)*np.exp(sPhiCoeff*t) + term2

def analyticalSolution3(icFunc, X, t):
    """ 
    """
    
    # Get the default parameters for the problem.
    params = PsiParameters()
        
    b = np.pi/params.lx 
    term2 = (1j*params.N*params.psi0*b*np.exp(-(b*X)**2)*(1 - 2*(b*X)**2)*
             (np.cos(params.omega*t - 1))/params.omega)
    
    return icFunc(X) - term2