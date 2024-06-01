"""
Student ID: 31827379
"""

import numpy as np
from scipy.integrate import quad_vec

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
    Q0 = params.psi0*b*(1 - 2*(b*X)**2)*np.exp(-(b*X)**2)
    
    term2 = -sPhiCoeff*Q0*(sPhiCoeff*np.sin(params.omega*t) + params.omega*
                           (np.cos(params.omega*t) - np.exp(sPhiCoeff*t)))
    term2 /= (params.omega**2 - params.N**2)
    
    return icFunc(X)*np.exp(sPhiCoeff*t) + term2
    
def analyticalSolution3(icFunc, sPhiCoeff,  X, t):
    """ 
    """
    
    # Get the default parameters for the problem.
    params = PsiParameters()
        
    b = np.pi/params.lx 
    Q0 = params.psi0*b*(1 - 2*(b*X)**2)*np.exp(-(b*X)**2)
    
    term2 = -sPhiCoeff*Q0*(np.cos(params.omega*t) - 1)/params.omega
    
    return icFunc(X) + term2



def integral4(t, X, u, params):
    
    b = np.pi/params.lx 
    a = X + u*t
    return params.psi0*b*(1 - 2*(a*b)**2)*np.sin(params.omega*t)*np.exp(-(a*b)**2)

def analyticalSolution4(icFunc, sPhiCoeff, u, X, t):
    """ 
    """
        
    # Get the default parameters for the problem.
    params = PsiParameters()
    
    # Numerically integrate for solution.
    res, err = quad_vec(integral4, 0, t, args=(X, u, params))
    
    return sPhiCoeff*res

def analyticalSolution5(icFunc, sPhiCoeff, u, X, t):
    """ 
    Full thing.
    """
    
    # Get the default parameters for the problem.
    params = PsiParameters()