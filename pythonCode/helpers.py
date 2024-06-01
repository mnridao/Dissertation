"""
Student ID: 31827379
"""

from scipy.stats import norm

import pythonCode.timeschemes as ts

def setScheme(key, args=()):
    
    if key == 1 or "ImpExRK":
        return ts.ExplicitImplicitRK2(*args)
    
    # Free to add more.

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