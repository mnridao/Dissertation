"""
Student ID: 31827379
"""

from scipy.stats import norm

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