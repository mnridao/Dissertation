"""
Student ID: 31827379
"""
import matplotlib.pyplot as plt

def defaultPlotter(grid):
    """ 
    """
    plt.figure(figsize=(10, 10))
    plt.plot(grid.X, grid.phi.real)
    plt.grid()
    plt.show(), plt.close()