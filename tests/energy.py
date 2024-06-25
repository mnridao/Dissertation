"""
Student ID: 31827379
"""

import numpy as np
import matplotlib.pyplot as plt

import pythonCode.plotters as plotters
import pythonCode.helpers as helpers
import pythonCode.analyticalSolutions as analy
from pythonCode.sourceterms import PsiParameters

def runDifferentTimestepsEnergy(c, Ndt, u, endtime, solver, 
                                icKey, alphaKey=3):
    """ 
    """
    # Add energy to solver.
    solver.addCustomEquation("energy", helpers.energyCustom, 
                              args=(solver,), 
                              store=True)
    
    # Add analytic solution to solver (not really needed here).
    solver.addCustomEquation("analytic", analy.analyticalSolution1, 
                              args=(helpers.piecewiseWaveIC, u, solver), 
                              nres=solver.model.grid.phi.shape[0],
                              store=False)
    
    dts = Ndt/solver.model.params.N

    energies = []
        
    for j, dtj in enumerate(dts):
        
        # Modify the scheme parameter.
        if alphaKey == 1:
            solver.scheme.alpha = 0.5
        if alphaKey == 2:
            solver.scheme.alpha = np.maximum(0.5, 1 - 1/c)
        if alphaKey == 3:
            solver.scheme.alpha = np.maximum(0.5, 1 - 1/(c + c*params.N*solver.dt))
                
        # Set new timestep.
        solver.setNewTimestep(dtj, endtime)
        
        # Set new dx.
        dxj = u*solver.dt/c
        solver.model.grid.setNewGridSpacing(dxj)
        
        # Reset initial condition.
        solver.model.grid.phi = helpers.setInitialCondition(icKey, solver, u)
        
        print(f"c={u*solver.dt/solver.model.grid.dx}, \tdt={solver.dt}, \tdx={solver.model.grid.dx}, \talpha={solver.scheme.alpha}")
                            
        # Run the solver.
        solver.run(u)
        
        # Append energy.
        energies.append(solver.getCustomData("energy"))        
        
    return energies
    
if __name__ == "__main__":
    
    # Setup solver.
    params = PsiParameters()
    params.N = 1
    endtime = 5/0.01
    
    # endtime = 1000
    dt = 1
    solver = helpers.setupSolver(xbounds=[0, 30], dx=0.005, 
                                 endtime=endtime, dt=dt, schemeKey=1, 
                                 sPhiCoeffFlag=True, sFuncFlag=False,
                                 params=params)
    
    # Initial condition.
    u = 0.01
        
    # Change scheme parameters.
    c = u*solver.dt/solver.model.grid.dx
              
    c = 3
    u = 0.01
    endtime = 15/u
    
    params.N = 1
    Ndt = np.array([1, 3, 5])
    # Ndt = np.array([5])
    
    # solver.rootname = 'plots//c3N1_'
    
    # # Uncomment to see plots.
    # solver.plotResults=True
    # solver.plotter = plotters.plotWithAnalytical1
    # solver.plotEveryNTimesteps=1
    # solver.addCustomEquation("analytic", analy.analyticalSolution1, 
    #                           args=(helpers.piecewiseWaveIC, u, solver), 
    #                           nres=solver.model.grid.phi.shape[0],
    #                           store=False)
    # solver.model.grid.phi = helpers.setInitialCondition("a2", solver, u)
    # solver.run(u)
    
    # solver.snapshots = [200, 210, 250]
    # energies = runDifferentTimestepsEnergy(c, Ndt, u, endtime, solver, "a2", alphaKey)
        
    #%%
    
    c = [1.5, 3, 5]
    alphaLabels=["max(0.5, 1-1/(c + c$\Delta$tN))", "max(0.5, 1-1/c)", "0.5"]
    
    # Setup solver.
    params = PsiParameters()
    params.N = 1
    u = 0.01
    endtime = 15/u

    solver = helpers.setupSolver(xbounds=[0, 30], dx=0.005, 
                                 endtime=endtime, dt=1, schemeKey=1, 
                                 sPhiCoeffFlag=True, sFuncFlag=False,
                                 params=params)

    # Ndt for all c.
    alphaKey = 3
    Ndt = np.array([1, 3, 5])
    
    fig, axs = plt.subplots(1, 3, figsize=(30, 10))
    
    allEnergy=[]
    for i, ci in enumerate(c):
        
        energies = runDifferentTimestepsEnergy(ci, Ndt, u, endtime, solver, 
                                                "a2", alphaKey)
        # energies = allEnergy[i]
        allEnergy.append(energies)
        
        # Plot the results.
        energyA = helpers.energy(solver.getCustomData('analytic')) 
        
        for j, energy in enumerate(energies):
            time = np.linspace(0, endtime, energy.shape[0])
            
            axs[i].plot(time, energy, label=f'N$\Delta$t={Ndt[j]}')
        axs[i].plot([0, endtime], [energyA]*2, 'k--')
        # plt.plot([0, endtime], [150]*2, 'k--')
        
        if i == 0:
            axs[i].set_ylabel("Energy", fontsize=20)
        
        axs[i].set_yscale("log")
        axs[i].set_xlabel("Time [s]", fontsize=20)
        axs[i].set_title(f"c={ci}, alpha={alphaLabels[alphaKey-1]}", fontsize=20)
        axs[i].set_xlim([time[0], time[-1]])
        axs[i].tick_params(axis='both', which='both', labelsize=15)
        axs[i].grid()
        axs[i].legend(fontsize=15)
    # plt.savefig("plots//energies_alpha3.png")
    
    #%% Plot the energies.
    
    # palette = sns.color_palette("gnuplot", 10)[::5]
    # sns.set_theme(context='paper', style='white', palette=palette, 
    #               rc={'xtick.bottom': True, 'ytick.left': True})
    
    energyA = helpers.energy(solver.getCustomData('analytic'))
    
    fig = plt.figure(figsize=(10, 10))
    for i, energy in enumerate(energies):
        time = np.linspace(0, endtime, energy.shape[0])
        plt.plot(time, energy, label=f'N$\Delta$t={Ndt[i]}')
    plt.plot([0, endtime], [energyA]*2, 'k--')
    # plt.plot([0, endtime], [150]*2, 'k--')
        
    plt.yscale("log")
    plt.xlabel("Time [s]", fontsize=20)
    plt.ylabel("Energy", fontsize=20)
    plt.title(f"{c=}", fontsize=20)
    plt.xlim([time[0], time[-1]])
    plt.tick_params(axis='both', which='both', labelsize=15)
    plt.grid()
    plt.legend(fontsize=15)
    # plt.savefig("plots//energyc1N1.png")
    plt.show()
    plt.close()
    
    # plotters.plotWithAnalytical()