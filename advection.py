###############################################################################
# Date Revised: 17 Sept 2025
# Author: Hilary Weller <h.weller@reading.ac.uk>
# Description: Basic tests of advection schemes
###############################################################################
import matplotlib.pyplot as plt
import numpy as np
import os
from advectionHelpers import cosBell, square_wave
from schemes import FTBS, CTCS, FTCS, Semi_Lagrange

# First make sure that the plots folder exists
if not os.path.exists("plots"):
    os.makedirs("plots")

def compareSchemeSnapshots(init, schemes, outFile,
                           nx = 40, nt = 10, u = 1, endTime = 0.4,
                           xLim = [0,1], yLim = [None,None]):
    '''
    Author: Hilary Weller <h.weller@reading.ac.uk>
    Compares snapshots of results of advection schemes using periodic boundaries
    Parameters:
        nx (int): The number of points in space
        init (callable): Function to set the initial conditions.
        schemes (list of callable): List of functions to evaluate the schemes.
        nt (int): The number of time steps.
        u (float): The wind speed
        endTime (float) : The final time
        outFile (string): The name of the output file.
        xLim (tuple): The minimum (xLim[0]) and maximim (xLim[1]) x-values of
                      the domain. xLim[1] is not part of the domain when periodic
        yLim (tuple): The minimum (yLim[0]) and maximum (yLim[1]) y values of
                      the graph. If [None,None] is specified then default limits
                      are used.
    Returns:
        None
    '''
    # Derived parameters
    dx = (xLim[1] - xLim[0])/nx
    print(dx)
    x = np.arange(xLim[0], xLim[1], dx)
    dt = endTime/nt
    c = nx*u*endTime/(nt*(xLim[1] - xLim[0])) # Courant number, dt*u/dx
    print('Courant number =', c)
    
    # Set the initial conditions
    phi0 = init(x)
    
    # Colours for each line
    cols = ['blue', 'orange', 'green', 'cyan', 'magenta', 'brown', 'purple',
            'pink', 'red', 'gray']
    
    # Find and plot the analytic solution
    phiA = init((x - u*endTime)%(xLim[1] - xLim[0]))
    plt.figure()
    plt.plot(x, phiA, 'k', label = "Analytic", linewidth=2)
    
    # Calculate the results of all of the schemes and plot
    for i in range(len(schemes)):
        # The solution, phi, is stored as a list of arrays at different time
        # steps, phiLevels. Initially just the initial conditions. The latest
        # time level array is always phiLevels[-1]
        phiLevels = [phi0.copy()]
        for n in range(nt):
            schemes[i](phiLevels, c)
        plt.plot(x, phiLevels[-1], cols[i], label = schemes[i].__name__)
    
    plt.xlim(xLim)
    plt.ylim(yLim)
    plt.xlabel('x', fontsize=15)
    plt.ylabel(r'$\phi$', fontsize=15)
    plt.legend(fontsize=13)
    plt.savefig("plots/"+outFile)

#Basic and new scheme lists
basic_schemes = [FTBS, CTCS, FTCS]
new_schemes = [FTBS, CTCS, Semi_Lagrange]

#Basic Scheme Comparison Plots

# compareSchemeSnapshots(lambda x: cosBell(x, a=0.1, b=0.5), basic_schemes, "basicSchemesCos.pdf",
#                        nx = 40, nt = 20, u = 1, endTime = 0.2, xLim = [0,1])

# compareSchemeSnapshots(lambda x: square_wave(x, a=0.1, b=0.5), basic_schemes, "basicSchemesSquare.pdf",
#                        nx = 40, nt = 20, u = 1, endTime = 0.2, xLim = [0,1])

# compareSchemeSnapshots(lambda x: cosBell(x, a=0.1, b=0.5), new_schemes, "newSchemesCos.pdf",
#                        nx = 40, nt = 20, u = 1, endTime = 0.2, xLim = [0,1])

# compareSchemeSnapshots(lambda x: square_wave(x, a=0.1, b=0.5), new_schemes, "newSchemesSquare.pdf",
#                        nx = 40, nt = 20, u = 1, endTime = 0.2, xLim = [0,1])

# compareSchemeSnapshots(lambda x: cosBell(x, a=0.1, b=0.5), new_schemes, "c1_12newSchemesCos.pdf",
#                        nx = 40, nt = 20, u = 2.8, endTime = 0.2, xLim = [0,1], yLim=[-0.05, 1.05])

# compareSchemeSnapshots(lambda x: square_wave(x, a=0.1, b=0.5), new_schemes, "c1_12newSchemesSquare.pdf",
#                        nx = 40, nt = 20, u = 2.8, endTime = 0.2, xLim = [0,1], yLim=[-0.05, 1.05])

compareSchemeSnapshots(lambda x: cosBell(x, a=0.1, b=0.5), new_schemes, "c1_12newSchemesCosdemo.pdf",
                       nx = 50, nt = 20, u = 1, endTime = 0.4, xLim = [0,1], yLim=[-0.05, 1.05])


