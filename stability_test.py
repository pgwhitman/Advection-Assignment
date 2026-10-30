#Description: Creates a graph of error against Courant number for different linear advection schemes

import numpy as np
import matplotlib.pyplot as plt

from advectionHelpers import cosBell, l2error
from schemes import FTBS, CTCS, Semi_Lagrange

def CourantErrors(init, schemes, outFile,
                           nx, nt = 10, u = 1, endTime = 0.4,
                           xLim = [0,1], yLim = [None,None]):
    '''
    Function to plot the error at a set time for a set of resolutions
    Parameters:
        nx (int): The number of points in space
        nt (int): The number of timesteps
        init (callable): Function to set the initial conditions.
        schemes (list of callable): List of functions to evaluate the schemes.
        nt (int): The number of time steps.
        u (array): The wind speed - varied in order to vary courant number
        endTime (float) : The final time
        outFile (string): The name of the output file.
        xLim (tuple): The maximum and minimun spacial dimensions
        yLim (tuple): The minimum (yLim[0]) and maximum (yLim[1]) y values of
                      the graph. If None is specified then default limits
                      are used.
    Returns:
        None
    ''' 

    # Colours for each line
    cols = ['blue', 'orange', 'green', 'cyan', 'magenta', 'brown', 'purple',
            'pink', 'red', 'gray', 'orange']

    # Initialise arrays
    courant = np.zeros(len(nx))
    l2errors = np.zeros((len(nx), len(schemes)))

    #Loop over Courant numbers/u
    for n in range(len(nx)): 
        # Derived parameters
        dx = (xLim[1] - xLim[0])/nx[n]
        x = np.arange(xLim[0], xLim[1], dx)
        dt = endTime/nt
        c = nx[n]*u*endTime/(nt*(xLim[1] - xLim[0])) # Courant number, dt*u/dx
        courant[n] = c
        print('Courant number =', c)

        # Set the initial conditions for all schemes
        phi0 = init(x)

        # Compute analytic solution at endtime
        phiA = init((x - u*endTime)%(xLim[1] - xLim[0]))

        #Loop over schemes
        for i in range(len(schemes)):
            # The solution for this scheme
            phiLevels = [phi0.copy()]
            # Loop through all time steps
            for t in range(nt):
                # Update phi
                schemes[i](phiLevels, c)
            # Calculate the l2 error
            l2errors[n, i] = l2error(phiLevels[-1], phiA)

    # Plot the error for this scheme
    plt.figure()
    for i in range(len(schemes)):
        plt.semilogy(courant, l2errors[:, i], cols[i], label = schemes[i].__name__)
    plt.ylim(yLim)
    plt.legend(fontsize=15)
    plt.xlabel(r'$c$', fontsize=13)
    plt.ylabel(r'$\ell_2$ error', fontsize=13)
    plt.savefig("plots/"+outFile)
    return None

new_schemes = [FTBS, CTCS, Semi_Lagrange]

nx = np.arange(15, 155, 10)

CourantErrors(lambda x: cosBell(x, a=0.1, b=0.5), new_schemes, "newSchemesCosStabilityErrorstighternonintdx.pdf",
                       nx = nx, nt = 50, u = 1, endTime = 1.0, xLim=[0,1], yLim=[10e-6,10e12])

