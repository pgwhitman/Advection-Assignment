#Description: Create a graph of errors against resolution while keeping the Courant
#number constant for different linear advection schemes.

import numpy as np
import matplotlib.pyplot as plt

from advectionHelpers import cosBell, l2error, sinwave
from schemes import FTBS, CTCS, Semi_Lagrange

def resolutionErrors(init, schemes, outFile,
                           nx, yLim = [None,None]):
    '''
    Function to plot the error at endTime for a set of resolutions for a courant number of 0.5
    Values of c, endTime and xLim are chosen to take fixed values such that nt = 2*nx  and hence 
    nt will always take integer values for integer nx, ensuring the endTime and Courant number can 
    remain fixed for any value of nx.
    Parameters:
        nx (array): The number of points in space
        init (callable): Function to set the initial conditions.
        schemes (list of callable): List of functions to evaluate the schemes.
        c (float): Courant number
        outFile (string): The name of the output file.
        yLim (tuple): The minimum (yLim[0]) and maximum (yLim[1]) y values of
                      the graph. If None is specified then default limits
                      are used.
    Returns:
        None
    ''' 

    #Constant Parameters
    u = 1
    endTime = 1.0
    c = 0.5
    xLim = [0, 1]
    L = xLim[1] - xLim[0]

    # Colours for each line
    cols = ['blue', 'orange', 'green', 'cyan', 'magenta', 'brown', 'purple',
            'pink', 'red', 'gray', 'orange']

    #Initialise arrays for spacing and errors
    delta_x = np.zeros(len(nx))
    l2errors = np.zeros((len(nx), len(schemes)))

    # Loop over different resolutions/spatial steps 
    for n in range(len(nx)):
        # Derived parameters
        dx = (xLim[1] - xLim[0])/nx[n]
        delta_x[n] = dx
        x = np.arange(xLim[0], xLim[1], dx)
        #Calculate nt - must take integer values
        nt = nx[n]*u*endTime/(c*(xLim[1] - xLim[0]))
        #Check that nt is an integer (accounting for the fact python will store as float)
        if not np.isclose(nt, round(nt)):
            raise ValueError(f"nt must be an integer, but calculated nt = {nt}") 
        nt = int(round(nt))
        print('Courant number =', c)

        # Set the initial conditions for all schemes
        phi0 = init(x, L)
        # Compute analytic solution at endtime
        phiA = init((x - u*endTime)%(xLim[1] - xLim[0]), L) 

        #Loop over schemes and store errors
        for i in range(len(schemes)):
            # Copy initial conditions
            phiLevels = [phi0.copy()]
            # Loop through all time steps
            for t in range(nt):
                # Update phi
                schemes[i](phiLevels, c)
            # Calculate and store the l2 error
            l2errors[n, i] = l2error(phiLevels[-1], phiA)

    # Plot the errors against resolution for each scheme
    plt.figure()
    for i in range(len(schemes)):
        #Calculate the gradient ~ order
        rate = np.polyfit(np.log(delta_x), np.log(l2errors[:, i]), 1)[0]
        plt.loglog(delta_x, l2errors[:, i], cols[i], label = f"{schemes[i].__name__} (n={rate:.2f})")
    plt.ylim(yLim)
    plt.legend(fontsize=13)
    plt.xlabel(r'$\Delta x$', fontsize=13)
    plt.ylabel(r'$\ell_2$ error', fontsize=13)
    plt.savefig("plots/"+outFile)
    return None


#Schemes used
new_schemes = [FTBS, CTCS, Semi_Lagrange]

#Run script
resolutionErrors(lambda x, L: sinwave(x, L), new_schemes, "newSchemesSinResolutionErrorsnew.pdf",
                       nx = (np.arange(40, 400, 10)))