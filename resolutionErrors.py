import numpy as np
import matplotlib.pyplot as plt

from advectionHelpers import cosBell, l2error
from schemesBasic import FTBS, CTCS, semi_lagrange

def resolutionErrors(init, schemes, outFile,
                           nx, nt = 10, u = 1, endTime = 0.4,
                           xLim = [0,1], yLim = [None,None]):
    '''
    Function to plot the error at a set time for a set of resolutions
    Parameters:
        nx (int): The number of points in space
        init (callable): Function to set the initial conditions.
        schemes (list of callable): List of functions to evaluate the schemes.
        nt (int): The number of time steps.
        u (float): The wind speed
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

    # Need to loop over resolutions and schemes  
    delta_x = np.zeros(len(nx))
    l2errors = np.zeros((len(nx), len(schemes)))
    for n in range(len(nx)): #Need to make sure this works with how nx is defined
        # Derived parameters
        dx = (xLim[1] - xLim[0])/nx[n]
        delta_x[n] = dx
        x = np.arange(xLim[0], xLim[1], dx)
        dt = endTime/nt
        c = nx[n]*u*endTime/(nt*(xLim[1] - xLim[0])) # Courant number, dt*u/dx
        print('Courant number =', c)
        # Set the initial conditions for all schemes
        phi0 = init(x)
        # Compute analytic solution at endtime
        phiA = init((x - u*endTime)%(xLim[1] - xLim[0]))
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
        plt.loglog(delta_x, l2errors[:, i], cols[i], label = schemes[i].__name__)
        rate = np.polyfit(np.log(delta_x), np.log(l2errors[:, i]), 1)[0]
        print(rate)
    plt.ylim(yLim)
    plt.legend(fontsize=15)
    plt.xlabel(r'$\Delta x$', fontsize=13)
    plt.ylabel(r'$\ell_2$ error', fontsize=13)
    plt.savefig("plots/"+outFile)
    return None

new_schemes = [FTBS, CTCS, semi_lagrange]

#nx = np.array([35, 32, 30, 25, 22, 20, 15, 10])
nx = (np.arange(1, 200, 1))
print(nx)

resolutionErrors(lambda x: cosBell(x, a=0.1, b=0.5), new_schemes, "newSchemesCosResolutionErrors.pdf",
                       nx, nt = 20, u = 1, endTime = 0.4, xLim=[0,1])

