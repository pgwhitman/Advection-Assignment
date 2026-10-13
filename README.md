Advection assignment. Hilary Weller <h.weller@reading.ac.uk>
Code to solve the advection equation.

run_all.py            : Run this to run all scripts to produce all plots in report

advection.py          : Run this to save plots of snapshots of the results
                        for schemes CTCS, FTBS, FTCS and Semi-Lagrangian
growingErrors.py      : Run this to create a graph of errors growing over time
                        for each scheme
resolutionErrors.py   : Run this to create a graph of errors against resolution 
                        for each scheme
stability_test.py     : Run this to create a graph of errors against Courant number
                        for each scheme

advectionHelpers.py   : Functions for initialising and calculating l2 error norm
schemes.py            : FTCS, CTCS, FTBS and Semi-Lagrangian advection schemes

