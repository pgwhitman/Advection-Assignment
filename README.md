Advection assignment. Hilary Weller <h.weller@reading.ac.uk>
Code to solve the advection equation.
advection.py          : Run this to save plots of snapshots of the results
                        for schemes CTCS, FTBS and FTCS
growingErrors.py      : Run this to create a graph of errors growing over time
                        for each scheme.

advectionHelpers.py   : Functions for initialising and calculating l2 error norm
schemesBasic.py    : FTCS, CTCS and FTBS advection schemes


Finalise this at the end:
Requires python 3.14.7

Requirements in requirement.txt

If dependencies required, run: 

python -m venv .venv

.venv\Scripts\activate

pip install -r requirements.txt
