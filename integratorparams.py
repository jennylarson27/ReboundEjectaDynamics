### File of Parameters for Integration ###
# Created 9/18/2017 #
from bodyparams import *

timeunits = 's'
dist      = 'm'
mass      = 'kg'

# additional effects
radpress = True       # True if using radiation pressure, else False

#motiontype = 'Bary'    # Choose either 'Bary' for barycentric motion or 'Helio' for heliocentric motion

dt   = 15.#/86400.            # time step (pick something smaller than a quarter of the smallest orbital period)
tmax = dt * 1e3             # maximum time to which we integrate

outerlim = 1.5e33

# non-axissymmetric gravity order
lmax = 2

integrator = 'ias15'   # Pick integrator to use (see REBOUND documentation about available integrators)

# plotting conditions
condit = '-rot'+str(int(Nparts))   # special conditions of simulation to include in file name
scale  = 90               # scaling factor of planets and binary for plotting
size   = 30                # size of system in plot
sizez  = 5             # size of zoomed in system
