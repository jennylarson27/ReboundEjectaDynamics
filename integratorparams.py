### File of Parameters for Integration ###
# Created 9/18/2017 #

time = 'day'
dist = 'km'
mass = 'kg'

# additional effects
radpress = True       # True if using radiation pressure, else False

#motiontype = 'Bary'    # Choose either 'Bary' for barycentric motion or 'Helio' for heliocentric motion

tmax = 15.             # maximum time to integrate to
dt   = 3e-4                # time step (pick something smaller than a quarter of the smallest orbital period)

outerlim = 1.5e30

# non-axissymmetric gravity order
lmax = 2

integrator = 'ias15'   # Pick integrator to use (see REBOUND documentation about available integrators)

# plotting conditions
condit = '-4rad'   # special conditions of simulation to include in file name
scale  = 90               # scaling factor of planets and binary for plotting
size   = 30                # size of system in plot
sizez  = 5             # size of zoomed in system