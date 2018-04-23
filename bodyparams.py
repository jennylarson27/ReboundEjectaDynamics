### File of Parameters for Setting up System Bodies ###
# Created 9/18/2017 #

import numpy as np

# constant for converting units
au   = 1.496e8   # km in au
days = 86400.  # seconds in day
G    = 6.67e-11 * (days**2 / (1000 * au))

# add target body
mtarg = 4.56e11 # kg
rtarg = 0.4     #km

# use non-axisymmetric gravity
asymgrav = False
a = 30
b = 20
c = 10

# add sun
msun = 2e30   #kg
rsun = 6.96e5 #km
possun = (-8.124427114108884E-01*-au, 1.173956933680922E+00*-au, 6.629167645464454E-02*-au)
velsun = (-1.490113894696643E-02*-au, -3.601506749485048E-03*-au, 7.876917977423489E-04*-au)

# add binary
binary = False          # True if using binary, else False
# If not including a binary (binary=False), set the following values to zero.
mbin = 3e9   #kg
rbin = .075  #km
abin = 1.2   #km
ebin = 0.
ibin = 0.
periapbin=0.
ascnodebin=0.
fbin = 180.

# add planets
planets = ['Mercury',
           'Venus',
           'Earth',
           'Mars',
           'Jupiter',
           'Saturn',
           'Uranus',
           'Neptune',
		   ]


# add particles
Nparts = 1e4

sizedist = True       # True if using size distribution, else False
# If not using size distribution (sizedist=False), set next four values to zero.
Res  = 100        # Resolution of distribution (larger number=more particle sizes between rmin and rmax)
rmin = 1e-9#/au    # minimum radius
rmax = 1e-5#/au    # maximum radius
p    = -3.        # power of distribution
rho  = 2e12       # kg/km**3

hpart    = 1e-6   # particle initial height (km)
theta    = 270.     # longitude around equator (0-360)
phi      = 90.    # latitude from pole (0-180)
beta     = 90.    # opening angle of ejecta cone
vpart    = .5 * np.sqrt((2 * G * mtarg) / rtarg)#33.67  # initial velocity relative to target body
