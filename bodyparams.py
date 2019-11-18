### File of Parameters for Setting up System Bodies ###
# Created 9/18/2017 #

import numpy as np

# Epoch = 2455355.5

# constant for converting units
au   = 1.496e11   # m in au
days = 86400.  # seconds in day
G    = 6.67e-11 #* (days**2) / (1000.**3)

# add target body
rtarg  = 75.#56.67*1e3#75.#400.     #km
mtarg  = 3e9#2e3 * (4./3.) * np.pi * (rtarg)**3#3e9#4.56e11 # kg

rotation = False
axtiltT = 0.      # axis of rotation beta
axdirT = 0.       # direction axis of rotation is tilted alpha
perT   = 12. * 360. #2.2593 * 60. * 60. # period of target body rotation (in seconds)

axtiltB = 5.      # axis of rotation beta
axdirB = 0.       # direction axis of rotation is tilted alpha
perB   = 2. * 60. * 60. # period of binary rotation (in seconds)

# use non-axisymmetric gravity
aT = 75.#56.67*1e3
bT = 75.#56.67*1e3
cT = 75.#56.67*1e3

aB = 400.
bB = 400.
cB = 400.

# add sun
msun = 2e30   #kg
#asys = 1.6446*au#2.9276*au  #1.6446
#rsun = 6.96e5 #km
possun = (-1.910547141790527E-01*-au, 1.207684573980673E+00*-au, 3.147139342751890E-02*-au)#(2.811386696837923E+00*-au, 1.665006888704921E+00*-au, 1.130150284142627E-01*-au) #(-1.910547141790527E-01*-au, 1.207684573980673E+00*-au, 3.147139342751890E-02*-au)
velsun = ((-1.732616227386132E-02*-au)/days, (-2.405596054545154E-03*-au)/days, (-1.029063688982998E-03*-au)/days)#((-3.098821517729514E-03*-au)/days, (6.981939328033965E-03*-au)/days, (4.638506014355656E-03*-au)/days)#((-1.732616227386132E-02*-au)/days, (-2.405596054545154E-03*-au)/days, (-1.029063688982998E-03*-au)/days)

# add binary
binary = True          # True if using binary, else False
# If not including a binary (binary=False), set the following values to zero.
mbin = 4.56e11#3e9   #kg
rbin = 400.#75.  #m
abin = 1.2e3   #m
ebin = 0.
ibin = 0.
periapbin=0.
ascnodebin=0.
fbin = 0.#np.radians(180.)


# impactor parameters - can be left blank to use defaults
vi   = 5.e3
di   = 35.
rhoi = 2e3
ex   = 2.



# add particles
Nparts = 1e4

sizedist = True       # True if using size distribution, else False
# If not using size distribution (sizedist=False), set next four values to zero.
Res  = Nparts/10.        # Resolution of distribution (larger number=more particle sizes between rmin and rmax)
rmin = 1e-4#/au    # minimum radius
rmax = 1e0#/au    # maximum radius
p    = -3.        # power of distribution
rho  = 2e3      # kg/m**3

hpart   = 1   # particle initial height (mm)
lon     = 90.     # longitude around equator (0-360)
lat     = 0.    # latitude from equator (-90-90) equator=0
beta    = 45.    # opening angle of ejecta cone
vinit   = .05#.7 * np.sqrt((2 * G * (mtarg)) / (hpart + bT))#33.67  # initial velocity relative to target body
print ('vinit', vinit)
