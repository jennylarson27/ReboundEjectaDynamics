# Debris Cloud Dynamics Rebound Package

This package allows for the simulation of debris clouds using N-body integrator Rebound. For specifics on the construction of this package, please reference Larson and Sarid (2021). This document should act as a manual for using this package. 

Necessary files/modules necessary to use package: 

- Python version of Rebound
- Python (any version works)
- Main script (binarytest.py)
- Input files (bodyparams.py; integratorparams.py)
- Effects files (additionaleffects.py; basicsolarsystem1015.py)
- Plotting script (plotejecta2.py)
 
	
	
# Simulation Set-up

All files must be saved in the same directory. The imput files bodyparams.py and bodyparams.py adjust the initial parameters of the simulation. The file bodyparams-sample.txt is a sample of the bodyparams.py input file with descriptions of each input variable. Similarly, integratorparams-sample.txt describes the input variables for integratorparams.py.

After the input files are set up with the initial parameters of the simulation, the simulation runs by calling the binarytest.py script with python in the terminal like so:

username:~/project_directory$ python binarytest.py

The output files with particle data at each time step will be saved to the project directory. It is possible to write a new plotting function to present the data; however, we have also included plotejecta2.py as a very simplistic plotting function that plots the particles in the x-y plane and the x-z plane. These two views may not be ideal for the scope of the user's study in which case we encourage them to use plotejecta2.py as a baseline from which to construct a new plotting function or to simply ignore plotejecta2.py and proceed with their own specific data analysis. For ease of use, plotejecta2-sample.txt explains the set-up of the script in order to plot the data. To run the plotejecta2.py script, call the file with python in the terminal:

username:~/project_directory$ python plotejecta2.py

We recommend plotting the data files in the same directory in which the simulation was run in order to keep track of which plots are associated with which simulation. After the plots are made, if the user wishes to save the plots in a separate directory, they can easily copy the plots to a separate directory (preferably named to reflect to which simulation those plots pertain).


# Function Descriptions: basicsolarsystem1015.py
This module contains functions to regulate the set-up and functions of running a basic simulation. 

**addpartcart** (sim, pos, vel, Nparts, [radmass, rho])
Adds particles to the simulation.

Parameters:
- sim: _Rebound_ simulation
- pos: (_tuple_) Particle position vectors
- vel: (_tuple_) Particle velocity vectors
- Nparts: (_int_) Number of particles in the simulation
- radmass [optional]: (_array_) Radii of particles
- rho [optional]: (_float_) Density of particles

Returns:
- sim: _Rebound_ simulation

**datasave** (sim, fname, Nplanets, Nparts)
Saves particle position and velocity data to .txt file.

Parameters:
- sim: _Rebound_ simulation
- fname: (_str_) Name of file where data will be saved
- Nplanets: (_int_) Number of non-particle objects in the system. Set value to 1 if the asteroid system does not contain a binary. Set value to 2 if there is a binary
- Nparts: (_int_) Number of particles in the system

Returns:
Saves file to current directory

**ejectaconpos** (mex, aT, bT, cT, h, lat, lon, beta, Nparts, tpos, mtarg, rhoT, shapemodel=False, vertfile=None)

Parameters:
- mex: (_float_) The total mass of material ejected from the impact site. Calculated based on the size of the resulting crater
- aT: (_float_) Semi-major axis of the target body
- bT: (_float_) First semi-minor axis of the target body
- cT: (_float_) Second semi-minor axis of the target body
- h: (_float_) Initial height of particles off the surface
- lat: (_float_) Latitude of impact location
- lon: (_float_) Longitude of impact location
- beta: (_float_) Angle of ejection from the vertical
- Nparts: (_float_) Number of particles in the system
- tpos: (_tuple_) Position vector of the target body
- mtarg: (_float_) Mass of the target body
- rhoT: (_float_) Target body bulk density
- shapemodel [optional]: (_boolean_) Default False. If True, will determine position of particles relative to a shape model rather than a sphere or ellipse
- vertfile [optional]: (_None_ or _string_) If shapemodel=False, then vertfile=None. This is the name of the vertices file used to define the shape model for the target body

Returns:
- pos: (_tuple_) Initial x-, y-, and z- components of the particle positions
- x: (_array_) x-coordinates of particles in the ejecta cone with the impact site at the origin
- y: (_array_) y-coordinates of particles in the ejecta cone with the impact site at the origin
- z: (_array_) z-coordinates of particles in the ejecta cone with the impact site at the origin
- thetaprime: (_array_) Angle defining the position of particles in the ejecta cone about the impact site
- r: (_array_) Distance from center of ejecta cone
- s: (_array_) Distance from particles to impact site
- pos0p: (_tuple_) x-, y-, and z- components of the impact site position
- vpos: (_tuple_) x-, y-, and z- components of the velocity unit vectors for particles

**ejectaconevel** (vpos, v0)
Calculates the initial velocity vector components of the particles

Parameters:
- vpos: (_tuple_) x-, y-, and z- components of the velocity unit vectors for particles
- v0: (_array_) Magnitude of initial particle velocities

Returns:
- vel: (_tuple_) x-, y-, and z- components of initial velocity vectors for particles

**rmdatasave** (rmx, rmy, rmz, fname)
Saves positions of removed particles. Particles that are removed are ones that have landed on the surface.

Parameters:
- rmx: (_array_) x-coordinates of removed particle positions
- rmy: (_array_) y-coordinates of removed particle positions
- rmz: (_array_) z-coordinates of removed particle positions
- fname: (_string_) File name to which removed particle positions are saved

Returns:
Saves .txt file to current directory

**rmland** (sim, Nparts, atarg, btarg, ctarg, landed, inttime, condit, axdir, axtilt, per, timestep, [abin, rotation, shapemodel, vertfile])
Removes particles that have landed on the surface of the target body and any potential binary components.

Parameters:
- sim: _Rebound_ simulation
- Nparts: (_float_) Number of particles in the system
- atarg: (_float_) Semi-major axis of the target body
- btarg: (_float_) First semi-minor axis of the target body
- ctarg: (_float_) Second semi-minor axis of the target body
- landed: (_float_) Total number of landed particles
- inttime: (_float_) Time step in the integration
- condit: (_string_) Name defining the conditions of the simulation. This becomes part of the file name where removed data are saved
- axdir: (_float_) Angle direction that the target body's axis of rotation points
- axtilt: (_float_) Angle at which the target body's axis of rotation is tilted from the vertical
- per: (_float_) Target body rotation period
- timestep: (_timestep_) Current time step in the simulation
- abin [optional]: (_float_) Semi-major axis of the binary component
- rotation [optional]: (_boolean_) True if target body is rotating. False if target body does not rotate.
- shapemodel [optional]: (_boolean_) True if implementing gravitational potential based on a shape model. False if using ellipsoidal or spherical gravity for the target body
- vertfile [optional]: (_None_ or _string_) File name for the target body shape model data

Returns:
- sim: _Rebound_ simulation
- landed: (_float_) Number of particles that have landed
- Nparts: (_float_) Number of particles in the system


**rotmatrix3d** (x, y, z, apha, phi, theta)
Rotates 3D array to produce a new 3D array

Parameters:
- x: (_array_) x-component of vector to be rotated
- y: (_array_) y-component of vector to be rotated
- z: (_array_) z-component of vector to be rotated
- alpha: (_float_) Angle to be rotated about the x-axis
- phi: (_float_) Angle to be rotated about the y-axis
- theta: (_float_) Angle to be rotated about the z-axis

Returns:
- xp: (_array_) x-component of rotated vector
- yp: (_array_) y-component of rotated vector
- zp: (_array_) z-component of rotated vector

**veldist** (sim, r, mtarg, rtarg, mex, mu, rhoT, [Cvps, Ctg, Ybar])
Calculates a particle velocity distribution based on scaling relations by Richardson (2011).

Parameters:
- sim: input _Rebound_ simulation
- r: (_array_) Initial particle distance from impact site
- mtarg: (_float_) Mass of target body
- rtarg: (_float_) Average radius of the target body
- mex: (_float_) Total mass excavated from the impact site
- mu: (_float_) Variable between 1/3 and 2/3 determined by whether impact is dominated by kinetic energy (mu=1/3) or momentum (mu=2/3)
- rhoT: (_float_) Bulk density of the target  body
- Cvps [optional]: (_float_) Proportionality constant
- Ctg [optional]: (_float_) Formation time constant
- Ybar[optional]: (_float_) Material strength of the target body


Returns:
- veff: (_array_) Particle speeds. Note - these are not the cartesian velocity vectors. To get the velocity vectors, multiply veff by the velocity unit vector.


# Function Descriptions: additionaleffects.py
This module contains functions pertaining to any additional effects necessary for different simulations.

_**Ellipsoidal Gravitational Potential**_

**addnetgrav** (sim, Mtarg, a1, b1, c1, Nplanets, Nparts, timestep, axdir, axtilt, [binary, rotation])
Adds ellipsoidal gravitational potential and subtracts the built-in spherical gravitational potential.

Parameters:
- sim:_Rebound_ simulation
- Mtarg: (_float_) Target body mass
- a1: (_float_) Semi-major axis of the target body
- b1: (_float_) First semi-minor axis of the target body
- c1: (_float_) Second semi-minor axis of the target body
- Nplanets: (_float_) Number of non-particle bodies in the system. This will be 1 if there is no binary and 2 if there is a binary component
- axdir: (_float_) aDirectional angle that the target body axis of rotation is tilted
- axtilt: (_float_) Angle that the target body axis of rotation is tilted from the vertical
- binary [optional]: (_boolean_) True if a binary component exists. False if there is no binary component.
- rotation [optional]: (_boolean_) True if the target body rotates. False if the target body is non-rotating.

Returns:
- sim: _Rebound_ simulation

**ellipgrav** (sim, surfpos, Mtarg, a1, b1, c1, [binary])

Returns:
- axellip
- ayellip
- azellip

**rmaddgrav** (a1, a2)

Returns:
- axnet
- aynet
- aznet



_**Shape Model Gravitational Potential**_

**shapemass** (vert, facet, M, [layers])

Returns:
- Mi
- Xg
- Yg
- Zg



_**Rotation of Primary Component**_

**rotpos** (sim, Nplanets, Nparts, axtilt, axdir, timestep, [per])

Returns:
- xp
- yp
- zp

**rotvel** (sim, per, lat, pos, axtilt, axdir)

Returns:
- vxrot
- vyrot
- vzrot



_**Binary Component**_

**binary** (sim, m, r, a, [e, i, periap, ascnode, f])

Returns:
- sim


_**Size Distribution**_

**partmass** (radii, rho)

Returns:
- mass

**partrad** (r, counts)

Returns:
- radii

**sizedist** (Nparts, Res, rmin, rmax, p)

Returns:
- r
- counts
- N


_**Radiation Pressure**_

**shadow** (sim, possun, Nplanets, Nparts, aview, bview, [rho])

Returns:
- accel

**solarradpress** (sim, Nplanets, possun, [rho])

Returns:
- acc
