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

Determines initial positions of particles in an ejecta cone above the surface.

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

**ellipgrav** (sim, surfpos, Mtarg, a1, b1, c1)

Calculates the acceleration due to a gravity field around an ellipsoidal body.

Parameters:
- sim: _Rebound_ simulation
- surfpos: (_tuple_) Cartesian coordinates (x, y, z) of the particle positions relative to the surface of the body
- Mtarg: (_float_) Mass of the body for which gravity is being calculated
- a1: (_float_) Semi-major axis of the body for which gravity is being calculated
- b1: (_float_) First semi-minor axis of the body for which gravity is being calculated
- c1: (_float_) Second semi-minor axis of the body for which gravity is being calculated

Returns:
- axellip: (_array_) x-component of the acceleration felt by particles due to the ellipsoidal gravity field
- ayellip: (_array_) y-component of the acceleration felt by particles due to the ellipsoidal gravity field
- azellip: (_array_) z-component of the acceleration felt by particles due to the ellipsoidal gravity field

**rmaddgrav** (a1, a2)

Calculates a1-a2. This is the net acceleration experienced by a particle around a body. This function is primarily used to subtract away the default spherical gravitational acceleration provided by _Rebound_ and to replace it with the ellipsoidal gravitational acceleration.

Parameters:
- a1: (_tuple_) Cartesian coordinates (x, y, z) of the acceleration to be added (typically the ellipsoidal gravitational acceleration)
- a2: (_tuple_) Cartesian coordinates (x, y, z) of the acceleration to be subtracted (typically the spherical gravitational acceleration)

Returns:
- axnet: (_array_) x-component of the net gravitational acceleration to be applied to each particle
- aynet: (_array_) y-component of the net gravitational acceleration to be applied to each particle
- aznet: (_array_) z-component of the net gravitational acceleration to be applied to each particle



_**Shape Model Gravitational Potential**_

**shapemass** (vert, facet, M, [layers])

Calculates the mass distribution within a shape model in order to calculate a gravity field based on the observed shape model.

Parameters:
- vert: (_string_) File name for the vertex locations that make up the body
- facet: (_string_) File name for the facet locations that make up the body

Returns:
- Mi: (_array_) Mass of each volume making up the body shape
- Xg: (_array_) x-component of volume locations in the shape model
- Yg: (_array_) y-component of volume locations in the shape model
- Zg: (_array_) z-component of volume locations in the shape model



_**Rotation of Primary Component**_

**rotpos** (sim, Nplanets, Nparts, axtilt, axdir, timestep, [per])

Calculates the positions of particles above the rotating body.

Parameters:
- sim: _Rebound_ simulation
- Nplanets: (_float_) Number of non-particle bodies in the simulation
- Nparts: (_float_) Number of particles in the simulation
- axtilt: (_float_) Angle that the axis of rotation tilts from the vertical
- axdir: (_float_) Direction the axis of rotation points away from the vertcal
- timestep: (_float_) Current time step in the simulation
- per [optional]: (_float_) Rotation period

Returns:
- xp: (_array_) x-coordinate of the particle locations post-rotation
- yp: (_array_) y-coordinate of the particle locations post-rotation
- zp: (_array_) z-coordinate of the particle locationa post-rotation

**rotvel** (sim, per, lat, pos, axtilt, axdir)

Calculates the rotational velocity that is added to the initial particle ejection velocity

Parameters:
- sim: _Rebound_ simulation
- per: (_float_) Rotation period
- lat: (_array_) Initial particle latitude
- pos: (_tuple_) Cartesian coordinates (x, y, and z) of the initial particle positions
- axtilt: (_float_) Angle that the axis of rotation tilts away from the vertical
- axdir: (_float_) Direction the axis of rotation points away from the vertical

Returns:
- vxrot: (_array_) x-component of the rotational velocity added to the initial particle velocities
- vyrot: (_array_) y-component of the rotational velocity added to the initial particle velocities
- vzrot: (_array_) z-component of the rotational velocity added to the initial particle velocities



_**Binary Component**_

**binary** (sim, m, r, a, [e, i, periap, ascnode, f])

Adds a secondary component to the system

Parameters:
- sim: _Rebound_ simulation
- m: (_float_) Mass of secondary body
- r: (_float_) Radius of the secondary
- a: (_float_) Distance of the secondary from the primary
- e [optional]: (_float_) Default is 0. Eccentricity of the secondary body
- i [optional]: ((_float_) Default is 0. Inclination of the secondary body
- periap [optional]: (_float_) Default is 0. Argument of periapse of the secondary body
- ascnode [optional]: (_float_) Default is 0. Ascending node of the secondary body
- f [optional]: (_float_) Default is 0. True anomaly of the secondary body

Returns:
- sim: _Rebound_ simulation


_**Size Distribution**_

**partmass** (radii, rho)

Calculates the mass of a particle based on the density and radius.

Parameters:
- radii: (_array_) Radii of the particles in the simulation
- rho: (_float_) Density of ejected particles

Returns:
- mass: (_array_) Mass of ejected particles

**partrad** (r, counts)

Creates array of particle radii.

Parameters:
- r: (_array_) Array of possible radius values
- counts: (_array_) Number particles assigned to each radius value

Returns:
- radii: (_array_) Radii for each particle

**sizedist** (Nparts, Res, rmin, rmax, p)

Calculates a size distribution of radii to be assigned to particles.

Parameters:
- Nparts: (_float_) Number of particles
- Res: (_int_) Resolution of size distribution. The larger this number is, the more unique radii will be produced. This should be less than the number of particles. If Res=Nparts, each particle will have a unique radius. It is recommended to select a value smaller than the number of particles.
- rmin: (_float_) Minimum particle radius
- rmax: (_float_) Maximum particle radius
- p: (_float_) Power of the power law distribution

Returns:
- r: (_array_) Possible radius values
- counts: (_array_) Exact number of particles with each possible radius value (Array must have the same shape as r)
- N: (_float_) Rounded number of particles with each possible radius value.


_**Radiation Pressure**_

**shadow** (sim, possun, Nplanets, Nparts, aview, bview, [rho])

Determines when particles are in the shadow of the body and adds a temporary acceleration to counteract the solar radiation pressure acceleration.

Parameters:
- sim: _Rebound_ simulation
- possun: (_tuple_) Cartesian coordinates (x, y, z) of the sun's position relative to the target body
- Nplanets: (_float_) Number of non-particle bodies in the system
- Nparts: (_float_) Number of particles in the system
- aview: (_float_) Semi-major axis of the target body
- bview: (_float_) First semi-minor axis of the target body
- rho [optional]: Density of the particles. Default is 2e3 kg/m^3.

Returns:
- accel: (_tuple_) Cartesian vector components (x, y, z) of the acceleration caused by radiation pressure including particles within the shadow (acceleration is zero for shadowed particles)

**solarradpress** (sim, Nplanets, possun, [rho])

Calculates the acceleration experienced by the particles due to radiation pressure from the sun.

Parameters:
- sim: _Rebound_ simulation
- Nplanets: (_float_) Number of non-particle bodies in the system
- possun: (_tuple_) Cartesian coordinates (x, y, z) of the sun's position relative to the target body
- rho [optional]: (_float_) The density of the particles. Default is 2e3 kg/m^3.

Returns:
- acc: (_tuple_) Cartesian vector components (x, y, z) of the acceleration caused by radiation pressure
