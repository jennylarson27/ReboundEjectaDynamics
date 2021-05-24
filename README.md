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

Returns:
- sim

**datasave** (sim, fname, Nplanets, Nparts)

Returns:
nothing

**ejectaconpos** (mex, aT, bT, cT, h, lat, lon, beta, tilt, axdir, Nparts, tpos, mtarg, rtarg, vi, a, mi, rhoT, rhoI, mu, K1, Ybar=0., shapemodel=False, vertfile=None)

Returns:
- pos
- x
- y
- z
- thetaprime
- r
- s
- pos0p
- vpos

**ejectaconevel** (vpos, pos0p, v0, [per, rotation])

Returns:
- vel

**rmdatasave** (rmx, rmy, rmz, fname)

Returns:
nothing

**rmland** (sim, Nparts, atarg, btarg, ctarg, landed, inttime, condit, axdir, axtilt, per, timestep, [abin, rotation, shapemodel, vertfile])

Returns:
- sim
- landed
- Nparts


**rotmatrix3d** (x, y, z, apha, phi, theta)
Rotates 3D array to produce a new 3D array

Parameters:
- x: (_array_ or _float_) x-component of vector to be rotated
- y: (_array_ or _float_) y-component of vector to be rotated
- z: (_array_ or _float_) z-component of vector to be rotated
- alpha: (_float_) Angle to be rotated about the x-axis
- phi: (_float_) Angle to be rotated about the y-axis
- theta: (_float_) Angle to be rotated about the z-axis

Returns:
- xp: (_array_ or _float_) x-component of rotated vector
- yp: (_array_ or _float_) y-component of rotated vector
- zp: (_array_ or _float_) z-component of rotated vector

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
- Cvps [optional]: (_float_)
- Ctg [optional]: (_float_)
- Ybar[optional]: (_float_) Material strength of the target body


Returns:
- veff: (_array_) Particle speeds. Note - these are not the cartesian velocity vectors. To get the velocity vectors, multiply veff by the velocity unit vector.


# Function Descriptions: additionaleffects.py
This module contains functions pertaining to any additional effects necessary for different simulations.

_**Ellipsoidal Gravitational Potential**_

**addnetgrav** (sim, Mtarg, a1, b1, c1, Nplanets, Nparts, timestep, axdir, axtilt, [binary, rotation])

Returns:
- sim

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
