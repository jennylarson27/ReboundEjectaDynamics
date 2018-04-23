### Main function for ejecta modeling ###
# Created 9/6/2017 #

import rebound
import basicsolarsystem as bss
import numpy as np
import matplotlib.pyplot as plt
import erroranalysis as ea
import additionaleffects as ae

from bodyparams import *
from integratorparams import *

sim = rebound.Simulation()
sim.units = (time, dist, mass)


# add target body
# asymmetric gravity?
if asymgrav == False:
	sim.add(r=rtarg,
			m=mtarg,
			x=0,
			y=0,
			z=0,
			vx=0,
			vy=0,
			vz=0,
			)

else:
	sim.add(r=rtarg,
			x=0,
			y=0,
			z=0,
			vx=0,
			vy=0,
			vz=0,
			)
	ar  = np.sqrt(a * b)
	rho = mtarg / ((4./3.) * np.pi * ar**3)



# add sun
sim.add(r=rsun,
		m=msun,
		x=possun[0],
		y=possun[1],
		z=possun[2],
		vx=velsun[0],
		vy=velsun[1],
		vz=velsun[2],
		)

Nplanets = 2

# Calculate hill radius
a = np.sqrt(possun[0] ** 2 + possun[1] ** 2 + possun[2] ** 2)
bound = a * ( (mtarg) / (3 * msun) )**(1./3.)


# add binary?
if binary == True:
	sim = ae.binary(sim,
					r=rbin,
					m=mbin,
					a=abin,
					e=ebin,
					i=ibin,
					periap=periapbin,
					ascnode=ascnodebin,
					f=fbin,
					)
	Nplanets += 1
	bound = a * ( (mbin + mtarg) / (3 * msun) )**(1./3.)


# add planets
for pl in planets:
	sim.add(pl, primary=sim.particles[1])
	Nplanets += 1




# add particles
Nparts = int(Nparts)

cpos = (0, 0, 0)#(-possun[0], -possun[1], -possun[2])
cvel = (0, 0, 0)#(-velsun[0], -velsun[1], -velsun[2])

pos, x1, y1, z1, unitv, r = bss.ejectaconepos(rtarg, hpart, theta, phi, beta, Nparts, cpos)
vel = bss.ejectaconevel(unitv, vpart, cvel)

if sizedist == True:

	# set up size distribution
	r, counts, N = ae.sizedist(Nparts, Res, rmin, rmax, p)
	radii = ae.partrad(r, counts, Nparts)

	sim = bss.addpartcart(sim, pos, vel, Nparts, radii)
else:
	sim = bss.addpartcart(sim, pos, vel, Nparts)


# add radiation pressure
if radpress == True:
	acc = ae.radforce(sim, Nparts, Nplanets, radii, rho)

	parts = sim.particles[Nplanets:]

	def addradforce(reb_sim):
		for i in np.arange(Nparts):
			x = parts[i].x - sim.particles[1].x
			y = parts[i].y - sim.particles[1].y
			z = parts[i].z - sim.particles[1].z

			dist  = np.sqrt(x ** 2 + y ** 2 + z ** 2)
			theta = np.arctan(y / x)
			phi   = np.arccos(z / dist)

			parts[i].ax += acc[i] * np.sin(phi) * np.cos(theta)
			parts[i].ay += acc[i] * np.sin(phi) * np.sin(theta)
			parts[i].az += acc[i] * np.cos(phi)

	sim.additional_forces = addradforce


# plot planet orbits
if binary == True:
	bel     = sim.particles[2]
	bela    = abin  # bel.a
	bele    = bel.e
	beli    = bel.inc
	belnode = bel.Omega
	belperi = bel.omega

	xbin, ybin, zbin = bss.orbplot(bela, bele, beli, belnode, belperi)


# define motion/times
motiontype = 'no'
if motiontype == 'Helio':
	body, sim = bss.helioorbits(sim)

elif motiontype == 'Bary':
	body, sim = bss.baryorbits(sim)
sim.calculate_com()
sim, times = bss.motion(sim, tmax, dt)


# set integrator
sim.integrator = integrator

E = []
T = []
inttime = 0

# set initial acceleration
parts = sim.particles[Nplanets:]
pax = []
pay = []
paz = []
for i in np.arange(Nparts):
	pax.append(parts[i].ax)
	pay.append(parts[i].ay)
	paz.append(parts[i].az)

# integrate
for i,time in enumerate(times):

	if asymgrav != False:
		for i in np.arange(Nparts):
			parts[i].ax = pax[i]
			parts[i].ay = pay[i]
			parts[i].az = paz[i]

		parts = sim.particles[Nplanets:]

		rprime, thetaprime, phiprime = ae.rthetphipart(sim, Nplanets)
		thetadeg = np.radians(np.arange(-180., 180.))
		phideg   = np.radians(np.arange(0., 180.))

		for theta in thetadeg:
			for phi in phideg:
				gx, gy, gz = ae.gravaccel(sim, a, b, c, theta, phi, thetaprime, phiprime, rprime, ar, rho, lmax, Nparts)

				def addgravforce(reb_sim):
					for i in np.arange(Nparts):
						parts[i].ax += gx[i]
						parts[i].ay += gy[i]
						parts[i].az += gz[i]

				sim.additional_forces = addgravforce


	print (time)

	if time == 0:
		pos0, vel0 = bss.init_posvel(sim, time)
		sim, rminds, landed, gone, binland, rmx, rmy, rmz = bss.rmparticles(sim, Nplanets, rtarg, outerlim, binr=rbin)


	# remove particles and update initial position
	sim, rminds, landed, gone, binland, rmx, rmy, rmz = bss.rmparticles(sim, Nplanets, rtarg, outerlim, landed, gone, binland, binr=rbin)
	pos0 = bss.updatepos(pos0, rminds)
	vel0 = bss.updatevel(vel0, rminds)
	totallost = landed + gone + binland
	print (landed, binland)

	# SAVE DATA
	bss.datasave(sim, 'particledata' + str(inttime) + condit + '.txt')
	bss.rmdatasave(rmx, rmy, rmz, 'rmparticledata' + str(inttime) + condit + '.txt')


	# plot results
#	if binary == True:
#		bss.plot3d(sim, 'didymos-' + str(inttime) + condit + '.png', Nplanets,
#				   rad=rtarg*scale, binary=rbin*scale, binorb=(xbin, ybin, zbin), size=size)

#		bss.plotxyxzplanes(sim, 'xyxzdidymos-' + str(hpart) + '-' + str(inttime) + condit + '.png',
#						   Nplanets,
#						   prad=radii*.0005,
#						   rad=rtarg*scale*3000,
#						   binary=rbin*scale,
#						   binorb=(xbin, ybin, zbin),
#						   #hill=bound,
#						   size=size,
#						   )

#		bss.plotxyxzplanes(sim, 'z-xyxzdidymos-' + str(hpart) + '-' + str(inttime) + condit + '.png',
#						   Nplanets,
#						   prad=radii*.0005,
#						   rad=rtarg*scale*3000,
#						   binary=rbin*scale,
#						   binorb=(xbin, ybin, zbin),
#						   #hill=bound,
#						   size=sizez,
#						   )


#	else:
#		bss.plot3d(sim, 'didymos-' + str(inttime) + condit + '.png',
#				   Nplanets,
#				   rad=rtarg * scale,
#				   size=size)

#		bss.plotxyxzplanes(sim, 'xyxzdidymos-' + str(hpart) + '-' + str(inttime) + condit + '.png',
#						Nplanets,
#						prad=radii*.005,
#						rad=rtarg * scale*3,
#						#hill=bound,
#						size=size,
#						)

#		bss.plotxyxzplanes(sim, 'z-xyxzdidymos-' + str(hpart) + '-' + str(inttime) + condit + '.png',
#						Nplanets,
#						prad=radii*.005,
#						rad=rtarg * scale*3,
#						#hill=bound,
#						size=sizez,
#						)




	E = ea.energycalc(sim, E)
	T.append(time)
	ea.xyplot(T, E, 'Timestep (days)', 'Total Energy', 'Total Energy at Each Timestep', 'EvsT.png')

	sim.integrate(time)
	inttime += 1


# change in energy
dE = ea.energydiff(E)
ea.xyplot(times, dE, 'Timestep (days)', 'dE', 'Change in Total Energy at Each Timestep', 'dEvsT.png')