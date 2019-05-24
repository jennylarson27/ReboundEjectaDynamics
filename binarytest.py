import rebound
import numpy as np
import time as ctime
import basicsolarsystem as bss
import matplotlib.pyplot as plt
import erroranalysis as ea
import additionaleffects as ae
import velocitydist as vd
import functools
import os

from bodyparams import *
from integratorparams import *



# add particles
def conesetup(*rest):
	# Set up cone
	global Nparts

	Nparts = int(Nparts)

	cpos = (0, 0, 0)
	cvel = (0, 0, 0)


	pos, x1, y1, z1, phiang, thetang, r, s = bss.ejectaconepos(rtarg,
															hpart,
															lat,
															lon,
															beta,
															Nparts,
															cpos,
															)



	# Start particle velocities with rotation


	# Produce size distribution and add particles
	if sizedist == True:

		# set up size distribution
		r, counts, N = ae.sizedist(Nparts,
								   Res,
								   rmin,
								   rmax,
								   p,
								   )

		radii = ae.partrad(r,
						   counts,
						   )


		# if veldist == 'power':
		# 	vpart = vd.power(vmin, vmax, Nparts, counts, Res, sizedist)
		#
		# elif veldist == 'zmodel':
		# 	vpart = vd.maxwellz(s, Z, alpha)
		#
		# else:
	vpart = np.full((Nparts), vmax)

	vel = bss.ejectaconevel(phiang,
							thetang,
							lat,
							lon,
							vpart,
							cvel,
							#Nparts,
							)

	vel = (vel[0], vel[1], vel[2])

	if sizedist == True:
		return pos, vel, Nparts, radii
	else:
		return pos, vel, Nparts


	# else:
	# 	if veldist == 'power':
	# 		vpart = vd.power(vmin, vmax, Nparts, counts, Res, sizedist)
	#
	# 	elif veldist == 'zmodel':
	# 		vpart = vd.maxwellz(s, Z, alpha)
	#
	# 	else:
	# 		vpart = np.full((Nparts), vmax)
	#
	# 	vel = bss.ejectaconevel(phiang,
	# 							thetang,
	# 							lat,
	# 							lon,
	# 							vpart,
	# 							cvel,
	# 							#Nparts,
	# 							)
	#
	# 	vel = (vel[0], vel[1], vel[2])
	#
	# 	return pos, vel, Nparts






def type1radpress(*rest):
	global Nplanets

	# add radiation pressure for type 1
	acc = ae.solarradpress(sim, Nplanets)

	def addradforce(reb_sim):

		partax = acc[0]
		partay = acc[1]
		partaz = acc[2]

		for i in np.arange(Nplanets, Nplanets+Nparts):

			p = sim.particles[int(i)]

			p.ax += partax[int(i-Nplanets)]
			p.ay += partay[int(i-Nplanets)]
			p.az += partaz[int(i-Nplanets)]

	sim.additional_forces = addradforce

	return sim, Nplanets



### Run all the functions ###

if sizedist == True:
	pos, vel, Nparts, radii = conesetup()

elif sizedist == False:
	pos, vel, Nparts = conesetup()




### For running less particles ###
def sim_setup(*rest):
	global Nparts
	global sim
	global pos
	global vel
	global radii
	global binx
	global biny
	global binz
	global Nplanets
	global landed
	global binland
	global rminds

	sim = rebound.Simulation()
	sim.units = (timeunits, dist, mass)

	Nparts_0 = Nparts

	if rotation == True:
		sim.add(m=mtarg,
				x=0,
				y=0,
				z=0,
				vx=0,
				vy=0,
				vz=0,
				)
	else:
		sim.add(r=rtarg,
				m=mtarg,
				x=0,
				y=0,
				z=0,
				vx=0,
				vy=0,
				vz=0,
				)
	Nplanets = 1

	if binary == True:
		if rotation == True:
			sim.add(m=mbin,
					a=abin,
					e=ebin,
					inc=ibin,
					primary=sim.particles[0],
					)
		else:
			sim.add(r=rbin,
					m=mbin,
					a=abin,
					e=ebin,
					inc=ibin,
					primary=sim.particles[0],
					)
		Nplanets += 1

	# add massless sun
	sim.add(x=possun[0],
			y=possun[1],
			z=possun[2],
			vx=velsun[0],
			vy=velsun[1],
			vz=velsun[2],
			)
	Nplanets += 1

	planets = ['Mercury',
			    'Venus',
			    'Earth',
			    'Mars',
			    'Jupiter',
			    'Saturn',
			    'Uranus',
			    'Neptune',
			   ]

	Nplanets += len(planets)

	if binary == True:
		#for planet in planets:
		sim.add(planets, primary=sim.particles[2])
	else:
		#for planet in planets:
		sim.add(planets, primary=sim.particles[1])

	if sizedist == True:
		pos, vel, Nparts, radii = conesetup()

		sim = bss.addpartcart(sim,
							  pos,
							  vel,
							  Nparts,
							  radii,
							  )

	elif sizedist == False:
		pos, vel, Nparts = conesetup()

		sim = bss.addpartcart(sim,
							  pos,
							  vel,
							  Nparts,
							  )

	if radpress == True:
		sim, Nplanets = type1radpress()

	inttime = 0
	sim.integrator = integrator
	landed = 0
	binland = 0
	rminds = []

	return sim







def sim_loop(time, *rest):
	global Nparts
	global sim
	global pos
	global vel
	global radii
	global binx
	global biny
	global binz
	global inttime
	global Nplanets
	global landed
	global binland
	global rminds


	# Update gravitational potentials
	if rotation == True:
		sim = ae.gravity(sim, mtarg, Nplanets, Nparts, aT, bT, cT, omegaT, axrotT, axdirT, time, binary=False)
		if binary == True:
			sim = ae.gravity(sim, mbin, Nplanets, Nparts, aB, bB, cB, omegaB, axrotB, axdirB, time, binary=True)

	sim.integrate(time)

	# Remove and count particles
	if rotation == True:
		sim, Nparts, rminds, landed, landx, landy, landz, simorigin = bss.rmlandedparts(sim,
																						Nplanets,
																						Nparts,
																						aT,
																						bT,
																						cT,
																						landed,
																						rminds,
																						)
	else:
		sim, Nparts, rminds, landed, landx, landy, landz, simorigin = bss.rmlandedparts(sim,
																						Nplanets,
																						Nparts,
																						aT,
																						bT,
																						cT,
																						landed,
																						rminds,
																						)
	if binary == True:
		if rotation == True:
			sim, Nparts, rminds, binland, binx, biny, binz = bss.rmbinparts(sim,
																			Nplanets,
																			Nparts,
																			aB,
																			bB,
																			cB,
																			binland,
																			rminds,
																			)
		else:
			sim, Nparts, rminds, binland, binx, biny, binz = bss.rmbinparts(sim,
																			Nplanets,
																			Nparts,
																			aB,
																			bB,
																			cB,
																			binland,
																			rminds,
																			)

	print ('Nparts =', Nparts)
	print ('time =', time)
	print (np.sqrt((sim.particles[1].x - sim.particles[0].x) ** 2 + (sim.particles[1].y - sim.particles[0].y) ** 2 + (sim.particles[1].z - sim.particles[0].z) ** 2))



	#	print np.sqrt((np.asarray(landx)-sim.particles[0].x)**2 + (np.asarray(landy)-sim.particles[0].y)**2 + (np.asarray(landz)-sim.particles[0].z)**2)
	#	print landx, landy, landz

	# SAVE DATA
	bss.datasave(sim, 'particledata' + str(inttime) + condit + '.txt', Nplanets, Nparts)
	bss.rmdatasave(landx, landy, landz, 'rmlandpartdata' + str(inttime) + condit + '.txt')

	if binary == True:
		bss.rmdatasave(binx, biny, binz, 'rmbinpartdata' + str(inttime) + condit + '.txt')

	inttime += 1
	print (' ')



### Main Run ###
sim_setup()

omegaT = 360. / perT
omegaB = 360. / perB

inttime = 0
N_out = (1 / dt) * tmax
times = np.linspace(0, tmax, N_out)
#processes = os.getenv('SLURM_NTASKS')

#if Parallel == False:
for i, time in enumerate(times):
	sim_loop(time)
	print('times =', ctime.time(), ctime.clock())
	np.savetxt('cputimes'+condit+'.txt', np.asarray([ctime.time(), ctime.clock()]), fmt='%.18e', delimiter=' ', newline='\n')

# elif Parallel == True:
# 	new_sim_loop = functools.partial(sim_loop, times)
# 	p = Pool(int(processes))
# 	p.map(new_sim_loop, times, chunksize=10)
# 	p.join()
# 	p.close()

