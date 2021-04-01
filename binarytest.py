import rebound
import numpy as np
import time as ctime
import basicsolarsystem1015 as bss
import matplotlib.pyplot as plt
#import erroranalysis as ea
import additionaleffects as ae
import functools
import os


from bodyparams import *
from integratorparams import *



global sim
sim = rebound.Simulation()


# add particles
def conesetup(*rest):
	# Set up cone
	global Nparts
	global sim

	Nparts = int(Nparts)

	cpos = (0, 0, 0)
	cvel = (0, 0, 0)


	pos, x1, y1, z1, thetaprime, r, s, pos0p, vpos = bss.ejectaconepos(mex, aT,
																			  bT,
																			  cT,
																			  hpart,
																			  lat,
																			  lon,
																			  beta,
																			  axtiltT,
																			  axdirT,
																			  Nparts,
																			  cpos,
																			  mtarg,
																			  rtarg,
																			  vi,
																			  a,
																			  mi,
																			  rhoT,
																			  rhoi,
																			  mu,
																			  K1,
																			  Ybar=0.,
															   )





	# Produce size distribution and add particles
	if sizedist == True:

		# set up size distribution
		rsize, counts, N = ae.sizedist(Nparts,
								   Res,
								   rmin,
								   rmax,
								   p,
								   )

		radii = ae.partrad(rsize,
						   counts,
						   )


	vpart = bss.veldist(sim, r, mtarg, rtarg, mex, mu, rhoT, Cvps, Ctg, Ybar)     #sim, r, mtarg, rtarg, vi, a, mi, rhoT, rhoi, mu, K1, Cvps=0., Ctg=0.8, Ybar=0.)
	print('vel=', vpart)
	#vpart = np.full(Nparts, vinit)

	vel = bss.ejectaconevel(vpos,
							pos0p,
							vpart,
							#vpos
							#thetaprime,
							#beta

		#opang,
		#					thetaprime,
		#					lat,
		#					lon,
		#					pos,
		#					vpart,
		#					cvel,
		#					axtiltT,
		#					axdirT,
		#					dt,
		#					perT,
		#					rotation
							)

	if rotation == True:
		vxrot, vyrot, vzrot = ae.rotvel(sim, perT, lat, pos, axtiltT, axdirB)
		vel = (vel[0] + vxrot, vel[1] + vyrot, vel[2] + vzrot)

	elif rotation == False:
		vel = (vel[0], vel[1], vel[2])

	velmag = np.sqrt(vel[0]**2 + vel[1]**2 + vel[2]**2)

	#print(np.min(velmag), np.max(velmag))

	if sizedist == True:
		return pos, vel, Nparts, radii
	else:
		return pos, vel, Nparts





def type1radpress(*rest):
	global Nplanets

	# add radiation pressure for type 1
	acc = ae.solarradpress(sim, Nplanets, possun, rho)

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

	return sim, Nplanets, acc



### Run all the functions ###



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
	global perT
	global axtiltT
	global axdirT

	#sim = rebound.Simulation()
	sim.units = (timeunits, dist, mass)

	if ellipsoid == True:
		sim.add(m=mtarg, r=aT, x=0,
				y=0,
				z=0,
				vx=0,
				vy=0,
				vz=0,
				hash="target"
				)
	else:
		sim.add(m=mtarg,
				r=aT,
				x=0,
				y=0,
				z=0,
				vx=0,
				vy=0,
				vz=0,
				hash="target"
				)
	Nplanets = 1

	if binary == True:
		if ellipsoid == True:
			sim.add(r=aB,
					m=mbin,
					a=abin,
					e=ebin,
					inc=ibin,
					primary=sim.particles[0],
					hash='bin'
					)
		else:
			sim.add(r=aB,
					m=mbin,
					a=abin,
					e=ebin,
					inc=ibin,
					primary=sim.particles[0],
					hash='bin'
					)
		Nplanets += 1
		#sim.remove(1)

	# add massless sun
	#sim.add(x=possun[0],
	#		y=possun[1],
	#		z=possun[2],
	#		vx=velsun[0],
	#		vy=velsun[1],
	#		vz=velsun[2],
	#		)
	#Nplanets += 1

	planets = ['Mercury',
			    'Venus',
			    'Earth',
			    'Mars',
			    'Jupiter',
			    'Saturn',
			    'Uranus',
			    'Neptune',
			   ]

#	Nplanets += len(planets)

#	if binary == True:
		#for planet in planets:
#		sim.add(planets, primary=sim.particles[2])
#	else:
		#for planet in planets:
#		sim.add(planets, primary=sim.particles[1])

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
		sim, Nplanets, acc = type1radpress()
		print('acc =', np.sqrt(acc[0] ** 2 + acc[1] ** 2 + acc[2] ** 2))




	inttime = 0
	sim.integrator = integrator
	landed = 0
	binland = 0
	rminds = []





	# Update gravitational potentials
	if ellipsoid == True:


		#ggx, ggy, ggz = ae.gravity(sim, mtarg, Nplanets, Nparts, aT, bT, cT, perT, axrotT, axdirT, binary=False)

		#ggx, ggy, ggz = ae.tiltrot(sim, perT, Nplanets, Nparts, axrotT, axdirT, mtarg, aT, bT, cT)
		#print('gravity =', np.sqrt(ggx**2+ggy**2+ggz**2))

		sim.gravity_ignore_terms = 1

		ps = sim.particles

		nontarg = Nparts
		if binary == True:
			nontarg += 1
		if radpress == True:
			nontarg += 1

		if nontarg >= 1:
			bodyrange = np.arange(1, nontarg)
			for j in bodyrange:
				def ellipsoidgrav(reb_sim):
					''' Actual function for ellipsoidal gravity'''
					global sim
					global perT
					global perB
					global binary
					global mtarg
					global axtiltT
					global axtiltB
					global a1


					if binary == True:
						per = perB
						axdir = axdirB
						axtilt = axtiltB
						a1 = aB
						b1 = bB
						c1 = cB
					else:
						per = perT
						axdir = axdirT
						axtilt = axtiltT
						a1 = aT
						b1 = bT
						c1 = cT

					dt = sim.dt

					omega = 360. / perT

					rot = omega * dt

					sim = ae.addnetgrav(sim, mtarg, a1, b1, c1, Nplanets, Nparts, dt, omega, timestep, axdir, axtilt, binary)



				sim.additional_forces = ellipsoidgrav

		if binary == True:
			bodies = np.arange(0, Nplanets + Nparts)
			exclude = [1]
			bodyrange = np.delete(bodies, exclude)

			ax = np.asarray([sim.particles[int(j)].ax for j in bodyrange])
			ay = np.asarray([sim.particles[int(j)].ay for j in bodyrange])
			az = np.asarray([sim.particles[int(j)].az for j in bodyrange])

			ggx, ggy, ggz = ae.gravity(sim, mbin, Nplanets, Nparts, aB, bB, cB, omegaB, axtiltB, axdirB,
									   binary=True)

			#sim.additional_forces = ellipsoidgrav

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
	global timestep
	global remaining


	#orb = sim.particles[1].calculate_orbit(primary=sim.particles[0])
	#print('target =', np.sqrt(sim.particles[4].x**2 + sim.particles[4].y**2 + sim.particles[4].z**2))
	#print('target2 =', np.sqrt(orb.x ** 2 + orb.y ** 2 + orb.z ** 2))
	#sim.status()



	sim.integrate(time)
	timestep += 1


	# Remove and count particles

	# if ellipsoid == True:
	# 	sim, Nparts, rminds, landed, landx, landy, landz, simorigin = bss.rmlandedparts(sim,
	# 	 																				Nparts,
	# 	 																				aT,
	# 	 																				bT,
	# 	 																				cT,
	# 	 																				landed,
	# 	 																				rminds,
	# 																					binary,
	# 																					rbin,
	# 																					axdirT,
	# 																					axtiltT,
	# 																					360. / perT,
	# 	 																				)
	#else:

	#sim, landed, Nparts, rminds = bss.rmparticles(sim, Nparts, landed, inttime, aT, aB, condit, rminds)

	sim, landed, Nparts = bss.rmland(sim,
 	 								 Nparts,
 	 								 aT,
 	 								 bT,
 	 								 cT,
 	 								 landed,
 	 								 inttime,
 	 								 condit,
 	 								 axdirT,
 	 								 axtiltT,
 	 								 perT,
 	 								 timestep,
 	 								 aB,
 	 								 rotation,
 	 						   )

		# sim, Nparts, rminds, landed, landx, landy, landz, simorigin = bss.rmlandedparts(sim,
	 	# 																				Nparts,
	 	# 																				aT,
	 	# 																				bT,
	 	# 																				cT,
	 	# 																				landed,
	 	# 																				rminds,
		# 																				binary,
		# 																				rbin,
	 	# 																				)

		#print(np.sqrt(np.asarray(landx) ** 2 + np.asarray(landy) ** 2 + np.asarray(landz) ** 2))

	# if binary == True:
	# 	sim, Nparts, rminds, binland, binx, biny, binz = bss.rmbinparts(sim,
	# 																	rbin,
	# 																	Nparts,
	# 																	binland,
	# 																	rminds,
	# 																	)

	print ('Nparts =', Nparts)

	print ('time =', time)



	#p = sim.particles[6]
	#print('particle 5 a =', np.sqrt(p.ax**2 + p.ay**2 + p.az**2))
	#print('particle 5 v =', np.sqrt(p.vx**2 + p.vy**2 + p.vz**2))
	#print('particle 5 r =', (np.sqrt(p.x ** 2 + p.y ** 2 + p.z ** 2))-rtarg)

	#	print np.sqrt((np.asarray(landx)-sim.particles[0].x)**2 + (np.asarray(landy)-sim.particles[0].y)**2 + (np.asarray(landz)-sim.particles[0].z)**2)
	#	print landx, landy, landz

	# SAVE DATA
	bss.datasave(sim, 'particledata' + str(inttime) + condit + '.txt', Nplanets, Nparts)


	#if binary == True:
	#	bss.rmdatasave(binx, biny, binz, 'rmbinpartdata' + str(inttime) + condit + '.txt')

	inttime += 1
	print (' ')



### Main Run ###

timestep = 1
sim = sim_setup()

omegaT = 360. / perT
omegaB = 360. / perB

inttime = 0
N_out = (1 / dt) * tmax
times = np.linspace(0, tmax, N_out)
#sim.move_to_com()
landed = 0

rminds = []


#print('aellip =', np.sqrt(p.ax**2 + p.ay**2 + p.az**2))
for i, time in enumerate(times):
	sim_loop(time)
	#print('mtarg =', mtarg)
	#print('targetbody =', sim.particles["target"])
#	print('times =', ctime.time(), ctime.clock())
	np.savetxt('cputimes'+condit+'.txt', np.asarray([ctime.time(), ctime.clock()]), fmt='%.18e', delimiter=' ', newline='\n')
