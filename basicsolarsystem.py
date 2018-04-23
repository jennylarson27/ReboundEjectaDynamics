### Basic Solar System Model ###
import rebound
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import csv
from mpl_toolkits.mplot3d import Axes3D
#from sympy.functions.special.delta_functions import Heaviside


def createsim(time, dist, mass):
	# Create rebound simulation
	sim = rebound.Simulation()
	sim.units = (time, dist, mass)

	return sim


def timecc(projdiam, velimpactor):
	tcc = projdiam / velimpactor

	return tcc


def addplanets(sim, planets, center=0):
	# Enter body parameters
	Nplanets = len(planets)
	for body in planets:
		sim.add(body, primary=sim.particles[center])

	return sim, Nplanets


def addplanetorbs(sim, orbsfile, Nplanets, center=0):
	with open(orbsfile) as csvfile:
		read = csv.reader(csvfile, delimiter=',')
		for row in read:
			sim.add(m=float(row[0]),
					r=float(row[1]),
					a=float(row[2]),
					e=float(row[3]),
					inc=float(row[4]),
					Omega=float(row[5]),
					omega=float(row[6]),
					f=float(row[7]),
					primary=sim.particles[center])
			Nplanets += 1

	return sim, Nplanets


def arand(a0, Nparts, da):
	a = 2 * da * np.random.random_sample(Nparts) + (a0 - da)
	return a


def erand(e0, Nparts, de):
	e = 2 * de * np.random.random_sample(Nparts) + (e0 - de)
	return e


def irand(i0, Nparts, di):
	i = 2 * di * np.random.random_sample(Nparts) + (i0 - di)
	return i


def frand(f0, Nparts, df):
	f = 2 * df * np.random.random_sample(Nparts) + (f0 - df)
	return f


def addparticles(sim, a0, e0, inc0, periap0, ascnode0, delta, Nparts):
	# Make sure Nparts is an int
	Nparts = int(Nparts)

	# Set variations per particle
	a = np.random.choice(np.arange(a0 - delta,
								   a0 + delta,
								   delta / 10.), Nparts)
	e = np.random.choice(np.arange(e0 - delta,
								   e0 + delta,
								   delta / 10.), Nparts)
	inc = np.random.choice(np.arange(inc0 - delta,
									 inc0 + delta,
									 delta / 10.), Nparts)
	ascnode = np.random.choice(np.arange(ascnode0 - delta,
										 ascnode0 + delta,
										 delta / 10.), Nparts)
	periap = np.random.choice(np.arange(periap0 - delta,
										periap0 + delta,
										delta / 10.), Nparts)

	for partnum in range(0, Nparts, 1):
		part = sim.add(a=a[partnum],
					   e=e[partnum],
					   inc=inc[partnum],
					   Omega=ascnode[partnum],
					   omega=periap[partnum],
					   hash=partnum)

	return sim, Nparts


def ejectaconepos(radius, h, theta, phi, beta, Nparts, tpos):
	'''Setting ejecta cone position'''

	r0 = radius + h

	# cone base
	x0 = radius * np.sin(np.radians(phi)) * np.cos(np.radians(theta))
	y0 = radius * np.sin(np.radians(phi)) * np.sin(np.radians(theta))
	z0 = radius * np.cos(np.radians(phi))

	# disk center
	x0p = r0 * np.sin(np.radians(phi)) * np.cos(np.radians(theta))
	y0p = r0 * np.sin(np.radians(phi)) * np.sin(np.radians(theta))
	z0p = r0 * np.cos(np.radians(phi))

	# disk radii
	rdisk = h * np.tan(np.radians(beta / 2.))

	r = rdisk * np.random.random_sample(Nparts)

	thetaprime = 360. * np.random.random_sample(Nparts)

	opang = np.degrees(np.arctan(r / h))

	mag1 = np.sqrt(h**2 + r**2)#np.sqrt(x1**2 + y1**2 + z1new**2)

	phiang  = opang * np.cos(np.radians(thetaprime))
	thetang = opang * np.sin(np.radians(thetaprime))

	x = mag1 * np.sin(np.radians(phi + phiang)) * np.cos(np.radians(theta + thetang)) + x0
	y = mag1 * np.sin(np.radians(phi + phiang)) * np.sin(np.radians(theta + thetang)) + y0
	z = mag1 * np.cos(np.radians(phi + phiang)) + z0


	# calculate unit vector
	xu = x - x0
	yu = y - y0
	zu = z - z0

	unitv = (xu / np.sqrt(xu**2 + yu**2 + zu**2), yu / np.sqrt(xu**2 + yu**2 + zu**2), zu / np.sqrt(xu**2 + yu**2 + zu**2))

	pos = (x + tpos[0], y + tpos[1], z + tpos[2])

	return pos, x, y, z, unitv, r


def ejectaconevel(unitv, v0, tvel):
	'''Setting the ejecta cone velocity'''

	vx = v0 * unitv[0] #* np.sin(np.radians(opang)) * np.cos(np.radians(thetaprime))
	vy = v0 * unitv[1] #* np.sin(np.radians(opang)) * np.sin(np.radians(thetaprime))
	vz = v0 * unitv[2] #* np.cos(np.radians(opang))

	vel = (vx + tvel[0], vy + tvel[1], vz + tvel[2])

	return vel





def addpartcart(sim, pos, vel, Nparts, radmass=0):
	'''Adding particles'''
	x = pos[0]
	y = pos[1]
	z = pos[2]

	vx = vel[0]
	vy = vel[1]
	vz = vel[2]

	# add particles
	if radmass.any() != 0:
		rad = radmass

		for partnum in range(Nparts):
			part = sim.add(r=rad[partnum],
						   x=x[partnum],
						   y=y[partnum],
						   z=z[partnum],
						   vx=vx[partnum],
						   vy=vy[partnum],
						   vz=vz[partnum])

	else:
		for partnum in range(Nparts):
			part = sim.add(x=x[partnum],
						   y=y[partnum],
						   z=z[partnum],
						   vx=vx[partnum],
						   vy=vy[partnum],
						   vz=vz[partnum])

	return sim


def helioorbits(sim):
	# Output orbits in Heliocentric coordinates
	bodies = sim.calculate_orbits(heliocentric=True)

	return bodies, sim


def baryorbits(sim):
	# Output orbits in Barycentric coordinates
	bodies = sim.calculate_orbits(barycentric=True)

	return bodies, sim


def motion(sim, tmax, dt):
	# sim.move_to_com()           # move to center of mass
	sim.dt = dt
	N_out = (1 / dt) * tmax
	# number of times position is taken during each interval

	times = np.linspace(0, tmax, N_out)

	return sim, times


def setintegrator(sim, integrator, times, dt, outlim, inlim, Nplanets):
	sim.integrator = integrator
	partfile = 'partorbits'
	for i, time in enumerate(times):
		print(time)
		filenm = partfile + str(int(time)) + '.txt'
		sim.integrate(time)  # sim.t+dt)
		sim, Nparts = rmparticles(sim, outlim, inlim, Nplanets)
		f = textorbits(sim, Nplanets, filenm, addbodies=None)
		partplot(sim, Nplanets, time)

	return sim


def rmparticles(sim, Nplanets, rad, outer, landed=0, gone=0, binland=0, binr=0):
	rminds = []
	rmx = []
	rmy = []
	rmz = []
	partind = Nplanets
	for p in sim.particles:
		if np.sqrt(p.x ** 2 + p.y ** 2 + p.z ** 2) <= rad:
			sim.remove(partind)
			rmx.append(p.x)
			rmy.append(p.y)
			rmz.append(p.z)
			rminds.append(partind)
			landed += 1
		if np.sqrt(p.x ** 2 + p.y ** 2 + p.z ** 2) > outer:
			sim.remove(partind)
			rminds.append(partind)
			gone += 1
		partind += 1

	if binr != 0:
		partind = Nplanets
		pbin = sim.particles[2]
		for p in sim.particles:
			if np.sqrt((p.x-pbin.x)**2 + (p.y-pbin.y)**2 + (p.z-pbin.z)**2) <= binr:
				sim.remove(partind)
				rminds.append(partind)
				binland += 1
			partind += 1

	return sim, rminds, landed, gone, binland, rmx, rmy, rmz


def updatepos(pos0, rminds):
	x0 = np.delete(pos0[0], rminds)
	y0 = np.delete(pos0[1], rminds)
	z0 = np.delete(pos0[2], rminds)

	pos0 = (x0, y0, z0)

	return pos0


def updatevel(vel0, rminds):
	vx0 = np.delete(vel0[0], rminds)
	vy0 = np.delete(vel0[1], rminds)
	vz0 = np.delete(vel0[2], rminds)

	vel0 = (vx0, vy0, vz0)

	return vel0


def calcorbits(sim, Nplanets):
	parts = sim.calculate_orbits()[Nplanets - 1:]
	aparts = np.array([o.a for o in parts])
	eparts = np.array([o.e for o in parts])
	iparts = np.array([o.inc for o in parts])

	planets = sim.calculate_orbits()[:Nplanets - 1]
	aplanets = np.array([o.a for o in planets])
	eplanets = np.array([o.e for o in planets])
	iplanets = np.array([o.inc for o in planets])

	return sim, aparts, eparts, iparts, aplanets, eplanets, iplanets


def textorbits(sim, Nplanets, partfile, addbodies=None):
	f = open(partfile, 'w')
	if addbodies == None:
		parts = sim.particles_ascii()  # sim.calculate_orbits()[Nplanets-1:]
	else:
		parts = sim.particles_ascii()  # sim.calculate_orbits()[Nplanets-1+addbodies:]
	f.write(parts)
	# for p in parts:
	#	f.write(p)
	f.close()
	return f


def init_posvel(sim, time):
	x = []
	y = []
	z = []

	vx = []
	vy = []
	vz = []

	if time == 0:
		for p in sim.particles:
			x.append(p.x)
			y.append(p.y)
			z.append(p.z)

			vx.append(p.vx)
			vy.append(p.vy)
			vz.append(p.vz)
	x = np.asarray(x)
	y = np.asarray(y)
	z = np.asarray(z)

	vx = np.asarray(vx)
	vy = np.asarray(vy)
	vz = np.asarray(vz)

	pos0 = (x, y, z)
	vel0 = (vx, vy, vz)

	return pos0, vel0


def veldispersion(sim, fname, Nplanets, pos0, vel0, rad):
	'''Plots delta V vs delta R

    Parameters
    ----------
    sim     : Rebound Simulation
    fname   : str; name of file
    Nplanets: int; number of non-particle objects in SS
    pos0    : tuple; initial particle positions (x, y, z)
    vel0    : tuple; initial particle velocities (vx, vy, vz)

    Returns
    -------
    '''

	x0 = pos0[0]
	y0 = pos0[1]
	z0 = pos0[2]

	vx0 = vel0[0]
	vy0 = vel0[1]
	vz0 = vel0[2]

	x = []
	y = []
	z = []

	vx = []
	vy = []
	vz = []

	for p in sim.particles:
		x.append(p.x)
		y.append(p.y)
		z.append(p.z)

		vx.append(p.vx)
		vy.append(p.vy)
		vz.append(p.vz)

	x = np.asarray(x)
	y = np.asarray(y)
	z = np.asarray(z)
	vx = np.asarray(vx)
	vy = np.asarray(vy)
	vz = np.asarray(vz)

	delx = np.asarray(x[Nplanets:]) - x0[Nplanets:]
	dely = np.asarray(y[Nplanets:]) - y0[Nplanets:]
	delz = np.asarray(z[Nplanets:]) - z0[Nplanets:]

	delvx = np.asarray(vx[Nplanets:]) - vx0[Nplanets:]
	delvy = np.asarray(vy[Nplanets:]) - vy0[Nplanets:]
	delvz = np.asarray(vz[Nplanets:]) - vz0[Nplanets:]

	# separ  = np.sqrt(delx**2 + dely**2 + delz**2)
	disper = np.sqrt(delvx ** 2 + delvy ** 2 + delvz ** 2) / 86400.

	dist0 = np.sqrt(x0[Nplanets:] ** 2 + y0[Nplanets:] ** 2 + z0[Nplanets:] ** 2)
	veloc0 = np.sqrt(vx0[Nplanets:] ** 2 + vy0[Nplanets:] ** 2 + vz0[Nplanets:] ** 2) / 86400.

	dist = np.sqrt(x[Nplanets:] ** 2 + y[Nplanets:] ** 2 + z[Nplanets:] ** 2)
	veloc = np.sqrt(vx[Nplanets:] ** 2 + vy[Nplanets:] ** 2 + vz[Nplanets:] ** 2) / 86400.

	separ = dist - rad
	# separation = np.where((dist - dist0) < 0, separ * -1)
	# dispersion = np.where((veloc - veloc0) < 0, disper * -1)
	for i in range(dist.shape[0]):
		if dist[i] - dist0[i] < 0:
			separ[i] = -separ[i]
		if veloc[i] - veloc0[i] < 0:
			disper[i] = -disper[i]

	# print separ, disper, veloc

	fig = plt.figure()
	plt.scatter(separ, disper)
	plt.xlabel('Height (km)')
	plt.ylabel(r'$\Delta$ v (km/s)')
	plt.xlim(0, 400)
	plt.ylim(-.2, 0)
	plt.savefig(fname)


def orbplot(a, e, i, ascnode, periap):
	'''Calculates cartesian coordinates for planetary orbits

    Parameters
    ----------
    a : float; semimajor axis
    e : float; eccentricity
    i : float; inclination
    ascnode : float; Omega, longitude of ascending node
    periap  : float; omega, argument of periapse

    Returns
    -------
    x : array; x-coordinate values for full orbit
    y : array; y-coordinate values for full orbit
    z : array; z-coordinate values for full orbit

    Revised
    -------
    19-06-2017 created function
    '''

	theta = np.arange(0, 360., 0.01)

	r = a * ((1 - e ** 2) / (1 + e * np.cos(np.radians(theta))))

	x = r * (np.cos(np.radians(ascnode)) * np.cos(np.radians(periap + theta)) -
			 np.sin(np.radians(ascnode)) * np.sin(np.radians(periap + theta)) * np.cos(np.radians(i)))
	y = r * (np.sin(np.radians(ascnode)) * np.cos(np.radians(periap + theta)) +
			 np.cos(np.radians(ascnode)) * np.sin(np.radians(periap + theta)) * np.cos(np.radians(i)))
	z = r * (np.sin(np.radians(periap + theta)) * np.sin(np.radians(i)))

	return x, y, z


def sunupdate(sim, xplanet, yplanet, zplanet, binary=False):
	'''Updates planet positions relative to sun for plotting

	Parameters
	----------
	sim     : REBOUND simulation
	xplanet : array; x coordinates of planetary orbits
	yplanet : array; y coordinates of planetary orbits
	zplanet : array; z coordinates of planetary orbits
	binary  :  bool; (default=False) determines if binary component exists

	Returns
	-------
	xp : array; updated x positions
	yp : array; updated y positions
	zp : array; updated z positions

	Revised
	-------
	20-06-2017 created function
	'''

	xs = sim.particles[1].x
	ys = sim.particles[1].y
	zs = sim.particles[1].z

	xp = xs + xplanet
	yp = ys + yplanet
	zp = zs + zplanet

	if binary != False:
		xp[0] = xplanet[0] - xs
		yp[0] = yplanet[0] - ys
		zp[0] = zplanet[0] - zs

	return xp, yp, zp


def plot3d(sim, fname, Nplanets, rad=None, binary=None, binorb=(0,0,0), size=10):
	x = []
	y = []
	z = []

	for p in sim.particles:
		x.append(p.x)
		y.append(p.y)
		z.append(p.z)

	x = np.asarray(x[Nplanets:])
	y = np.asarray(y[Nplanets:])
	z = np.asarray(z[Nplanets:])

	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	# ax = plt.axes(projection='3d')
	ax.scatter(x, y, z)
	if rad != None:
		ax.scatter(0, 0, 0, s=np.pi * rad ** 2, c='red', alpha=0.75)
		if binary != None:
			xbin = binorb[0]
			ybin = binorb[1]
			zbin = binorb[2]
			ax.plot(xbin, ybin, zbin)
			bin = sim.particles[2]
			bx = bin.x
			by = bin.y
			bz = bin.z
			ax.scatter(bx, by, bz, s=np.pi * binary ** 2, c='magenta', alpha=0.75)
	ax.set_xlim(-size, size)
	ax.set_ylim(-size, size)
	ax.set_zlim(-size, size)
	ax.set_xlabel('x (km)')
	ax.set_ylabel('y (km)')
	ax.set_zlabel('z (km)')
	plt.savefig(fname)


def plotxyplane(sim, fname, Nplanets, xplanet=None, yplanet=None, rad=None, binary=None, binorb=(0,0,0), hill=None, size=10):
	'''Plots positions on xy plane

	Parameters
	----------
	sim      : REBOUND simulation
	fname    :   str; save as name of plot
	Nplanets :   int; number of planets
	xplanet  : array; x coordinate of planet orbits
	yplanet  : array; y coordinate of planet orbits
	rad      : float; (default=None) radius of central body
	binary   : float; (default=None) radius of binary
	hill     : float; (default=None) radius of hill sphere boundary
	size     : float; (default=10) size of positive x/y scale

	Returns
	-------
	fig : figure in xy plane
	'''

	fig, ax = plt.subplots()
	if hill != None:
		circle = plt.Circle((0, 0), radius=hill, linestyle='dashed', facecolor='None', edgecolor='black')
		ax.add_artist(circle)
	x = []
	y = []

	for p in sim.particles:
		x.append(p.x)
		y.append(p.y)

	x = np.asarray(x[Nplanets:])
	y = np.asarray(y[Nplanets:])

	ax.scatter(x, y)

	if rad != None:
		ax.scatter(0, 0, s=np.pi * rad ** 2, c='red', alpha=0.75)

	if binary != None:
		ax.plot(binorb[0], binorb[1])
		bin = sim.particles[2]
		bx = bin.x
		by = bin.y
		ax.scatter(bx, by, s=np.pi * binary ** 2, c='magenta', alpha=0.75)


	if xplanet != None and yplanet != None:
		ax.plot(xplanet, yplanet)

	ax.set_xlim(-size, size)
	ax.set_ylim(-size, size)
	ax.set_xlabel('x (km)')
	ax.set_ylabel('y (km)')
	plt.savefig(fname)


def plotxzplane(sim, fname, Nplanets, xplanet=None, zplanet=None, rad=None, binary=None, binorb=(0,0,0), hill=None, size=10):
	'''Plots positions on xz plane

		Parameters
		----------
		sim      : REBOUND simulation
		fname    :   str; save as name of plot
		Nplanets :   int; number of planets
		xplanet  : array; x coordinates for planet orbits
		zplanet  : array; z coordinates for planet orbits
		rad      : float; (default=None) radius of central body
		binary   : float; (default=None) radius of binary
		hill     : float; (default=None) radius of hill sphere boundary
		size     : float; (default=10) size of positive x/z scale

		Returns
		-------
		fig : figure in xz plane
		'''

	fig, ax = plt.subplots()
	if hill != None:
		circle = plt.Circle((0, 0), radius=hill, linestyle='dashed', facecolor='None', edgecolor='black')
		ax.add_artist(circle)

	x = []
	z = []

	for p in sim.particles:
		x.append(p.x)
		z.append(p.z)

	x = np.asarray(x[Nplanets:])
	z = np.asarray(z[Nplanets:])

	ax.scatter(x, z)
	if rad != None:
		ax.scatter(0, 0, s=np.pi * rad ** 2, c='red', alpha=0.75)

	if binary != None:
		ax.plot(binorb[0], binorb[2])
		bin = sim.particles[2]
		bx = bin.x
		bz = bin.z
		ax.scatter(bx, bz, s=np.pi * binary ** 2, c='magenta', alpha=0.75)


	if xplanet != None and zplanet != None:
		ax.plot(xplanet, zplanet)

	ax.set_xlim(-size, size)
	ax.set_ylim(-size, size)
	ax.set_xlabel('x (km)')
	ax.set_ylabel('z (km)')
	plt.savefig(fname)



def plotxyxzplanes(sim, fname, Nplanets, xplanet=None, yplanet=None, zplanet=None, rad=None, prad=1, binary=None, binorb=(0,0,0), hill=None, size=10):

	fig, ax = plt.subplots(2)
	if hill != None:
		circle = plt.Circle((0, 0), radius=hill, linestyle='dashed', facecolor='None', edgecolor='black')
		ax[0].add_artist(circle)
		ax[1].add_artist(circle)
	x = []
	y = []
	z = []

	for p in sim.particles:
		x.append(p.x)
		y.append(p.y)
		z.append(p.z)

	x = np.asarray(x[Nplanets:])
	y = np.asarray(y[Nplanets:])
	z = np.asarray(z[Nplanets:])

	ax[0].scatter(x, y)
	ax[1].scatter(x, z)

	if rad != None:
		ax[0].scatter(0, 0, s=np.pi * rad ** 2, c='red', alpha=0.75)
		ax[1].scatter(0, 0, s=np.pi * rad ** 2, c='red', alpha=0.75)

	if binary != None:
		ax[0].plot(binorb[0], binorb[1])
		ax[1].plot(binorb[0], binorb[2])

		bin = sim.particles[2]
		bx = bin.x
		by = bin.y
		bz = bin.z

		ax[0].scatter(bx, by, s=np.pi * binary ** 2, c='black', alpha=0.75)
		ax[1].scatter(bx, bz, s=np.pi * binary ** 2, c='black', alpha=0.75)


	if xplanet != None and yplanet != None:
		xplanet = x
		yplanet = y
		ax[0].plot(xplanet, yplanet, s=np.pi * prad**2, c='black')

	if xplanet != None and zplanet != None:
		xplanet = x
		zplanet = z
		ax[1].plot(xplanet, zplanet, s=np.pi * prad**2, c='black')

	ax[0].set_xlim(-0, size)
	ax[0].set_ylim(-size, 5)
	ax[1].set_xlim(-0, size)
	ax[1].set_ylim(-10, 10)

	ax[0].set_xlabel('x (km)')
	ax[0].set_ylabel('y (km)')
	ax[1].set_xlabel('x (km)')
	ax[1].set_ylabel('z (km)')

	plt.savefig(fname)


def datasave(sim, fname, binary=False):
	x = []
	y = []
	z = []
	vx = []
	vy = []
	vz = []

	for p in sim.particles:
		x.append(p.x)
		y.append(p.y)
		z.append(p.z)
		vx.append(p.vx)
		vy.append(p.vy)
		vz.append(p.vz)
	x = np.asarray(x[1:])
	y = np.asarray(y[1:])
	z = np.asarray(z[1:])
	vx = np.asarray(vx[1:])
	vy = np.asarray(vy[1:])
	vz = np.asarray(vz[1:])

	data = np.zeros((6, len(x)))

#	names = []
#	names.append('Sun')
#	if binary == True:
#		names.append('Binary')
#	for plan in planets:
#		names.append(plan)
#	for n in range(Nparts):
#		names.append('part' + str(int(n)))

	#data[0, :] = names
	data[0, :] = x
	data[1, :] = y
	data[2, :] = z
	data[3, :] = vx
	data[4, :] = vy
	data[5, :] = vz

	np.savetxt(fname, data, fmt='%.18e', delimiter=' ', newline='\n')




def rmdatasave(rmx, rmy, rmz, fname):
	data = np.zeros((3, len(rmx)))

	data[0, :] = rmx
	data[1, :] = rmy
	data[2, :] = rmz

	np.savetxt(fname, data, fmt='%.18e', delimiter=' ', newline='\n')









def posplot(sim, time, fname, crop=(-1e-1, 1e-1, -1e-1, 1e-1), centered=0):
	'''centered=False defaults to the sun as the center; centered=particle index centers plot on a certain particle'''

	x = []
	y = []
	z = []
	if centered == 0:
		for p in sim.particles:
			x.append(p.x)
			y.append(p.y)
			z.append(p.z)

		fig, ax = plt.subplots(2, 2, sharex=True, sharey=True)
		plt.subplot(2, 2, 3)
		plt.scatter(x, y)
		plt.scatter(x[centered], y[centered], color='red')
		plt.xlabel('x [AU]')
		plt.ylabel('y [AU]')

		plt.subplot(2, 2, 1)
		plt.scatter(x, z)
		plt.scatter(x[centered], z[centered], color='red')
		plt.xlabel('x [AU]')
		plt.ylabel('z [AU]')

		plt.subplot(2, 2, 4)
		plt.scatter(z, y)
		plt.scatter(z[centered], y[centered], color='red')
		plt.xlabel('y [AU]')
		plt.ylabel('z [AU]')

	else:
		xmin = crop[0]
		xmax = crop[1]
		ymin = crop[2]
		ymax = crop[3]
		cx = sim.particles[centered].x
		cy = sim.particles[centered].y
		cz = sim.particles[centered].z
		for p in sim.particles:
			x.append(p.x - cx)
			y.append(p.y - cy)
			z.append(p.z - cz)

		fig, ax = plt.subplots(2, 2, sharex=True, sharey=True)
		plt.subplot(2, 2, 3)
		plt.scatter(x[centered + 1:], y[centered + 1:])
		plt.scatter(x[centered], y[centered], color='purple')
		plt.scatter(x[:centered], y[:centered], color='red')
		plt.xlabel('x [AU]')
		plt.ylabel('y [AU]')
		#		plt.xlim(xmin, xmax)
		#		plt.ylim(ymin, ymax)

		plt.subplot(2, 2, 1)
		plt.scatter(x[centered + 1:], z[centered + 1:])
		plt.scatter(x[centered], z[centered], color='purple')
		plt.scatter(x[:centered], z[:centered], color='red')
		plt.xlabel('x [AU]')
		plt.ylabel('z [AU]')
		#		plt.xlim(xmin, xmax)
		#		plt.ylim(ymin, ymax)

		plt.subplot(2, 2, 4)
		plt.scatter(z[centered + 1:], y[centered + 1:])
		plt.scatter(z[centered], y[centered], color='purple')
		plt.scatter(z[:centered], y[:centered], color='red')
		plt.xlabel('y [AU]')
		plt.ylabel('z [AU]')
	#		plt.xlim(xmin, xmax)
	#		plt.ylim(ymin, ymax)

	plt.savefig(fname)


def partplot(sim, Nplanets, time):
	cx = sim.particles[Nplanets].x
	cy = sim.particles[Nplanets].y
	cz = sim.particles[Nplanets].z
	cr = sim.particles[Nplanets].r

	# x-y
	fig = plt.figure()
	ax1, ax2, ax3 = plt.subplot(3)
	ax.set_xlabel('x [AU]')
	ax.set_ylabel('y [AU]')

	targ = patches.Circle((0, 0), cr, facecolor='darkgray', edgecolor='black')
	ax.add_patch(targ)

	for i, p in enumerate(sim.particles[Nplanets:]):
		px = p.x - cx
		py = p.y - cy
		ax.plot(px, py, marker='.', color='blue')

	plt.savefig('partposx-y' + str(int(time)) + '.png')

	# x-z
	fig = plt.figure()
	ax = plt.subplot(111)
	ax.set_xlabel('x [AU]')
	ax.set_ylabel('z [AU]')

	targ = patches.Circle((0, 0), cr, facecolor='darkgray', edgecolor='black')
	ax.add_patch(targ)

	for i, p in enumerate(sim.particles[Nplanets:]):
		px = p.x - cx
		pz = p.z - cz
		ax.plot(px, pz, marker='.', color='blue')

	plt.savefig('partposx-z' + str(int(time)) + '.png')

	# y-z
	fig = plt.figure()
	ax = plt.subplot(111)
	ax.set_xlabel('y [AU]')
	ax.set_ylabel('z [AU]')

	targ = patches.Circle((0, 0), cr, facecolor='darkgray', edgecolor='black')
	ax.add_patch(targ)

	for i, p in enumerate(sim.particles[Nplanets:]):
		py = p.y - cy
		pz = p.z - cz
		ax.plot(py, pz, marker='.', color='blue')

	plt.savefig('partposy-z' + str(int(time)) + '.png')


def ploti(aparts, iparts, Nparts, tmax):
	plt.ion()
	plt.clf()

	plt.scatter(aparts, iparts, linewidth='0')
	plt.xlabel('Semimajor Axis (AU)')
	plt.ylabel('Inclination')
	# plt.title('Inclination Versus Semimajor Axis for '+str(Nparts)+' Particles Over '+str(tmax)+' yrs')
	plt.savefig('ivsa_' + str(int(Nparts)) + 'part' + str(int(tmax)) + 'yr.png')


def plote(aparts, eparts, Nparts, tmax):
	plt.ion()
	plt.clf()

	plt.scatter(aparts, eparts, linewidth='0')
	plt.xlabel('Semimajor Axis (AU)')
	plt.ylabel('Eccentricity')
	# plt.title('Eccentricity Versus Semimajor Axis for '+str(Nparts)+' Particles Over '+str(tmax)+' yrs')
	plt.savefig('evsa_' + str(int(Nparts)) + 'part' + str(int(tmax)) + 'yr.png')


def heavisidefunc(x):
	'''
	Heaviside function

	Parameters
	----------
	x : array; values to apply heaviside function to

	Returns
	-------
	heav : array; heaviside function applied to x values
	'''

	heav = 0.5 * (np.sign(x) + 1)

	return heav


def corrD(r, pos, Nparts):
	'''
	Calculates correlation estimation C(r)

	Parameters
	----------
	r      : array; r values
	pos    : tuple; cartesian coordinates of each particle
	Nparts :   int; number of particles in set

	Returns
	-------
	corr : array; correlation estimations for each r value
	'''

	x = pos[0]
	y = pos[1]
	z = pos[2]

	for i in range(len(x)):
		for j in range(len(y)):
			if i < j:
				dist = np.sqrt((x[i] - x[j]) ** 2 + (y[i] - y[j]) ** 2 + (z[i] - z[j]) ** 2)

	csum = np.zeros(len(r))
	for k in range(len(r)):
		csum[k] = np.sum(heavisidefunc(r[k] - dist))

	corr = (2. / (Nparts * (Nparts - 1))) * csum

	return corr
