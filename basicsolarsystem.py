### Basic Solar System Model ###
import rebound
import numpy as np



def addpartcart(sim, pos, vel, Nparts, radmass=np.zeros(3)):
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

		for p in np.arange(0, Nparts, 1):
			p = int(p)
			sim.add(r=rad[p],
					x=x[p],#np.asarray([xval for xval in pos[0]]),
					y=y[p],#np.asarray([yval for yval in pos[1]]),
					z=z[p],#np.asarray([zval for zval in pos[2]]),
					vx=vx[p],#np.asarray([vxval for vxval in vel[0]]),
					vy=vy[p],#np.asarray([vyval for vyval in vel[1]]),
					vz=vz[p],#np.asarray([vzval for vzval in vel[2]]),
					)

	else:
		for p in np.arange(0, Nparts, 1):
			p = int(p)
			sim.add(x=x[p],#np.asarray([xval for xval in pos[0]]),
					y=y[p],#np.asarray([yval for yval in pos[1]]),
					z=z[p],#np.asarray([zval for zval in pos[2]]),
					vx=vx[p],#np.asarray([vxval for vxval in vel[0]]),
					vy=vy[p],#np.asarray([vyval for vyval in vel[1]]),
					vz=vz[p],#np.asarray([vzval for vzval in vel[2]]),
					)

	return sim




def datasave(sim, fname, Nplanets, Nparts):

	x = np.asarray([sim.particles[j].x for j in range(0, Nplanets + Nparts)])
	y = np.asarray([sim.particles[j].y for j in range(0, Nplanets + Nparts)])
	z = np.asarray([sim.particles[j].z for j in range(0, Nplanets + Nparts)])

	vx = np.asarray([sim.particles[j].vx for j in range(0, Nplanets + Nparts)])
	vy = np.asarray([sim.particles[j].vy for j in range(0, Nplanets + Nparts)])
	vz = np.asarray([sim.particles[j].vz for j in range(0, Nplanets + Nparts)])

	data = np.zeros((6, len(x)))


	data[0, :] = x
	data[1, :] = y
	data[2, :] = z
	data[3, :] = vx
	data[4, :] = vy
	data[5, :] = vz

	np.savetxt(fname, data, fmt='%.18e', delimiter=' ', newline='\n')






def ejectaconepos(radius, h, lat, lon, beta, Nparts, tpos):
	'''Setting ejecta cone position'''

	r0 = radius + h

	# cone base
	x0 = radius * np.cos(np.radians(lat)) * np.cos(np.radians(lon))
	y0 = radius * np.cos(np.radians(lat)) * np.sin(np.radians(lon))
	z0 = radius * np.sin(np.radians(lat))


	# disk radius
	rdisk  = h * np.tan(np.radians(beta))
	rinner = h * np.tan(np.radians((beta) - 10.))

	r = (rdisk - rinner) * np.random.random_sample(Nparts) + rinner

	thetaprime = 360. * np.random.random_sample(Nparts)

	opang = np.degrees(np.arctan(r / h))

	mag1 = np.sqrt(h**2 + r**2)#np.sqrt(x1**2 + y1**2 + z1new**2)

	phiang  = opang * np.cos(np.radians(thetaprime))
	thetang = opang * np.sin(np.radians(thetaprime))

	x = mag1 * np.cos(np.radians(lat + phiang)) * np.cos(np.radians(lon + thetang)) + x0
	y = mag1 * np.cos(np.radians(lat + phiang)) * np.sin(np.radians(lon + thetang)) + y0
	z = mag1 * np.sin(np.radians(lat + phiang)) + z0


	pos = (x + tpos[0], y + tpos[1], z + tpos[2])

	s = np.sqrt((x-x0)**2 + (y-y0)**2 + (z-z0)**2)   # distance from particles to cone base

	return pos, x, y, z, phiang, thetang, r, s




def ejectaconevel(phiang, thetang, lat, lon, v0, tvel):
	'''Setting the ejecta cone velocity'''

	#vx = []
	#vy = []
	#vz = []

	#for i in range(Nparts):
	vx = (v0 * np.cos(np.radians(lat + phiang)) * np.cos(np.radians(lon + thetang))) #* np.sin(np.radians(opang)) * np.cos(np.radians(thetaprime))
	vy = (v0 * np.cos(np.radians(lat + phiang)) * np.sin(np.radians(lon + thetang))) #* np.sin(np.radians(opang)) * np.sin(np.radians(thetaprime))
	vz = (v0 * np.sin(np.radians(lat + phiang))) #* np.cos(np.radians(opang))

	vel = (np.asarray(vx) + tvel[0], np.asarray(vy) + tvel[1], np.asarray(vz) + tvel[2])

	return vel






def rmbinparts(sim, Nplanets, Nparts, abin, bbin, cbin, binland, rminds):
	binx = []
	biny = []
	binz = []

	bodies = np.arange(0, Nplanets + Nparts)
	exclude = [1]
	bodyrange = np.delete(bodies, exclude)

	xp = np.asarray([sim.particles[int(j)].x for j in bodyrange])
	yp = np.asarray([sim.particles[int(j)].y for j in bodyrange])
	zp = np.asarray([sim.particles[int(j)].z for j in bodyrange])

	cx = sim.particles[1].x
	cy = sim.particles[1].y
	cz = sim.particles[1].z

	x = xp - cx
	y = yp - cy
	z = zp - cz

	lat = np.arctan(z / x)
	lon = np.arctan(y / x)

	x0 = abin * np.cos(lat) * np.cos(lon)
	y0 = bbin * np.cos(lat) * np.sin(lon)
	z0 = cbin * np.sin(lat)

	distancerm = np.sqrt(x ** 2 + y ** 2 + z ** 2)
	surfacerm  = np.sqrt(x0**2 + y0**2 + z0**2)

	inds = np.where(distancerm <= surfacerm)[0]

	Nparts -= len(inds)
	if len(inds) != 0:
		binx = [sim.particles[int(i)].x for i in inds]
		biny = [sim.particles[int(i)].y for i in inds]
		binz = [sim.particles[int(i)].z for i in inds]

		for i in inds[::-1]:
			sim.remove(int(i))
			rminds.append(int(i))
		binland += len(inds)
	# binx = np.where((distancerm <= surfacerm), xp)
	# biny = np.where((distancerm <= surfacerm), yp)
	# binz = np.where((distancerm <= surfacerm), zp)
	#
	# partind = np.where((distancerm <= surfacerm))
	#
	# for p in partind:
	# 	sim.remove(p)
	# 	rminds.append(p)

	Nparts -= len(binx)

	return sim, Nparts, rminds, binland, binx, biny, binz






def rmdatasave(rmx, rmy, rmz, fname):
	data = np.zeros((3, len(rmx)))

	data[0] = rmx
	data[1] = rmy
	data[2] = rmz

	np.savetxt(fname, data, fmt='%.18e', delimiter=' ', newline='\n')





def rmlandedparts(sim, Nplanets, Nparts, atarg, btarg, ctarg, landed, rminds):
	landx = []
	landy = []
	landz = []

	xp = np.asarray([sim.particles[int(j)].x for j in range(1, Nplanets + Nparts)])
	yp = np.asarray([sim.particles[int(j)].y for j in range(1, Nplanets + Nparts)])
	zp = np.asarray([sim.particles[int(j)].z for j in range(1, Nplanets + Nparts)])

	cx = sim.particles[0].x
	cy = sim.particles[0].y
	cz = sim.particles[0].z

	x = xp - cx
	y = yp - cy
	z = zp - cz

	lat = np.arctan(z / x)
	lon = np.arctan(y / x)

	x0 = atarg * np.cos(lat) * np.cos(lon)
	y0 = btarg * np.cos(lat) * np.sin(lon)
	z0 = ctarg * np.sin(lat)

	distancerm = np.sqrt(x ** 2 + y ** 2 + z ** 2)
	surfacerm  = np.sqrt(x0**2 + y0**2 + z0**2)
	inds = np.where(distancerm <= surfacerm)[0]

	Nparts -= len(inds)
	if len(inds) != 0:
		landx = [sim.particles[int(i)].x for i in inds]
		landy = [sim.particles[int(i)].y for i in inds]
		landz = [sim.particles[int(i)].z for i in inds]

		for i in inds[::-1]:
			sim.remove(int(i))
			rminds.append(int(i))
		landed += len(inds)

	return sim, Nparts, rminds, landed, landx, landy, landz, sim.particles[0].x


