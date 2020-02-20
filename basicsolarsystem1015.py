### Basic Solar System Model ###
import rebound
import numpy as np



def addpartcart(sim, pos, vel, Nparts, radmass=np.zeros(3), rho=2e3):
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

		partmass = rho * (4./3.) * np.pi * rad**3

		for p in np.arange(0, Nparts, 1):
			p = int(p)
			sim.add(r=rad[p],
					#m=partmass[p],
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


#	print(x)
	np.savetxt(fname, data, fmt='%.18e', delimiter=' ', newline='\n')



def ejectaconepos(aT, bT, cT, h, lat, lon, beta, tilt, axdir, Nparts, tpos, diskwidth=5):
	'''Setting ejecta cone position'''

	# cone base
	#x0p = aT * np.cos(np.radians(lat)) * np.cos(np.radians(lon))
	#y0p = bT * np.cos(np.radians(lat)) * np.sin(np.radians(lon))
	#z0p = cT * np.sin(np.radians(lat))

	# calculate cone base on surface
	curv = (aT ** 2) / np.sqrt(aT ** 2
							   - (aT ** 2 - cT ** 2) * np.sin(np.radians(lat)) ** 2
							   - (aT ** 2 - bT ** 2) * np.cos(np.radians(lat)) ** 2 * np.sin(np.radians(lon)) ** 2)

	ee = (aT ** 2 - bT ** 2) / aT ** 2
	ex = (aT ** 2 - cT ** 2) / aT ** 2

	x0p = curv * np.cos(np.radians(lat)) * np.cos(np.radians(lon))
	y0p = curv * (1 - ee) * np.cos(np.radians(lat)) * np.sin(np.radians(lon))
	z0p = curv * (1 - ex) * np.sin(np.radians(lat))

	r0 = np.sqrt(x0p**2 + y0p**2 + z0p**2)


	#x0 = r0*np.cos(np.radians(lat)) * np.cos(np.radians(lon))
	#y0 = r0*np.cos(np.radians(lat)) * np.sin(np.radians(lon))
	#z0 = r0*np.sin(np.radians(lat))
	#x0 = r0 * np.sin(np.radians(tilt)) * np.cos(np.radians(axdir))
	#y0 = r0 * np.sin(np.radians(tilt)) * np.sin(np.radians(axdir))
	#z0 = r0 * np.cos(np.radians(tilt))


	# disk radius
	rdisk  = h * np.tan(np.radians((beta) + (diskwidth/2.)))
	rinner = h * np.tan(np.radians((beta) - (diskwidth/2.)))

	r = (rdisk - rinner) * np.random.random_sample(Nparts) + rinner

	thetaprime = 360. * np.random.random_sample(Nparts)

	opang = np.degrees(np.arctan(r / h))

	mag1 = np.sqrt(h**2 + r**2)#np.sqrt(x1**2 + y1**2 + z1new**2)

	phiang  = opang * np.sin(np.radians(thetaprime))
	thetang = opang * np.cos(np.radians(thetaprime))

	#x = mag1 * np.cos(np.radians(lat + phiang)) * np.cos(np.radians(lon + thetang)) + x0
	#y = mag1 * np.cos(np.radians(lat + phiang)) * np.sin(np.radians(lon + thetang)) + y0
	#z = mag1 * np.sin(np.radians(lat + phiang)) + z0

	x = mag1 * np.sin(np.radians(opang)) * np.cos(np.radians(thetaprime))
	y = mag1 * np.sin(np.radians(opang)) * np.sin(np.radians(thetaprime))
	z = mag1 * np.cos(np.radians(opang))

	x1 = x * np.cos(np.radians(lon)) * np.cos(np.radians(90-lat)) \
		 - y * np.sin(np.radians(lon)) \
		 + z * np.cos(np.radians(lon)) * np.sin(np.radians(90-lat))

	y1 = x * np.sin(np.radians(lon)) * np.cos(np.radians(90-lat)) \
		 + y * np.cos(np.radians(lon)) \
		 + z * np.sin(np.radians(lon)) * np.sin(np.radians(90-lat))

	z1 = -x * np.sin(np.radians(90-lat)) \
		 + z * np.cos(np.radians(90-lat))

	x1 += x0p
	y1 += y0p
	z1 += z0p

	pos = (x1 + tpos[0], y1 + tpos[1], z1 + tpos[2])

	s = np.sqrt((x-x0p)**2 + (y-y0p)**2 + (z-z0p)**2)   # distance from particles to cone base

	return pos, x, y, z, opang, thetaprime, r, s




def ejectaconevel(opang, thetaprime, lat, lon, pos, v0, tvel, axtilt, axdir, dt, per=86400., rotation=False):
	'''Setting the ejecta cone velocity'''

	#for i in range(Nparts):
	#vx = (v0 * np.sin(np.radians(lat + phiang)) * np.cos(np.radians(lon + thetang))) #* np.sin(np.radians(opang)) * np.cos(np.radians(thetaprime))
	#vy = (v0 * np.sin(np.radians(lat + phiang)) * np.sin(np.radians(lon + thetang))) #* np.sin(np.radians(opang)) * np.sin(np.radians(thetaprime))
	#vz = (v0 * np.cos(np.radians(lat + phiang))) #* np.cos(np.radians(opang))

	vx = (v0 * np.sin(np.radians(opang)) * np.cos(
		np.radians(thetaprime)))  # * np.sin(np.radians(opang)) * np.cos(np.radians(thetaprime))
	vy = (v0 * np.sin(np.radians(opang)) * np.sin(
		np.radians(thetaprime)))  # * np.sin(np.radians(opang)) * np.sin(np.radians(thetaprime))
	vz = (v0 * np.cos(np.radians(opang)))  # * np.cos(np.radians(opang))


	vx1 = vx * np.cos(np.radians(lon)) * np.cos(np.radians(90-lat)) \
		  - vy * np.sin(np.radians(lon)) \
		  + vz * np.cos(np.radians(lon)) * np.sin(np.radians(90-lat))

	vy1 = vx * np.sin(np.radians(lon)) * np.cos(np.radians(90-lat)) \
		  + vy * np.cos(np.radians(lon)) \
		  + vz * np.sin(np.radians(lon)) * np.sin(np.radians(90-lat))

	vz1 = -vx * np.sin(np.radians(90-lat)) \
		  + vz * np.cos(np.radians(90-lat))



	if rotation == True:
		r = np.sqrt(pos[0] ** 2 + pos[1] ** 2 + pos[2] ** 2)
		omega = (2. * np.pi) / per
		vmag = r * omega

		#vxrot = vmag * np.cos(np.radians(tilt)) * np.cos(np.radians(axdir))
		#vyrot = vmag * np.cos(np.radians(tilt)) * np.sin(np.radians(axdir))
		#vzrot = vmag * np.sin(np.radians(tilt))

		#vx = np.asarray(vx) + vxrot
		#vy = np.asarray(vy) + vyrot
		#vz = np.asarray(vz) + vzrot

		theta = np.radians(axdir)
		phi = np.radians(axtilt)
		alpha = np.radians(omega * dt)

		phix = phi * np.sin(theta)
		phiy = phi * np.cos(theta)

		vx1, vy1, vz1 = rotmatrix3d(vx1, vy1, vz1, alpha, phiy, phix)

		vx1 = np.asarray(vx1)  # + vxrot
		vy1 = np.asarray(vy1)  # + vyrot
		vz1 = np.asarray(vz1)  # + vzrot

	vel = (np.asarray(vx1) + tvel[0], np.asarray(vy1) + tvel[1], np.asarray(vz1) + tvel[2])

	return vel



#
# def ejectaconepos(aT, bT, cT, h, lat, lon, beta, tilt, axdir, Nparts, tpos, diskwidth=10.):
# 	'''Setting ejecta cone position'''
#
# 	# cone base
# 	x0p = aT * np.cos(np.radians(lat)) * np.cos(np.radians(lon))
# 	y0p = bT * np.cos(np.radians(lat)) * np.sin(np.radians(lon))
# 	z0p = cT * np.sin(np.radians(lat))
#
# 	r0 = np.sqrt(x0p**2 + y0p**2 + z0p**2)
#
# 	x0 = r0 * np.sin(np.radians(tilt)) * np.cos(np.radians(axdir))
# 	y0 = r0 * np.sin(np.radians(tilt)) * np.sin(np.radians(axdir))
# 	z0 = r0 * np.cos(np.radians(tilt))
#
#
# 	# disk radius
# 	rdisk  = h * np.tan(np.radians(beta + (diskwidth/2.)))
# 	rinner = h * np.tan(np.radians(beta - (diskwidth/2.)))
#
# 	r = (rdisk - rinner) * np.random.random_sample(Nparts) + rinner
#
# 	thetaprime = 360. * np.random.random_sample(Nparts)
#
# 	opang = np.degrees(np.arctan(r / h))
#
# 	mag1 = h #+ np.sqrt(x0p**2 + y0p**2 + z0p**2) #np.sqrt(h**2 + r**2)#np.sqrt(x1**2 + y1**2 + z1new**2)
#
# 	phiang  = opang * np.sin(np.radians(thetaprime))
# 	thetang = opang * np.cos(np.radians(thetaprime))
#
# 	x = mag1 * np.cos(np.radians(lat + phiang)) * np.cos(np.radians(lon + thetang)) + x0
# 	y = mag1 * np.cos(np.radians(lat + phiang)) * np.sin(np.radians(lon + thetang)) + y0
# 	z = mag1 * np.sin(np.radians(lat + phiang)) + z0
#
#
# 	pos = (x + tpos[0], y + tpos[1], z + tpos[2])
#
# 	s = np.sqrt((x-x0)**2 + (y-y0)**2 + (z-z0)**2)   # distance from particles to cone base
#
# 	return pos, x, y, z, opang, thetaprime, r, s
#
#
#
#
# def ejectaconevel(phiang, thetang, lat, lon, pos, v0, tvel, tilt, axdir, per=86400., rotation=False):
# 	'''Setting the ejecta cone velocity'''
#
# 	vx = (v0 * np.cos(np.radians(lat + phiang)) * np.cos(np.radians(lon + thetang))) #* np.sin(np.radians(opang)) * np.cos(np.radians(thetaprime))
# 	vy = (v0 * np.cos(np.radians(lat + phiang)) * np.sin(np.radians(lon + thetang))) #* np.sin(np.radians(opang)) * np.sin(np.radians(thetaprime))
# 	vz = (v0 * np.sin(np.radians(lat + phiang))) #* np.cos(np.radians(opang))
#
#
# 	if rotation == True:
# 		r = np.sqrt(pos[0] ** 2 + pos[1] ** 2 + pos[2] ** 2)
# 		omega = (2. * np.pi) / per
# 		vmag = r * omega
#
# 		vxrot = vmag * np.sin(np.radians(tilt)) * np.sin(np.radians(axdir))
# 		vyrot = vmag * np.sin(np.radians(tilt)) * np.cos(np.radians(axdir))
# 		vzrot = vmag * np.cos(np.radians(tilt))
#
# 		vx = np.asarray(vx) + vxrot
# 		vy = np.asarray(vy) + vyrot
# 		vz = np.asarray(vz) + vzrot
#
# 	vel = (np.asarray(vx) + tvel[0], np.asarray(vy) + tvel[1], np.asarray(vz) + tvel[2])
#
# 	return vel
#
#
#



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

	lat = np.arctan(z / np.sqrt(x**2 + y**2))
	lon = np.arctan(y / x)

	x0 = abin * np.cos(lat) * np.cos(lon)
	y0 = bbin * np.cos(lat) * np.sin(lon)
	z0 = cbin * np.sin(lat)

	distancerm = np.sqrt(x ** 2 + y ** 2 + z ** 2)
	surfacerm  = np.sqrt(x0**2 + y0**2 + z0**2)

	inds = np.where(distancerm <= surfacerm)[0]

	Nparts -= len(inds)
	if len(inds) != 0:
		inds += 1
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




def rmlandedparts(sim, Nplanets, Nparts, atarg, btarg, ctarg, landed, rminds, axdir=0, axtilt=0, omega=0):
	landx = []
	landy = []
	landz = []

	dt = sim.dt

	bodies = np.linspace(1, Nparts - 2, Nparts-1)

	# particle positions relative to gravitational body
	xp = np.asarray([sim.particles[int(i)].x for i in bodies])
	yp = np.asarray([sim.particles[int(i)].y for i in bodies])
	zp = np.asarray([sim.particles[int(i)].z for i in bodies])

	theta = np.radians(180. - axdir)
	phi = np.radians(axtilt)
	alpha = np.radians(-omega * dt)

	phix = phi * np.cos(theta)
	phiy = phi * np.sin(theta)

	xtilt, ytilt, ztilt = rotmatrix3d(xp, yp, zp, alpha, phiy, phix)

	cx = sim.particles[0].x
	cy = sim.particles[0].y
	cz = sim.particles[0].z

	x = xtilt - cx
	y = ytilt - cy
	z = ztilt - cz

	# calculate point on surface
	lat = np.arcsin(z / np.sqrt(x**2 + y**2 + z**2))
	lon = np.arccos(x / (np.sqrt(x**2 + y**2 + z**2) * np.cos(lat)))

	curv = (atarg**2) / np.sqrt(atarg**2
								- (atarg**2-ctarg**2) * np.sin(lat)**2
								- (atarg**2-btarg**2) * np.cos(lat)**2 * np.sin(lon)**2)

	ee = (atarg**2 - btarg**2) / atarg**2
	ex = (atarg**2 - ctarg**2) / atarg**2

	x0 = curv * np.cos(lat) * np.cos(lon) - cx
	y0 = curv * (1 - ee) * np.cos(lat) * np.sin(lon) - cy
	z0 = curv * (1 - ex) * np.sin(lat) - cz


	# calculate distances
	distancerm = np.sqrt(x ** 2 + y ** 2 + z ** 2)

	A = ((np.cos(lat))**2 * (np.cos(lon))**2) / atarg**2
	B = ((np.cos(lat))**2 * (np.sin(lon)**2)) / btarg**2
	C = ((np.sin(lat))**2) / ctarg**2
	#x0 = atarg * np.cos(lon) * np.cos(lat)
	#y0 = btarg * np.sin(lon) * np.cos(lat)
	#z0 = ctarg * np.sin(lat)

	surfacerm  = np.sqrt((A + B + C) ** -1)
	print(surfacerm)

	#surfacerm = np.sqrt(x0 ** 2 + y0 ** 2 + z0 ** 2)

	inds = np.where(distancerm[:] <= surfacerm[:])[0]

	print('surf dist =', distancerm - surfacerm)

	Nparts -= len(inds)
	if len(inds) != 0:
		inds += Nplanets
		landx = [sim.particles[int(i)].x for i in inds]
		landy = [sim.particles[int(i)].y for i in inds]
		landz = [sim.particles[int(i)].z for i in inds]

		for i in inds[::-1]:
			sim.remove(int(i))
			rminds.append(int(i))
		landed += len(inds)

	return sim, Nparts, rminds, landed, landx, landy, landz, sim.particles[0].x


def rotmatrix3d(x, y, z, alpha, phi, theta):
	'''rotates x, y, z in 3D. make sure angles are in radians.'''
	xp = x * (np.cos(alpha) * np.cos(phi)) \
		 + y * (-np.sin(alpha) * np.cos(theta) + np.cos(alpha) * np.sin(phi) * np.sin(theta)) \
		 + z * (np.sin(alpha) * np.sin(theta) + np.cos(alpha) * np.cos(theta) * np.sin(phi))

	yp = x * (np.sin(alpha) * np.cos(phi)) \
		 + y * (np.cos(alpha) * np.cos(theta) + np.sin(alpha) * np.sin(theta) * np.sin(phi)) \
		 + z * (-np.cos(alpha) * np.sin(theta) + np.sin(alpha) * np.cos(theta) * np.sin(phi))

	zp = x * (-np.sin(phi)) \
		 + y * (np.sin(theta) * np.cos(phi)) \
		 + z * (np.cos(theta) * np.cos(phi))

	return xp, yp, zp

#
# def rmlandedparts(sim, Nplanets, Nparts, atarg, btarg, ctarg, landed, rminds):
# 	landx = []
# 	landy = []
# 	landz = []
#
# 	xp = np.asarray([sim.particles[int(j)].x for j in range(Nplanets, Nplanets + Nparts)])
# 	yp = np.asarray([sim.particles[int(j)].y for j in range(Nplanets, Nplanets + Nparts)])
# 	zp = np.asarray([sim.particles[int(j)].z for j in range(Nplanets, Nplanets + Nparts)])
#
# 	cx = sim.particles[0].x
# 	cy = sim.particles[0].y
# 	cz = sim.particles[0].z
#
# 	x = xp - cx
# 	y = yp - cy
# 	z = zp - cz
#
# 	lat = np.arctan(z / np.sqrt(x**2 + y**2))
# 	lon = np.arctan(y / x)
#
# 	x0 = atarg * np.cos(lat) * np.cos(lon)
# 	y0 = btarg * np.cos(lat) * np.sin(lon)
# 	z0 = ctarg * np.sin(lat)
#
# 	distancerm = np.sqrt(x ** 2 + y ** 2 + z ** 2)
# 	surfacerm  = np.sqrt(x0 ** 2 + y0 ** 2 + z0 ** 2)
# 	inds = np.where(distancerm[:] <= surfacerm[:])[0]
#
# 	Nparts -= len(inds)
# 	if len(inds) != 0:
# 		inds += Nplanets
# 		landx = [sim.particles[int(i)].x for i in inds]
# 		landy = [sim.particles[int(i)].y for i in inds]
# 		landz = [sim.particles[int(i)].z for i in inds]
#
# 		for i in inds[::-1]:
# 			sim.remove(int(i))
# 			rminds.append(int(i))
# 		landed += len(inds)
#
# 	return sim, Nparts, rminds, landed, landx, landy, landz, sim.particles[0].x


def veldist(rtarg, pos, vmin, mtarg, hpart=1., di=1., vi=5., rhoi=2e3, rhot=2e3, ex=2):
	''' Defines velocity distribution based on Housen et al. (1983)'''

	g = (6.67e-11 * mtarg) / (hpart + rtarg) ** 2

	R = 1.16 * (rhoi / rhot) ** (1. / 3.) * di ** .78 * vi ** .44 * g ** -.22

	r    = np.sqrt(pos[0] ** 2 + pos[1] ** 2 + pos[2] ** 2)
	rmin = np.min(r)

	k = (vmin / np.sqrt(g * R)) * (rmin / R) ** ex

	v = k * np.sqrt(g * R) * (r / R) ** -ex
	#v = np.full(int(len(pos[0])), vmin)

	return v