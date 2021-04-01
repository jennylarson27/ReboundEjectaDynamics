### Basic Solar System Model ###
import rebound
import numpy as np
import matplotlib.pyplot as plt
from rebound import hash as h



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
					hash=int(p+1)
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
					hash=int(p+1)
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



def ejectaconepos(mex, aT, bT, cT, h, lat, lon, beta, tilt, axdir, Nparts, tpos, mtarg, rtarg, vi, a, mi, rhoT, rhoI, mu, K1, Ybar=0.):
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

	pos0p = (x0p, y0p, z0p)

	r0 = np.sqrt(x0p**2 + y0p**2 + z0p**2)


	#x0 = r0*np.cos(np.radians(lat)) * np.cos(np.radians(lon))
	#y0 = r0*np.cos(np.radians(lat)) * np.sin(np.radians(lon))
	#z0 = r0*np.sin(np.radians(lat))
	#x0 = r0 * np.sin(np.radians(tilt)) * np.cos(np.radians(axdir))
	#y0 = r0 * np.sin(np.radians(tilt)) * np.sin(np.radians(axdir))
	#z0 = r0 * np.cos(np.radians(tilt))


	# disk radius
	G = 6.67e-11

	g = (G * mtarg) / r0 ** 2

	# transient crater radius

	# Rg = ((3 * K1 * mi) / (np.pi * rhoT)) ** (1. / 3.) \
	# 	 * (((g * a) / vi ** 2) * (rhoT / rhoI) ** (-1. / 3.)
	# 		+ (Ybar / (rhoT * vi ** 2)) ** ((2 + mu) / 2)) ** (-mu / (2 + mu))

	#Rg = K1 * (rhoT / mi ) ** (-1./3.) * (rhoT / rhoI) ** ((2 + mu - 6 * 0.4) / (6 + 3 * mu)) * ((g * a) / vi ** 2) ** (-mu / (2 + mu))
	Rg = ((3. / np.pi) * (mex / rhoT)) ** (1. / 3.)
	print('rg=', Rg)

	#rdisk  = h * np.tan(np.radians(beta))
	#rinner = h * np.tan(np.radians(beta - diskwidth))


	r = (Rg - 1e-0) * np.random.random_sample(Nparts) + 1e-0


	thetaprime = 360. * np.random.random_sample(Nparts)

	radial = np.arctan(r / (r0 + h))
	# make disk shape
	x = r * np.cos(np.radians(thetaprime)) #(r0 + h) * np.cos(np.radians(thetaprime)) * np.cos(radial)
	y = r * np.sin(np.radians(thetaprime)) #r0 + h) * np.sin(np.radians(thetaprime)) * np.cos(radial)
	z = h #(r0 + h) * np.sin(radial)



	#opang = np.degrees(np.arctan(r / h))
	#mag1 = np.sqrt(h**2 + r**2)#np.sqrt(x1**2 + y1**2 + z1new**2)
	#phiang  = opang * np.sin(np.radians(thetaprime))
	#thetang = opang * np.cos(np.radians(thetaprime))
	#x = mag1 * np.cos(np.radians(lat + phiang)) * np.cos(np.radians(lon + thetang)) + x0
	#y = mag1 * np.cos(np.radians(lat + phiang)) * np.sin(np.radians(lon + thetang)) + y0
	#z = mag1 * np.sin(np.radians(lat + phiang)) + z0
	#x = mag1 * np.sin(np.radians(opang)) * np.cos(np.radians(thetaprime))
	#y = mag1 * np.sin(np.radians(opang)) * np.sin(np.radians(thetaprime))
	#z = mag1 * np.cos(np.radians(opang))

	# set velocity position
	vxpos = np.sin(np.radians(beta)) * np.cos(np.radians(thetaprime))
	vypos = np.sin(np.radians(beta)) * np.sin(np.radians(thetaprime))
	vzpos = np.cos(np.radians(beta))


	# positon rotation
	latx = lat * np.sin(np.radians(lon))
	laty = lat * np.cos(np.radians(lon))
	#x1, y1, z1 = rotmatrix3d(x, y, z, np.radians(lon), np.radians(laty), np.radians(latx))
	x1 = x * np.cos(np.radians(lon)) * np.cos(np.radians(90-lat)) \
		 - y * np.sin(np.radians(lon)) \
	 	 + z * np.cos(np.radians(lon)) * np.sin(np.radians(90-lat))

	y1 = x * np.sin(np.radians(lon)) * np.cos(np.radians(90-lat)) \
	 	 + y * np.cos(np.radians(lon)) \
	 	 + z * np.sin(np.radians(lon)) * np.sin(np.radians(90-lat))

	z1 = -x * np.sin(np.radians(90-lat)) \
	 	 + z * np.cos(np.radians(90-lat))



	# velocity rotation
	vx1 = vxpos * np.cos(np.radians(lon)) * np.cos(np.radians(90 - lat)) \
		 - vypos * np.sin(np.radians(lon)) \
		 + vzpos * np.cos(np.radians(lon)) * np.sin(np.radians(90 - lat))

	vy1 = vxpos * np.sin(np.radians(lon)) * np.cos(np.radians(90 - lat)) \
		 + vypos * np.cos(np.radians(lon)) \
		 + vzpos * np.sin(np.radians(lon)) * np.sin(np.radians(90 - lat))

	vz1 = -vxpos * np.sin(np.radians(90 - lat)) \
		 + vzpos * np.cos(np.radians(90 - lat))

	# position
	x1 += x0p
	y1 += y0p
	z1 += z0p

	pos = (x1 + tpos[0], y1 + tpos[1], z1 + tpos[2])

	s = np.sqrt((x-x0p)**2 + (y-y0p)**2 + (z-z0p)**2)   # distance from particles to cone base

	# # velocity
	# vx1 += x0p
	# vy1 += y0p
	# vz1 += z0p

	vpos = (vx1, vy1, vz1)

	return pos, x, y, z, thetaprime, r, s, pos0p, vpos



def ejectaconevel(vpos, pos0p, v0, per=86400., rotation=False):
	'''Sets up ejecta cone velocities 2/24/20'''

	# impact position on surface
	x0p = pos0p[0]
	y0p = pos0p[1]
	z0p = pos0p[2]

	vxhat = vpos[0]
	vyhat = vpos[1]
	vzhat = vpos[2]

	# multiply unit vector by velocity
	vx = v0 * vxhat
	vy = v0 * vyhat
	vz = v0 * vzhat

	vel = (vx, vy, vz)

	return vel



def NOPEejectaconevel(opang, thetaprime, lat, lon, pos, v0, tvel, axtilt, axdir, dt, per=86400., rotation=False):
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



# def rmbinparts(sim, Nplanets, Nparts, abin, bbin, cbin, binland, rminds):
# 	binx = []
# 	biny = []
# 	binz = []
#
# 	bodyrange = np.linspace(0, Nparts+1, Nparts+2)
#
# 	xp = np.asarray([sim.particles[int(j)].x for j in bodyrange])
# 	yp = np.asarray([sim.particles[int(j)].y for j in bodyrange])
# 	zp = np.asarray([sim.particles[int(j)].z for j in bodyrange])
#
# 	pdist = np.sqrt(xp**2 + yp**2 + zp**2)
#
# 	#print('alldistance =', pdist)
#
# 	cx = xp[1]
# 	cy = yp[1]
# 	cz = zp[1]
# 	#print('i hate this =', pdist[1])
#
# 	#print('bindistance =',np.sqrt(cx**2 + cy**2 + cz**2))
#
#
# 	x = xp[2:] - cx
# 	y = yp[2:] - cy
# 	z = zp[2:] - cz
#
# 	lat = np.arctan(z / np.sqrt(x**2 + y**2))
# 	lon = np.arctan(y / x)
#
# 	x0 = abin * np.cos(lat) * np.cos(lon)
# 	y0 = bbin * np.cos(lat) * np.sin(lon)
# 	z0 = cbin * np.sin(lat)
#
# 	distancerm = np.sqrt(x ** 2 + y ** 2 + z ** 2)
# 	surfacerm  = np.sqrt(x0**2 + y0**2 + z0**2)
#
# 	inds = np.where(distancerm <= surfacerm)[0]
#
# 	Nparts -= len(inds)
# 	if len(inds) != 0:
# 		inds += 2
# 		binx = [sim.particles[int(i)].x for i in inds]
# 		biny = [sim.particles[int(i)].y for i in inds]
# 		binz = [sim.particles[int(i)].z for i in inds]
#
# 		for i in inds[::-1]:
# 			sim.remove(int(i))
# 			rminds.append(int(i))
# 		binland += len(inds)
# 	# binx = np.where((distancerm <= surfacerm), xp)
# 	# biny = np.where((distancerm <= surfacerm), yp)
# 	# binz = np.where((distancerm <= surfacerm), zp)
# 	#
# 	# partind = np.where((distancerm <= surfacerm))
# 	#
# 	# for p in partind:
# 	# 	sim.remove(p)
# 	# 	rminds.append(p)
#
# 	Nparts -= len(binx)
#
# 	return sim, Nparts, rminds, binland, binx, biny, binz


def rmbinparts(sim, rbin, Nparts, binland, rminds):
	'''DO NOT USE THIS'''
	binx = []
	biny = []
	binz = []

	dt = sim.dt

	number = 1

	bodies = np.linspace(number, Nparts, Nparts)


	# particle positions relative to gravitational body
	xp = np.asarray([sim.particles[int(i)].x for i in bodies])
	yp = np.asarray([sim.particles[int(i)].y for i in bodies])
	zp = np.asarray([sim.particles[int(i)].z for i in bodies])



	cx = xp[0]   #sim.particles[0].x
	cy = yp[0]   #sim.particles[0].y
	cz = zp[0]   #sim.particles[0].z

	x = xp[1:] - cx
	y = yp[1:] - cy
	z = zp[1:] - cz

	# calculate distances
	distancerm = np.sqrt(x ** 2 + y ** 2 + z ** 2)

	surfacerm  = rbin       #np.sqrt((A + B + C) ** -1)
	print('bin surface=', surfacerm)

	inds = np.where(distancerm - surfacerm <= 0)[0] + number #np.where(distancerm[:] <= surfacerm[:])[0]

	print('binary dist =', distancerm - surfacerm)

	Nparts -= len(inds)
	if len(inds) != 0:
		inds += 1#Nplanets
		binx = [sim.particles[int(i)].x for i in inds]
		biny = [sim.particles[int(i)].y for i in inds]
		binz = [sim.particles[int(i)].z for i in inds]

		for i in inds[::-1]:
			sim.remove(int(i))
			rminds.append(int(i))
		binland += len(inds)

	return sim, Nparts, rminds, binland, binx, biny, binz



def rmdatasave(rmx, rmy, rmz, fname):
	'''USE THIS'''
	data = np.zeros((3, len(rmx)))

	data[0] = rmx
	data[1] = rmy
	data[2] = rmz

	np.savetxt(fname, data, fmt='%.18e', delimiter=' ', newline='\n')



def rmparticles(sim, Nparts, landed, inttime, aT, aB, condit, rminds):
	'''remove based on names 10/21/2020'''

	if aB == 0:
		binary = False
		indshift = 1
	else:
		binary = True
		indshift = 2

	nhash = [p.hash for p in sim.particles]
	print(nhash)



	cx = sim.particles[h('target')].x
	cy = sim.particles[h('target')].y
	cz = sim.particles[h('target')].z


	landx, landy, landz = ([], [], [])




	print('len', len(rminds))

	#if len(remaining) == 0:
	leftparts = np.linspace(1, Nparts, Nparts)
		#leftparts = [1,3,5]
	if len(rminds) != 0:
		removeinds = np.asarray(rminds) - 1
		leftparts = np.delete(leftparts, removeinds)#remaining
		#print('left=', leftparts)


	print(leftparts)



	# for n in reversed(leftparts):
	#
	# 	n = int(n)
		#print(n)

		#nhash = sim.particles[n].hash
		#print(nhash)


	x = [sim.particles[h(n)].x for n in leftparts]
	y = [sim.particles[h(n)].y for n in leftparts]
	z = [sim.particles[h(n)].z for n in leftparts]

	xT = x - cx
	yT = y - cy
	zT = z - cz

	dist = np.sqrt(xT ** 2 + yT ** 2 + zT ** 2)

	surf = aT

	inds = np.where(dist <= surf)[0]

	for i in reversed(inds):
		if i > 1:
			landx.append(sim.particles[i].x)
			landy.append(sim.particles[i].y)
			landz.append(sim.particles[i].z)

			landed += 1
			Nparts -= 1

			sim.remove(i)

	if binary == True:
		bx = sim.particles[h('bin')].x
		by = sim.particles[h('bin')].y
		bz = sim.particles[h('bin')].z

		blandx, blandy, blandz = ([], [], [])

		xB = x - bx
		yB = y - by
		zB = z - bz

		bdist = np.sqrt(xB ** 2 + yB ** 2 + zB ** 2)

		bsurf = aB
		indsb = np.where(bdist <= bsurf)[0]

		for i in reversed(indsb):
			if i > 1:
				blandx.append(sim.particles[i].x)
				blandy.append(sim.particles[i].y)
				blandz.append(sim.particles[i].z)

				landed +=1
				Nparts -=1

				sim.remove(i)




	if dist <= surf:
		print('rm=', n)

		landx.append(x)
		landy.append(y)
		landz.append(z)

		sim.remove(hash=n)

		landed += 1
		Nparts -= 1
		rminds.append(n)




	if binary == True:
		xB = x - bx
		yB = y - by
		zB = z - bz

		bdist = np.sqrt(xB ** 2 + yB ** 2 + zB ** 2)

		bsurf = aB

		if bdist <= bsurf:
			print('rmb=', n)
			blandx.append(x)
			blandy.append(y)
			blandz.append(z)

			sim.remove(hash=n)
			landed += 1
			Nparts -= 1
			rminds.append(n)




	rmdatasave(landx, landy, landz, 'rmlandpartdata' + str(inttime) + condit + '.txt')

	if binary == True:
		rmdatasave(blandx, blandy, blandz, 'rmblandpartdata' + str(inttime) + condit + '.txt')

	print('landed=', landed)
	#print('rminds',rminds)
	#print(sim.particles[h('bin')])

	return sim, landed, Nparts, rminds





def rmland(sim, Nparts, atarg, btarg, ctarg, landed, inttime, condit, axdir, axtilt, per, timestep, abin=0, rotation=False):
	'''USE THIS'''
	if abin > 0:
		binary = True
		number = 2
	else:
		binary = False
		number = 1

	bodies = np.linspace(0, sim.N-1, sim.N)

	xp = np.asarray([sim.particles[int(i)].x for i in bodies])
	yp = np.asarray([sim.particles[int(i)].y for i in bodies])
	zp = np.asarray([sim.particles[int(i)].z for i in bodies])


	# rotate particles back to where they should be above the surface
	if rotation == True:
		omega = 360. / per
		dt = sim.dt * timestep
		theta = np.radians(180. - axdir)
		phi = np.radians(axtilt)
		alpha = np.radians(-omega * dt)

		phix = phi * np.cos(theta)
		phiy = phi * np.sin(theta)

		xtilt, ytilt, ztilt = rotmatrix3d(xp, yp, zp, alpha, phiy, phix)

		# distance from center of body
		cx = xtilt - sim.particles[0].x
		cy = ytilt - sim.particles[0].y
		cz = ztilt - sim.particles[0].z

	elif rotation == False:
		cx = xp - xp[0] #sim.particles[0].x
		cy = yp - yp[0] #sim.particles[0].y
		cz = zp - zp[0] #sim.particles[0].z

	# calculate point on surface
	lat = np.arcsin(cz / np.sqrt(cx ** 2 + cy ** 2 + cz ** 2))
	lon = np.arccos(cx / (np.sqrt(cx ** 2 + cy ** 2 + cz ** 2) * np.cos(lat)))

	A = ((np.cos(lat)) ** 2 * (np.cos(lon)) ** 2) / atarg ** 2
	B = ((np.cos(lat)) ** 2 * (np.sin(lon) ** 2)) / btarg ** 2
	C = ((np.sin(lat)) ** 2) / ctarg ** 2

	surf = np.sqrt((A + B + C) ** -1)
	dist = np.sqrt(cx ** 2 + cy ** 2 + cz ** 2)
	#surf = atarg



	print('zero', np.sqrt(sim.particles[0].x**2 + sim.particles[0].y**2 + sim.particles[0].z**2))
	print('dist1', np.sqrt(xp ** 2 + yp ** 2 + zp ** 2))

	# if binary == False:
	# 	inds = np.where((dist - surf) <= 0)[0]
	# 	inds += number
	#
	# elif binary == True:
		# bsurf = abin
		# bx = xp - xp[1] #sim.particles[1].x
		# by = yp - yp[1] #sim.particles[1].y
		# bz = zp - zp[1] #sim.particles[1].z
		# bdist = np.sqrt(bx ** 2 + by ** 2 + bz ** 2)
		# print('binary=', bdist)
		# #np.concatenate((np.where(dist <= surf)[0], np.where(bdist <= bsurf)[0]))
		# indsb = np.where(bdist <= bsurf)[0]
		#inds += number - 1
		#indsb += number - 1

	landx = []
	landy = []
	landz = []
	landbx, landby, landbz = ([], [], [])

	inds = np.where(dist <= surf)[0]

	#print('shape', inds)
	if inds.shape != 0:
		for i in reversed(inds):
			if i > 1:
				print('ind',i)
				landx.append(sim.particles[int(i)].x)
				landy.append(sim.particles[int(i)].y)
				landz.append(sim.particles[int(i)].z)

				sim.remove(int(i))
				landed += 1
				Nparts -= 1
		print('land =', landed)

	rmdatasave(landx, landy, landz, 'rmlandpartdata' + str(inttime) + condit + '.txt')


	if binary == True:
		#print indsb

		bodies = np.linspace(0, sim.N - 1, sim.N)

		xp = np.asarray([sim.particles[int(i)].x for i in bodies])
		yp = np.asarray([sim.particles[int(i)].y for i in bodies])
		zp = np.asarray([sim.particles[int(i)].z for i in bodies])

		bsurf = abin
		bx = xp - xp[1]  # sim.particles[1].x
		by = yp - yp[1]  # sim.particles[1].y
		bz = zp - zp[1]  # sim.particles[1].z
		bdist = np.sqrt(bx ** 2 + by ** 2 + bz ** 2)
		print('binary=', bdist)
		# np.concatenate((np.where(dist <= surf)[0], np.where(bdist <= bsurf)[0]))
		indsb = np.where(bdist <= bsurf)[0]

		if indsb.shape != 0:
			for ib in reversed(indsb):
				if ib > 1:
					landbx.append(sim.particles[int(ib)].x)
					landby.append(sim.particles[int(ib)].y)
					landbz.append(sim.particles[int(ib)].z)

					sim.remove(int(ib))

					landed += 1
					Nparts -= 1

		rmdatasave(landbx, landby, landbz, 'rmblandpartdata' + str(inttime) + condit + '.txt')

	return sim, landed, Nparts


def rmlandedparts(sim, Nparts, atarg, btarg, ctarg, landed, rminds, binary=False, rbin=0, axdir=0, axtilt=0, omega=0):
	'''DO NOT USE'''
	landx = []
	landy = []
	landz = []

	dt = sim.dt

	if binary == True:
		number = 2
		bodies = np.linspace(number, Nparts+1, Nparts+1)
	else:
		number = 1
		bodies = np.linspace(number, Nparts, Nparts)


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

	cx = 0#sim.particles[0].x
	cy = 0#sim.particles[0].y
	cz = 0#sim.particles[0].z

	x = xtilt - cx
	y = ytilt - cy
	z = ztilt - cz

	# calculate point on surface
	lat = np.arcsin(z / np.sqrt(x**2 + y**2 + z**2))
	lon = np.arccos(x / (np.sqrt(x**2 + y**2 + z**2) * np.cos(lat)))

	#curv = (atarg**2) / np.sqrt(atarg**2
	#							- (atarg**2-ctarg**2) * np.sin(lat)**2
	#							- (atarg**2-btarg**2) * np.cos(lat)**2 * np.sin(lon)**2)

	#ee = (atarg**2 - btarg**2) / atarg**2
	#ex = (atarg**2 - ctarg**2) / atarg**2

	#x0 = curv * np.cos(lat) * np.cos(lon) - cx
	#y0 = curv * (1 - ee) * np.cos(lat) * np.sin(lon) - cy
	#z0 = curv * (1 - ex) * np.sin(lat) - cz


	# calculate distances
	distancerm = np.sqrt(x ** 2 + y ** 2 + z ** 2)

	A = ((np.cos(lat))**2 * (np.cos(lon))**2) / atarg**2
	B = ((np.cos(lat))**2 * (np.sin(lon)**2)) / btarg**2
	C = ((np.sin(lat))**2) / ctarg**2
	#x0 = atarg * np.cos(lon) * np.cos(lat)
	#y0 = btarg * np.sin(lon) * np.cos(lat)
	#z0 = ctarg * np.sin(lat)

	surfacerm  = np.sqrt((A + B + C) ** -1)
	print('surface=', surfacerm)

	#surfacerm = np.sqrt(x0 ** 2 + y0 ** 2 + z0 ** 2)

	# set up binary case
	if binary == True:
		bx = sim.particles[1].x
		by = sim.particles[1].y
		bz = sim.particles[1].z
		print('bindist=', np.sqrt(bx**2 + by**2 + bz**2))

		binx = xp - bx
		biny = yp - by
		binz = zp - bz

		bdist = np.sqrt(binx**2 + biny**2 + binz**2)
		#print('bdist =', bdist)

	inds = []

#	for i in range(len(distancerm)):
#		if distancerm[i] - surfacerm[i] <= 0:
#			inds.append(i+number)
#		if binary == True:
#			if bdist[i] - rbin <= 0:
#				inds.append(i)#+number)

		#inds = np.where(distancerm - surfacerm <= 0 or bdist - rbin <= 0)[0] + number

	inds = np.where(distancerm - surfacerm <= 0)[0] #+ number #np.where(distancerm[:] <= surfacerm[:])[0]

	print('surf dist =', distancerm - surfacerm)
	print('indslen =', len(inds))

	Nparts -= len(inds)
	if len(inds) != 0:
		landx = [xp[i] for i in inds] #[sim.particles[int(i)].x for i in inds]
		landy = [yp[i] for i in inds] #[sim.particles[int(i)].y for i in inds]
		landz = [zp[i] for i in inds] #[sim.particles[int(i)].z for i in inds]

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


def veldist(sim, r, mtarg, rtarg, mex, mu, rhoT, Cvps=0., Ctg=0.8, Ybar=0.): #vi, a, mi, rhoT, rhoI, mu, K1, Cvps=0., Ctg=0.8, Ybar=0.):
	''' Defines velocity distribution based on Richardson (2007)'''

	G = sim.G

	g = (G * mtarg) / rtarg ** 2

	# transient crater radius

	# Rg = ((3. / np.pi) * (K1 * (mi / rhoT)) * (((g * a) / vi ** 2) * (rhoT / rhoI) ** (-1./3.)
	# 										   + (Ybar / (rhoT * vi ** 2)) ** ((2 + mu) / 2)) ** (-3 * mu / (2 + mu))) ** (1. / 3.)

#	Rg = K1 * (rhoT / mi) ** (-1. / 3.) * (rhoT / rhoI) ** ((2 + mu - 6 * 0.4) / (6 + 3 * mu)) * (
#																								 (g * a) / vi ** 2) ** (
#																								 -mu / (2 + mu))

	Rg = ((3. / np.pi) * (mex / rhoT)) ** (1. / 3.)

	#print(Rg)

	# velocity distribution (y-axis)

	Cvpg = (np.sqrt(2) / Ctg) * (mu / (mu + 1))

	ve = Cvpg * np.sqrt(g * Rg) * (r / Rg) ** (-1 / mu)


	veff = np.sqrt(ve ** 2 - Cvpg ** 2 * (g * r) - Cvps ** 2 * (Ybar / rhoT))




	# g = (6.67e-11 * mtarg) / (hpart + rtarg) ** 2   # constant
	#
	# R = 1.16 * (rhoi / rhot) ** (1. / 3.) * di ** .78 * vi ** .44 * g ** -.22   # constant
	#
	# #r    = np.sqrt(pos[0] ** 2 + pos[1] ** 2 + pos[2] ** 2)    # dependent on particle position
	# rmin = np.min(r)
	# rmax = np.max(r)
	# print(len(r))
	#
	# vmin = .01 * v0
	# vmax = 2 * v0
	#
	# k = (vmin / np.sqrt(g * R)) * (rmax / R) ** -ex   # constant
	#
	#
	# ex = np.log(vmax / vmin) / np.log((rmin * rmax) / R**2)
	# print(ex)
	#
	# v = (vmin) * ((rmin * r) / R**2) ** ex  #k * np.sqrt(g * R) * (r / R) ** -ex     #k * np.sqrt(g * R) *
	# #v = np.full(int(len(pos[0])), vmin)
	#
	# v = (v / np.max(v)) * vmax*1.3

	# plot things

	# fig, ax = plt.subplots()
	# fig.subplots_adjust(hspace=.4, bottom=.25, left=.25)
	#
	# ax.scatter(r, veff / np.max(veff))  # / np.max(veff)))
	# #ax.plot(r,np.zeros(len(r)), color='black', marker='_')
	# #ax.set_xscale('log')
	# ax.set_yscale('log')
	# ax.set_ylabel(r'$v_{eff}$ [m/sec]')
	# ax.set_xlabel(r'$R_{g}$ [m]')
	# plt.savefig('veltest.png')
	# plt.close()


	data = np.zeros((2, len(r)))

	data[0, :] = r
	data[1, :] = veff
	np.savetxt('veldistinit.txt', data, fmt='%.18e', delimiter=' ', newline='\n')

	return veff