### Package of additional effects such as size distribution, radiation pressure, binaries, etc ###

import rebound
import numpy as np
import math as ma
import basicsolarsystem1015 as bss
import csv



### Ellipsoidal Gravitational Acceleration ###

def spheregrav(sim, surfpos, Mtarg, Nparts, binary=False):
	''' ***Might not need this afterall*
	Calculates gravity of spherical body

	Currently in progress
	Oct 4, 2019 - testing
	'''

	G = sim.G

	inds = np.linspace(0, Nparts, Nparts)

	if binary == True:
		gravbod = 1
		bodies	= np.delete(inds, gravbod, 0)
	else:
		gravbod = 0
		bodies	= np.delete(inds, gravbod, 0)


	c = sim.particles[gravbod]

	### THIS IS WHERE WE PUT IN THE ROTATION CALCULATION ###
	# particle positions relative to gravitational body
	xsurf = surfpos[0] #- c.x #np.asarray([sim.particles[int(i)].x for i in bodies]) - c.x
	ysurf = surfpos[1] #- c.y #np.asarray([sim.particles[int(i)].y for i in bodies]) - c.y
	zsurf = surfpos[2] #- c.z #np.asarray([sim.particles[int(i)].z for i in bodies]) - c.z

	xp = xsurf[2:]
	yp = ysurf[2:]
	zp = zsurf[2:]

	#updated 4/17/20

	r = np.sqrt(xp**2 + yp**2 + zp**2)
	#print('surfpos=', r)

	amag = -(G * Mtarg) / r**2

	#print('amag =', amag[1])

	ax = amag * (xp / r)
	ay = amag * (yp / r)
	az = amag * (zp / r)

	return ax, ay, az


def ellipgrav(sim, surfpos, Mtarg, Nparts, a1, b1, c1, dt, omega, axdir, axtilt, binary=False):
	''' Calculates gravity of ellipsoidal body

	Sept 17, 2019 - in progress
	Oct 4, 2019 - testing
	Oct 9, 2019 - still writing and testing
	'''

	G = sim.G

	inds = np.linspace(0, Nparts-1, Nparts)

	if binary == True:
		gravbod = 1
		bodies = np.delete(inds, gravbod, 0)
	else:
		gravbod = 0
		bodies = np.delete(inds, gravbod, 0)

	c = sim.particles[gravbod]

	# particle positions relative to gravitational body
	xp = np.asarray([sim.particles[int(i)].x for i in bodies]) - c.x
	yp = np.asarray([sim.particles[int(i)].y for i in bodies]) - c.y
	zp = np.asarray([sim.particles[int(i)].z for i in bodies]) - c.z

	theta = np.radians(180. - axdir)
	phi   = np.radians(180. - axtilt)
	#alpha = np.radians(180. - omega * dt)	# degrees rotated per timestep

	phix = phi * np.sin(theta)
	phiy = phi * np.cos(theta)

	### THIS IS WHERE WE WANT TO PUT IN THE ROTATION FUNCTION THING FOR POSITIONS ###
	# These are the particle positions relative to the surface
	xsurf, ysurf, zsurf = surfpos #bss.rotmatrix3d(xp, yp, zp, alpha, phiy, phix)
	xtilt = xsurf[2:]
	ytilt = ysurf[2:]
	ztilt = zsurf[2:]


	r = np.sqrt(xtilt ** 2 + ytilt ** 2 + ztilt ** 2)

	# unit vector that points to the center of the body from the particle position
	xhat = xtilt / r
	yhat = ytilt / r
	zhat = ztilt / r


	# Rotate hat vector back to global frame
	#theta = np.radians(axdir)
	#phi   = np.radians(axtilt)
	#alpha = np.radians(omega * dt)

	#phix = phi * np.sin(theta)
	#phiy = phi * np.cos(theta)

	#xhatp, yhatp, zhatp = bss.rotmatrix3d(xhat, yhat, zhat, alpha, phiy, phix)


	# calculate gravitational acceleration from shell potential
	C1 = (a1 ** 2) / (b1 ** 2)
	C2 = (a1 ** 2) / (c1 ** 2)

	apos = np.sqrt(xtilt ** 2 + C1 * ytilt ** 2 + C2 * ztilt ** 2)

	h = np.sqrt(a1**2 - b1**2)
	k = np.sqrt(a1**2 - c1**2)

	aellip = -(G * Mtarg) / np.sqrt(np.abs((apos**2 - h**2) * (apos**2 - k**2)))
	#print('denom', apos**2 - h**2, apos**2 - k**2)


	axellip = aellip * xhat
	ayellip = aellip * yhat
	azellip = aellip * zhat
	#print('axellip =', axellip[1])
	#print('ayellip =', ayellip[1])
	#print('azellip =', azellip[1])

	#print('aellip =', -np.sqrt(axellip[1]**2 + ayellip[1]**2 + azellip[1]**2))

	return axellip, ayellip, azellip



def rmaddgrav(a1, a2):
	''' Subtracts additional gravitational acceleration due to default or crater/shape model
	Sept 17, 2019 - in progress
	Oct 4, 2019 - testing '''

	ax1 = a1[0]
	ay1 = a1[1]
	az1 = a1[2]

	ax2 = a2[0]
	ay2 = a2[1]
	az2 = a2[2]

	# remove additional gravity
	axnet = ax1 - ax2
	aynet = ay1 - ay2
	aznet = az1 - az2

	return axnet, aynet, aznet



def addnetgrav(sim, Mtarg, a1, b1, c1, Nplanets, Nparts, dt, omega, timestep, axdir, axtilt, binary=False):
	''' Add net gravitational acceleration to simulation
	Sept 17, 2019 - in progress
	Oct 4, 2019 - testing '''

	per = omega * 360.

	# Put in calculation relative to surface here
	xsurf, ysurf, zsurf = rotpos(sim, per, Nplanets, Nparts, axtilt, axdir, timestep)

	asphere = ellipgrav(sim, (xsurf, ysurf, zsurf), Mtarg, Nparts, a1, a1, a1, dt, omega, axdir, axtilt, binary)#spheregrav(sim, (xsurf, ysurf, zsurf), Mtarg, Nparts, binary)

	aellip	= ellipgrav(sim, (xsurf, ysurf, zsurf), Mtarg, Nparts, a1, b1, c1, dt, omega, axdir, axtilt, binary)

	axnet, aynet, aznet = rmaddgrav(aellip, asphere)#aellip[0], aellip[1], aellip[2]  #


	inds = np.linspace(0, Nparts-1, Nparts)

	if binary == True:
		gravbod = 1
		bodies = np.delete(inds, gravbod, 0)
	else:
		gravbod = 0
		bodies = np.delete(inds, gravbod, 0)

	c = sim.particles[gravbod]

	xp = np.asarray([sim.particles[int(i)].x for i in bodies]) - c.x
	yp = np.asarray([sim.particles[int(i)].y for i in bodies]) - c.y
	zp = np.asarray([sim.particles[int(i)].z for i in bodies]) - c.z

	# Centripetal force due to rotation
	#dist = np.sqrt(xp**2 + yp**2)
	#arot = omega ** 2 * dist
	#rotxhat = xp / dist
	#rotyhat = yp / dist
### Don't forget to tilt the xhat and yhat
	#axrot = arot * rotxhat
	#ayrot = arot * rotyhat

	correction = 1
	# add net acceleration to given accelerations
	for i in bodies:
		p = sim.particles[int(i)]

		if binary == True:
			if int(i) == 0:
				correction = 0
			else:
				correction = 1
		#print(correction)
		p.ax += (axnet[int(i)-correction] )#+ axrot[int(i)-correction])
		p.ay += (aynet[int(i)-correction] )#+ ayrot[int(i)-correction])
		p.az += aznet[int(i)-correction]



	#ax = np.asarray([sim.particles[int(i)].ax for i in bodies])
	#ay = np.asarray([sim.particles[int(i)].ay for i in bodies])
	#az = np.asarray([sim.particles[int(i)].az for i in bodies])


	return sim

#########################################################################









def rotvel(sim, per, lat, pos, axtilt, axdir,):# M, a1, b1, c1, binary=False):
	''' This function works! Updated 4/18/2020. '''
	#global sim

	dt = sim.dt

	omega = (2 * np.pi) / per

	rot = omega * dt   # degrees rotated in one time step at the equator

	#bodyrange = np.arange(Nparts+3)
	#print bodyrange
	#exclude = [1]
	#bodyrange = np.delete(bodies, exclude)

	x = pos[0] # np.asarray([sim.particles[int(j)].x for j in bodyrange])
	y = pos[1] # np.asarray([sim.particles[int(j)].y for j in bodyrange])
	z = pos[2] # np.asarray([sim.particles[int(j)].z for j in bodyrange])

	r = np.sqrt(x ** 2 + y ** 2 + z ** 2)

	# Particles at the surface are rotating at omega as well.
	# linear velocity:
	vlin = omega * r * np.cos(np.radians(lat))

	#sans tilt:
	lon = np.arccos(x / r)
	vrot = (-vlin * np.cos(np.radians(lon)), vlin * np.sin(np.radians(lon)), 0.)

	if axtilt == 0.:
		vxrot = vrot[0]
		vyrot = vrot[1]
		vzrot = vrot[2]

	else:
		# body is tilted some amount about both the x and y axes
		theta = axtilt * np.sin(np.radians(axdir))
		phi   = axtilt * np.cos(np.radians(axdir))

		# the initial rotation velocity is tilted because of the tilted rotation axis.
		vxrot = vrot[0] * np.cos(np.radians(phi)) \
				+ vrot[1] * np.sin(np.radians(phi)) * np.sin(np.radians(theta)) \
				+ vrot[2] * np.sin(np.radians(phi)) * np.cos(np.radians(theta))

		vyrot = vrot[1] * np.cos(np.radians(theta)) \
				- vrot[2] * np.sin(np.radians(theta))

		vzrot = -vrot[0] * np.sin(np.radians(phi)) \
				+ vrot[1] * np.cos(np.radians(phi)) * np.sin(np.radians(theta)) \
				+ vrot[2] * np.cos(np.radians(phi)) * np.cos(np.radians(theta))

	#print(vxrot, vyrot, vzrot)
	#print('vlin=',vlin)

	return vxrot, vyrot, vzrot	    # return the tilted velocities. these are in the global frame



def rotpos(sim, per, Nplanets, Nparts, axtilt, axdir, timestep):
	'''defines the fake location where the particles are in relation to the surface based on rotation'''
	#global sim

	dt = sim.dt

	omega = 360. / per

	rot = omega * dt   # degrees rotated in one time step at the equator

	bodyrange = np.arange(0, Nplanets + Nparts)
	#exclude = [1]
	#bodyrange = np.delete(bodies, exclude)

	x = np.asarray([sim.particles[int(j)].x for j in bodyrange])
	y = np.asarray([sim.particles[int(j)].y for j in bodyrange])
	z = np.asarray([sim.particles[int(j)].z for j in bodyrange])

	# body is tilted some amount about both the x and y axes
	theta = axtilt * np.sin(np.radians(axdir))
	phi   = axtilt * np.cos(np.radians(axdir))

	# positions shifted to vertical axis
	xvp = x * np.cos(np.radians(-phi)) \
		  + y * np.sin(np.radians(-phi)) * np.sin(np.radians(-theta)) \
		  + z * np.sin(np.radians(-phi)) * np.cos(np.radians(-theta))

	yvp = y * np.cos(np.radians(-theta)) \
		  - z * np.sin(np.radians(-theta))

	zvp = -x * np.sin(np.radians(-phi)) \
		  + y * np.cos(np.radians(-phi)) * np.sin(np.radians(-theta)) \
		  + z * np.cos(np.radians(-phi)) * np.cos(np.radians(-theta))

	# how far the body has rotated total
	rottot = -rot * timestep

	# locations relative to the surface
	xp = xvp * np.cos(np.radians(rottot)) \
		 - yvp * np.sin(np.radians(rottot))

	yp = xvp * np.sin(np.radians(rottot)) \
		 + yvp * np.cos(np.radians(rottot))

	zp = zvp

	return xp, yp, zp



### SHAPE MODEL ###
def shapemass(vert, facet, M, layers=5):
	x = []
	y = []
	z = []
	v1 = []
	v2 = []
	v3 = []

	# read in data files
	with open(vert) as vertices:
		vertice = csv.reader(vertices, delimiter='\t')
		rownum = 0
		for row in vertice:
			#print(float(row[2]))
			#print(row)
			x.append(float(row[0]))
			y.append(float(row[1]))
			z.append(float(row[2]))
			rownum += 1
			#print(rownum)

	with open(facet) as face:
		faces = csv.reader(face, delimiter='\t')

		for row in faces:
			v1.append(float(row[0]))
			v2.append(float(row[1]))
			v3.append(float(row[2]))

	x = np.asarray(x)
	y = np.asarray(y)
	z = np.asarray(z)
	v1 = np.asarray(v1) - 1
	v2 = np.asarray(v2) - 1
	v3 = np.asarray(v3) - 1

	N = len(v1)  # number of faces
	facerange = np.linspace(0, N, N)
	Atri = np.empty(shape=(layers, N))
	GiFx = np.empty(shape=(layers, N))
	GiFy = np.empty(shape=(layers, N))
	GiFz = np.empty(shape=(layers, N))

	v1x = np.empty(shape=(layers, N))
	v1y = np.empty(shape=(layers, N))
	v1z = np.empty(shape=(layers, N))

	v2x = np.empty(shape=(layers, N))
	v2y = np.empty(shape=(layers, N))
	v2z = np.empty(shape=(layers, N))

	v3x = np.empty(shape=(layers, N))
	v3y = np.empty(shape=(layers, N))
	v3z = np.empty(shape=(layers, N))


	for layer in np.linspace(0, layers-1, layers):
		for i in v1:
			# Define coordinates of the vertex of each face
			v1x[int(layer), :] = ((layer + 1) * x[int(i)]) / layers
			v1y[int(layer), :] = ((layer + 1) * y[int(i)]) / layers
			v1z[int(layer), int(i)] = ((layer + 1) * z[int(i)]) / layers

		for j in v2:
			v2x[int(layer), :] = ((layer + 1) * x[int(j)]) / layers
			v2y[int(layer), :] = ((layer + 1) * y[int(j)]) / layers
			v2z[int(layer), :] = ((layer + 1) * z[int(j)]) / layers

		for k in v3:
			v3x[int(layer), :] = ((layer + 1) * x[int(k)]) / layers
			v3y[int(layer), :] = ((layer + 1) * y[int(k)]) / layers
			v3z[int(layer), :] = ((layer + 1) * z[int(k)]) / layers


	# Calculate area of the faces
	ax = v2x - v1x
	ay = v2y - v1y
	az = v2z - v1z

	bx = v3x - v1x
	by = v3y - v1y
	bz = v3z - v1z

	cx = v3x - v2x
	cy = v3y - v2y
	cz = v3z - v2z

	axax = ax * ax	#np.multiply(ax, ax)
	bxbx = bx * bx	#np.multiply(bx, bx)
	cxcx = cx * cx	#np.multiply(cx, cx)
	ayay = ay * ay	#np.multiply(ay, ay)
	byby = by * by	#np.multiply(by, by)
	cycy = cy * cy	#np.multiply(cy, cy)
	azaz = az * az	#np.multiply(az, az)
	bzbz = bz * bz	#np.multiply(bz, bz)
	czcz = cz * cz	#np.multiply(cz, cz)

	# Surface area; Heron Method
	a = np.sqrt(axax + ayay + azaz)	 # triangle edge value a
	b = np.sqrt(bxbx + byby + bzbz)	 # triangle edge value b
	c = np.sqrt(cxcx + cycy + czcz)	 # triangle edge value c

	p = (a + b + c) / 2.  # semi-perimeter
	pa = p - a
	pb = p - b
	pc = p - c

	Atri = np.sqrt(p * pa * pb * pc)  # area of each triangle
	Alayer = np.sum(Atri, axis=1)
	Atot = np.sum(Atri)  # total surface area of body


	# Barycenter of faces
	XgF = (v1x + v2x + v3x) / 3.  # median x
	YgF = (v1y + v2y + v3y) / 3.  # median y
	ZgF = (v1z + v2z + v3z) / 3.  # median z



	Xg = np.empty(shape=(layers, N))
	Yg = np.empty(shape=(layers, N))
	Zg = np.empty(shape=(layers, N))

	# Barycenter of Sub-Tetrahedron (prisms 1 through layers-1)
	for layer in np.linspace(0, layers-2, layers-1):
		Xg[int(layer)] = (XgF[int(layer)] + XgF[int(layer+1)]) / 2.
		Yg[int(layer)] = (YgF[int(layer)] + YgF[int(layer+1)]) / 2.
		Zg[int(layer)] = (ZgF[int(layer)] + ZgF[int(layer+1)]) / 2.

	# Barycenter of inner sub-tetrahedron
	Xg[-1] = (v1x[-1] + v2x[-1] + v3x[-1]) / 4.
	Yg[-1] = (v1y[-1] + v2y[-1] + v3y[-1]) / 4.
	Zg[-1] = (v1z[-1] + v2z[-1] + v3z[-1]) / 4.


	# Tetrahedron volumes
	H = np.sqrt(XgF ** 2 + YgF ** 2 + ZgF ** 2)   # height of tetrahedron
	Volall = np.empty(shape=(layers, N))
	for l in np.linspace(0, layers-1, layers):
		Volall[int(l)] = (Alayer[int(l)] * H[int(l)]) / 3


	Vfin = np.empty(shape=(layers, N))
	Vfin[-1] = Volall[-1]
	for l in np.arange(layers-2, -1, -1):
		Vfin[l] = Volall[l] - np.sum(Volall[-1:l:-1])

	Vlayer = np.sum(Vfin, axis=1)
	Vtotal = np.sum(Vfin)




	# # Mass of each sub-face
	# MiA = np.empty(shape=(layers, N))
	# for l in np.linspace(0, layers-1, layers):
	#	MiA[int(l)] = (M / Atot) * Atri[int(l)]
	#
	# MiAprova = np.sum(MiA, axis=1)   # Mass of each layer
	# MprovaAF = np.sum(MiAprova)	   # Total mass overall


	# Mass related to each volume
	Mi = (M / Vtotal) * Vfin

	return Mi, Xg, Yg, Zg






	# Do not use the below functions.

	#Ignore binary stuff for now. ce dan deal with that later when we know what we are doing.
	# if binary == True:
	#	bodies = np.arange(0, Nplanets + Nparts)
	#	exclude = [1]
	#	bodyrange = np.delete(bodies, exclude)
	#
	#	x = np.asarray([sim.particles[int(j)].x for j in bodyrange])
	#	y = np.asarray([sim.particles[int(j)].y for j in bodyrange])
	#	z = np.asarray([sim.particles[int(j)].z for j in bodyrange])
	#
	#	cx = sim.particles[1].x
	#	cy = sim.particles[1].y
	#	cz = sim.particles[1].z


	#else:

	# Get particle locations
	# bodyrange = np.arange(Nplanets, Nplanets + Nparts)
	#
	# x = np.asarray([sim.particles[int(j)].x for j in bodyrange])
	# y = np.asarray([sim.particles[int(j)].y for j in bodyrange])
	# z = np.asarray([sim.particles[int(j)].z for j in bodyrange])
	#
	#
	#
	# xtilt = np.radians(axtilt) * np.sin(np.radians(axdir + 180))
	# ytilt = np.radians(axtilt) * np.cos(np.radians(axdir + 180))
	#
	# xp = x*np.cos(rot)*np.cos(xtilt) \
	#	 + y*(np.sin(rot)*np.cos(ytilt)-np.cos(rot)*np.sin(xtilt)*np.sin(ytilt)) \
	#	 + z*(np.sin(rot)*np.cos(ytilt)+np.cos(rot)*np.sin(xtilt)*np.cos(ytilt))
	#
	# yp = -x*np.sin(rot)*np.cos(xtilt) \
	#	 + y*(np.cos(rot)*np.cos(ytilt)+np.sin(rot)*np.sin(xtilt)*np.sin(ytilt)) \
	#	 + z*(np.sin(rot)*np.sin(ytilt)+np.cos(rot)*np.sin(xtilt)*np.cos(ytilt))
	#
	# zp = -x*np.sin(xtilt) \
	#	 - y*np.cos(xtilt)*np.sin(ytilt) \
	#	 + z*np.cos(xtilt)*np.cos(ytilt)
	#
	# # Locate lat/lon:
	# lon = np.degrees(np.arctan(yp / xp))
	# lat = np.degrees(np.sin(zp / np.sqrt(xp**2 + yp**2 + zp**2)))
	#
	# # Point on the surface
	# # xsurf = aT * np.cos(np.radians(lon)) * np.cos(np.radians(lat))
	# # ysurf = bT * np.sin(np.radians(lon)) * np.cos(np.radians(lat))
	# # zsurf = cT * np.sin(np.radians(lat))
	#
	# C1 = (a1 ** 2) / (b1 ** 2)
	# C2 = (a1 ** 2) / (c1 ** 2)
	#
	# apos = np.sqrt(xp ** 2 + C1 * yp ** 2 + C2 * zp ** 2)
	#
	# h = np.sqrt(a1**2 - b1**2)
	# k = np.sqrt(a1**2 - c1**2)
	#
	# pot = (sim.G * M) / np.sqrt((apos**2 - h**2) * (apos**2 - k**2))
	#
	# gx = -1*pot * (xp / apos)#np.sqrt(xp**2 + (a1 / b1)**2 * yp**2 + (a1 / c1) * zp**2)
	# gy = -1*pot * (yp / apos)#np.sqrt(xp**2 + (a1 / b1)**2 * yp**2 + (a1 / c1) * zp**2)
	# gz = -1*pot * (zp / apos)#np.sqrt(xp**2 + (a1 / b1)**2 * yp**2 + (a1 / c1) * zp**2)
	#
	# #gx = sim.G * M * (xsurf / apos) * (1 / np.sqrt((apos ** 2 - h ** 2) * (apos ** 2 - k ** 2)))
	# #gy = sim.G * M * (ysurf / apos) * (1 / np.sqrt((apos ** 2 - h ** 2) * (apos ** 2 - k ** 2)))
	# #gz = sim.G * M * (zsurf / apos) * (1 / np.sqrt((apos ** 2 - h ** 2) * (apos ** 2 - k ** 2)))
	#
	# gx += omega ** 2 * xp
	# gy += omega ** 2 * yp
	# #print out rotation phase
	#
	#
	# # Rotate grav vector back to global coordinates
	# gxtilt = np.radians(axtilt) * np.sin(np.radians(axdir))
	# gytilt = np.radians(axtilt) * np.cos(np.radians(axdir))
	#
	# ggx = gx * np.cos(rot + 180.) * np.cos(gxtilt) \
	#	 + gy * (np.sin(rot + 180.) * np.cos(gytilt) - np.cos(rot + 180.) * np.sin(gxtilt) * np.sin(gytilt)) \
	#	 + gz * (np.sin(rot + 180.) * np.cos(gytilt) + np.cos(rot + 180.) * np.sin(gxtilt) * np.cos(gytilt))
	#
	# ggy = -gx * np.sin(rot + 180.) * np.cos(gxtilt) \
	#	 + gy * (np.cos(rot + 180.) * np.cos(gytilt) + np.sin(rot + 180.) * np.sin(gxtilt) * np.sin(gytilt)) \
	#	 + gz * (np.sin(rot + 180.) * np.sin(gytilt) + np.cos(rot + 180.) * np.sin(gxtilt) * np.cos(gytilt))
	#
	# ggz = -gx * np.sin(gxtilt) \
	#	 - gy * np.cos(gxtilt) * np.sin(gytilt) \
	#	 + gz * np.cos(gxtilt) * np.cos(gytilt)
	#
	# #ax += ggx
	# #ay += ggy
	# #az += ggz
	# print('shapeggx', ggx.shape)
	# return ggx, ggy, ggz





def globalvector(sim, tilt, gx, gy, gz, Nplanets, Nparts, binary=False):
	'''In progress
	Inputs the gravitational acceleration vectors in the local frame.
	Outputs gravitational acceleration vectors in global frame.'''

	if binary == True:
		bodies = np.arange(0, Nplanets + Nparts)
		exclude = [1]
		bodyrange = np.delete(bodies, exclude)

		x = np.asarray([sim.particles[int(j)].x for j in bodyrange])
		y = np.asarray([sim.particles[int(j)].y for j in bodyrange])
		z = np.asarray([sim.particles[int(j)].z for j in bodyrange])

		cx = sim.particles[1].x
		cy = sim.particles[1].y
		cz = sim.particles[1].z


	else:
		x = np.asarray([sim.particles[int(j)].x for j in range(Nplanets, Nplanets + Nparts)])
		y = np.asarray([sim.particles[int(j)].y for j in range(Nplanets, Nplanets + Nparts)])
		z = np.asarray([sim.particles[int(j)].z for j in range(Nplanets, Nplanets + Nparts)])

		cx = sim.particles[0].x
		cy = sim.particles[0].y
		cz = sim.particles[0].z

	x -= cx
	y -= cy
	z -= cz

	gmag = np.sqrt(gx ** 2 + gy ** 2 + gz ** 2)

	# adjust tilt
	phi0 = np.degrees(np.arcsin(gmag / gz))
	phi  = phi0 - tilt

	theta = np.degrees(np.arctan(gy / gx))


	# Define globally
	# This is where things need to be fixed. The gravitational vectors relative to the target body will stay the same.
	# Just need to put them back in the locations given in global coordinates
	ggx = gmag * np.cos(np.radians(phi)) * np.sin(np.radians(theta))
	ggy = gmag * np.cos(np.radians(phi)) * np.cos(np.radians(theta))
	ggz = gmag * np.sin(np.radians(phi))

	return ggx, ggy, ggz


def gravity(sim, M, Nplanets, Nparts, atarg, btarg, ctarg, omega, tilt, axdir, binary=False):
	'''In progress
	Inputs gravitational acceleration vectors
	Applies acceleration vectors to system objects'''

	gx, gy, gz = localvector(sim, M, Nplanets, Nparts, atarg, btarg, ctarg, omega, tilt, axdir, binary)

	ggx, ggy, ggz = globalvector(sim, tilt, gx, gy, gz, Nplanets, Nparts, binary)

	return ggx, ggy, ggz





def localvector(sim, M, Nplanets, Nparts, atarg, btarg, ctarg, omega, tilt, axdir, binary=False):
	'''In progress
	Inputs particle positions and ellipsoid axes.
	Outputs gravitational acceleration vector felt by each particle.'''

	h = np.sqrt(atarg**2 - btarg**2)
	k = np.sqrt(atarg**2 - ctarg**2)

	if binary == True:
		bodies = np.arange(0, Nplanets + Nparts)
		exclude = [1]
		bodyrange = np.delete(bodies, exclude)

		xp = np.asarray([sim.particles[int(j)].x for j in bodyrange])
		yp = np.asarray([sim.particles[int(j)].y for j in bodyrange])
		zp = np.asarray([sim.particles[int(j)].z for j in bodyrange])

		cx = sim.particles[1].x
		cy = sim.particles[1].y
		cz = sim.particles[1].z

	else:
		xp = np.asarray([sim.particles[int(j)].x for j in range(Nplanets, Nplanets + Nparts)])
		yp = np.asarray([sim.particles[int(j)].y for j in range(Nplanets, Nplanets + Nparts)])
		zp = np.asarray([sim.particles[int(j)].z for j in range(Nplanets, Nplanets + Nparts)])

		cx = sim.particles[0].x
		cy = sim.particles[0].y
		cz = sim.particles[0].z

	x = xp - cx
	y = yp - cy
	z = zp - cz

	lat = np.arctan(z / np.sqrt(x**2 + y**2))
	lon = np.arctan(y / x)

	x0 = atarg * np.cos(lat) * np.cos(lon)
	y0 = btarg * np.cos(lat) * np.sin(lon)
	z0 = ctarg * np.sin(lat)

	dist = np.sqrt((x - x0) ** 2 + (y - y0) ** 2 + (z - z0) ** 2)

	dx = dist * (np.sin(np.radians(tilt)) * np.cos(np.radians(axdir)) )#- np.cos(np.radians(omega * time)))
	dy = dist * (np.sin(np.radians(tilt)) * np.sin(np.radians(axdir)) )#- np.sin(np.radians(omega * time)))
	dz = dist * (np.cos(np.radians(tilt)))

	xloc = x + dx
	yloc = y + dy
	zloc = z + dz

	C1 = (atarg**2) / (atarg**2 - h**2)
	C2 = (atarg**2) / (atarg**2 - k**2)

	apos = np.sqrt(xloc**2 + C1 * yloc**2 + C2 * zloc**2)


	gx = sim.G * M * (xloc / apos) * (1 / np.sqrt((apos ** 2 - h ** 2) * (apos ** 2 - k ** 2)))
	gy = sim.G * M * (yloc / apos) * (1 / np.sqrt((apos ** 2 - h ** 2) * (apos ** 2 - k ** 2)))
	gz = sim.G * M * (zloc / apos) * (1 / np.sqrt((apos ** 2 - h ** 2) * (apos ** 2 - k ** 2)))

	gx += omega ** 2 * xloc
	gy += omega ** 2 * yloc

	return gx, gy, gz


# Do not use the above functions.



### Binary Component ###
def binary(sim, m, r, a, e=0, i=0, periap=0, ascnode=0, f=0):
	''' ***This is ok!***
	Adds binary component to system (orbits central body)

	Parameters
	----------
	sim	: REBOUND simulation
	m	: float; binary mass
	r	: float; binary radius
	a	: float; distance between binaries
	e	: float; (default=0) binary eccentricity
	i	: float; (default=0) binary inclination
	periap	: float; (default=0) argument of periapsis for binary
	ascnode : float; (default=0) longitude of ascending node for binary
	f	: float; (default=0) true anomaly of binary

	Returns
	-------
	sim : REBOUND simulation
	'''

	sim.add(m=m, r=r, a=a, e=e, inc=i, omega=periap, Omega=ascnode, f=f, primary=sim.particles[0])

	return sim





### Particle Size Distribution ###
def partmass(radii, rho):
	''' ***This is ok!***
	Calculates mass of particles

	Parameters
	----------
	radii : array; particle radii
	rho   : float; particle density

	Returns
	-------
	mass : array; particle masses
	'''

	mass = (4. / 3.) * np.pi * radii ** 3 * rho

	return mass


def partrad(r, counts):
	''' ***This is ok!***
	Creates array of particle radii

	Parameters
	----------
	r      : array; possible radii
	counts : array; number of particles for each diameter
	Nparts :   int; total number of particles

	Returns
	-------
	radii : array; particle radii
	'''

	# if sum(counts) != Nparts:
	#   counts[0] += (Nparts - sum(counts))

	# posrad = D / 2.
	rad = []

	for i in range(len(counts)):
		numparts = int(counts[i])

		for n in range(numparts):
			rad.append(r[i])

	radii = np.asarray(rad)
   # print (len(radii))
	return radii


def sizedist(Nparts, Res, rmin, rmax, p):
	''' ***This is ok!***
	Calculates power law size distribution for particles. Based on O'Brian (2003).

	Parameters
	----------
	Nparts :   int; number of particles
	Res    :   int; resolution of power law
	rmin   : float; minimum particle radius
	rmax   : float; maximum particle radius
	p      : float; power of power law distribution

	Returns
	-------
	r     : array; particle radii
	counts: array; number of particles that have each radius size
	'''

	pow = -p
	C = Nparts * (1 / (rmax ** (1 - pow) - rmin ** (1 - pow)))
	const = 1 - Nparts * ((rmax ** (1 - pow)) / (rmax ** (1 - pow) - rmin ** (1 - pow)))
	r = np.linspace(rmin, rmax, int(Res))
	N = -(C * r ** (1 - pow)) + const

	counts = np.round(N, 0)

	# print (counts/np.sum(counts)) * Nparts
	return r, counts, N



### Radiation Pressure ###
def poyntingrobertson(sim, r, a, rho, L=3.828e26, msun=2e30):
	''' ***Do not use this***
	Calculates radiation pressure on particle sizes

	Parameters
	----------
	sim  : REBOUND simulation
	r    : array; particle radii
	a    : float; distance to sun (km)
	rho  : float; density of particles
	L    : float; defaults to luminosity of the sun
	msun : float; defaults to mass of sun

	Returns
	-------
	acc : array; radiation pressure acceleration for each particle size
	'''

	km   = 1e3  # m in km
	days = 86400.  # seconds in day

	G = sim.G  # 6.67e-11 * km**3 * days**2

	Lum = L * days ** 3 * km ** -2

	c = 3e8 * (days / km)

	fpr = ((Lum * r ** 2) / (4 * c ** 2)) * np.sqrt((G * msun) / a ** 5)

	acc = fpr / ((4. / 3.) * np.pi * r ** 3 * rho)

	return acc

def radforce(sim, Nparts, Nplanets, r, rho, L=3.828e26, msun=2e30):
	''' ***Do not use!***
	Calculates Poynting-Robertson Drag on particle sizes

	Parameters
	----------
	sim	: REBOUND simulation
	Nparts	: int; number of particles
	Nplanets: int; number of planets
	r	: array; particle diameters
	rho	: float; particle density
	L	: float; defaults to luminosity of the sun
	msun	: float; defaults to mass of sun

	Returns
	-------
	acc : array; accelerations of each particle
	'''

	for i in np.arange(Nparts):
		p = sim.particles[int(i + Nplanets)]
		x = -(p.x - sim.particles[1].x)
		y = -(p.y - sim.particles[1].y)
		z = -(p.z - sim.particles[1].z)

	dist = np.sqrt(x ** 2 + y ** 2 + z ** 2)
	theta = np.arctan(y / x)
	phi = np.arccos(z / dist)

	acc = poyntingrobertson(r, dist, rho, L, msun)

	return acc






def solarradpress(sim, Nplanets, possun, rho=2e3):
	'''USE THIS EQUATION'''

	SBconst = 5.67e-8     # Stephen-Boltzmann Constant
	Tsun = 5.778e3	      # Kelvin
	Rsun = 6.95508e8      # meters Sun radius

	tsi = 1.36e3	      # W/m**2
	eps = 1.5
	c = 3e8		      # m/s
	Msun = 2e30	      #kg
	G = sim.G

	ax = []
	ay = []
	az = []

	sunind = Nplanets-1

	sun = sim.particles[int(sunind)]
	sunx = sun.x
	suny = sun.y
	sunz = sun.z


	for p in sim.particles[Nplanets:]:
		dist = np.sqrt((p.x-possun[0])**2 + (p.y-possun[1])**2 + (p.z-possun[2])**2)

		Ksc = 1.366e3 * 1.496e11**2

		amean = (3 * Ksc) / (4 * np.pi * c * rho * dist ** 2 * p.r)

		#I = 3.8e26 / (4 * np.pi * dist ** 2)
		#p = (2 * I) / c
		#F = p * np.pi * p.r ** 2

		#amean = (3 * 3.8e26) / (4 * np.pi * c * rho * dist ** 2 * p.r)
	  #  print(amean)
		# *** Note: Make sure to finalize the equation for solar constant at varying radii.
		# It's some form of inverse square law.
		# integrate to solve. Use (1AU, solar constant) as the point to solve the constant


		#Fs = SBconst * Tsun**4 * (Rsun**2 / dist**2)

		#beta = (3 * eps * tsi * dist**2) / (4 * c * G * Msun * rho * p.r)
		#amean = .1 * ((G * Msun) / (dist**2)) * (beta)
		#amean = ((G * Msun) / (dist**2)) - (3 * eps * tsi) / (4 * c * rho * p.r)
		ax.append(amean * ((p.x-possun[0]) / dist))
		ay.append(amean * ((p.y-possun[1]) / dist))
		az.append(amean * ((p.z-possun[2]) / dist))

	acc = (np.asarray(ax), np.asarray(ay), np.asarray(az))

	return acc



def shadow(sim, possun, Nplanets, Nparts, aview, bview, rho=2e3):
	'''Adding a shadow to the target body
	Sets up condition for adding the solar radiation force.'''


	xp = np.asarray([sim.particles[int(j)].x for j in range(Nplanets, Nplanets + Nparts)])
	yp = np.asarray([sim.particles[int(j)].y for j in range(Nplanets, Nplanets + Nparts)])
	zp = np.asarray([sim.particles[int(j)].z for j in range(Nplanets, Nplanets + Nparts)])

	xs = possun[0]
	ys = possun[1]
	zs = possun[2]

	xc = 0.
	yc = 0.
	zc = 0.

	xy = np.sqrt(xs ** 2 + ys ** 2)

	# z greater than 0
	if zs >= 0 and xs >= 0 and ys >= 0:
		theta = 2 * np.pi - np.arctan(ys / xs)
		phi = (2 * np.pi) - (np.arctan(zs / xs) + (np.pi / 2))

	elif zs >= 0 and xs >= 0 and ys < 0:
		theta = np.arctan(ys / xs)
		phi = (2 * np.pi) - (np.arctan(zs / xs) + (np.pi / 2))

	elif zs >= 0 and xs < 0 and ys >= 0:
		theta = np.pi - np.arctan(ys / xs)
		phi = (2 * np.pi) - (np.arctan(zs / xs) + (np.pi / 2))

	elif zs >= 0 and xs < 0 and ys < 0:
		theta = np.pi + np.arctan(ys / xs)
		phi = (2 * np.pi) - (np.arctan(zs / xs) + (np.pi / 2))

	# z less than 0
	elif zs < 0 and xs >= 0 and ys >= 0:
		theta = 2 * np.pi - np.arctan(ys / xs)
		phi = np.arctan(zs / xs)

	elif zs < 0 and xs >= 0 and ys < 0:
		theta = np.arctan(ys / xs)
		phi = np.arctan(zs / xs)

	elif zs < 0 and xs < 0 and ys >= 0:
		theta = np.pi - np.arctan(ys / xs)
		phi = np.arctan(zs / xs)

	elif zs < 0 and xs < 0 and ys < 0:
		theta = np.pi + np.arctan(ys / xs)
		phi = np.arctan(zs / xs)



	xpprime = xp * np.cos(theta) \
			  - yp * np.sin(theta)
	ypprime = xp * np.sin(theta) * np.cos(phi) \
			  + yp * np.cos(theta) * np.cos(phi) \
			  - zp * np.sin(phi)
	zpprime = xp * np.sin(theta) * np.sin(phi) \
			  + yp * np.cos(theta) * np.sin(phi) \
			  + zp * np.cos(phi)

	acc = solarradpress(sim, Nplanets, possun, rho)
	axall = acc[0]
	ayall = acc[1]
	azall = acc[2]

	ax = np.where((xpprime ** 2 / aview ** 2) + (ypprime ** 2 / bview ** 2) > 1,
					 axall,
					 0)

	ay = np.where((xpprime ** 2 / aview ** 2) + (ypprime ** 2 / bview ** 2) > 1,
					 ayall,
					 0)

	az = np.where((xpprime ** 2 / aview ** 2) + (ypprime ** 2 / bview ** 2) > 1,
					 azall,
					 0)

	accel = (ax, ay, az)

	return accel
