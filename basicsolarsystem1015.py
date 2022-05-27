### Basic Solar System Model ###
import rebound
import numpy as np
import matplotlib.pyplot as plt
import csv
from scipy import interpolate
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



def ejectaconepos(mex, aT, bT, cT, h, lat, lon, beta, Nparts, tpos, mtarg, rhoT, shapemodel=False, vertfile=None):
	'''Setting ejecta cone position'''

	# cone base
	#x0p = aT * np.cos(np.radians(lat)) * np.cos(np.radians(lon))
	#y0p = bT * np.cos(np.radians(lat)) * np.sin(np.radians(lon))
	#z0p = cT * np.sin(np.radians(lat))

	# calculate cone base on surface
	if shapemodel == True:
		x = []
		y = []
		z = []
		# read in data files
		with open(vertfile) as vertices:
			vertice = csv.reader(vertices, delimiter='\t')
			rownum = 0
			for row in vertice:
				#print(float(row[2]))
				#print(row)
				x.append(float(row[0]))
				y.append(float(row[1]))
				z.append(float(row[2]))
				rownum += 1
		vertx = 1e3 * np.asarray(x)
		verty = 1e3 * np.asarray(y)
		vertz = 1e3 * np.asarray(z)
		#print(vertx)

		#vlon = np.degrees(np.arccos(vertx / np.sqrt(vertx ** 2 + verty ** 2)))
		#vlat = np.degrees(np.arcsin(vertz / np.sqrt(vertx ** 2 + verty ** 2 + vertz ** 2)))

		p = np.sqrt(vertx ** 2 + verty ** 2)
		#print('p',p)

		vlat = np.degrees(np.arctan(vertz / p))
		vlon = np.degrees(np.arctan(verty / vertx))   #np.arccos(cx / (np.sqrt(cx ** 2 + cy ** 2 + cz ** 2) * np.cos(lat)))
		#print('vlat',vlat)

		points = (vlon, vlat)
		vr = np.sqrt(vertx ** 2 + verty ** 2 + vertz ** 2)
		point = (lon, lat)
		hav = (.5 * (1 - np.cos(np.radians(lat - vlat))) + .5 * (1 - np.cos(np.radians(lon - vlon))) * np.cos(np.radians(lat)) * np.cos(np.radians(vlat)))
		#inter = interpolate.interp2d(vlon, vlat, vr)
		rbase = vr[int(np.argmin(hav))]	  #inter(lon, lat)
		print('rbase',rbase)

		x0p = rbase * np.cos(np.radians(lat)) * np.cos(np.radians(lon))
		y0p = rbase * np.cos(np.radians(lat)) * np.sin(np.radians(lon))
		z0p = rbase * np.sin(np.radians(lat))
		print(x0p)

	else:
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
	#	 * (((g * a) / vi ** 2) * (rhoT / rhoI) ** (-1. / 3.)
	#		+ (Ybar / (rhoT * vi ** 2)) ** ((2 + mu) / 2)) ** (-mu / (2 + mu))

	#Rg = K1 * (rhoT / mi ) ** (-1./3.) * (rhoT / rhoI) ** ((2 + mu - 6 * 0.4) / (6 + 3 * mu)) * ((g * a) / vi ** 2) ** (-mu / (2 + mu))
	Rg = ((3. / np.pi) * (mex / rhoT)) ** (1. / 3.)
	#print('rg=', Rg)

	#rdisk	= h * np.tan(np.radians(beta))
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
	#phiang	 = opang * np.sin(np.radians(thetaprime))
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



def ejectaconevel(vpos, v0):
	'''Sets up ejecta cone velocities 2/24/20'''

	# impact position on surface
	#x0p = pos0p[0]
	#y0p = pos0p[1]
	#z0p = pos0p[2]

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




def rmdatasave(rmx, rmy, rmz, fname):
	'''USE THIS'''
	data = np.zeros((3, len(rmx)))

	data[0] = rmx
	data[1] = rmy
	data[2] = rmz

	np.savetxt(fname, data, fmt='%.18e', delimiter=' ', newline='\n')







def rmland(sim, Nparts, atarg, btarg, ctarg, landed, inttime, condit, axdir, axtilt, per, timestep, savefile, rotation=False, shapemodel=False, vertfile=None, bincent=False):
	'''USE THIS'''


	bodies = np.linspace(0, sim.N-1, sim.N)

	xp = np.asarray([sim.particles[int(i)].x for i in bodies])
	yp = np.asarray([sim.particles[int(i)].y for i in bodies])
	zp = np.asarray([sim.particles[int(i)].z for i in bodies])
	#print('xp =',xp)

	if bincent == True:
		centx = xp[1]
		centy = yp[1]
		centz = zp[1]
	else:
		centx = xp[0]
		centy = yp[0]
		centz = zp[0]


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

		if bincent == True:
			centx = xtilt[1]
			centy = ytilt[1]
			centz = ztilt[1]
		else:
			centx = xtilt[0]
			centy = ytilt[0]
			centz = ztilt[0]

		# distance from center of body
		cx = xtilt - centx#sim.particles[0].x
		cy = ytilt - centy#sim.particles[0].y
		cz = ztilt - centz#sim.particles[0].z
		#print('center', np.sqrt(sim.particles[0].x**2 + sim.particles[0].y**2 + sim.particles[0].z**2))

	else:
		cx = xp - centx #sim.particles[0].x
		cy = yp - centy #sim.particles[0].y
		cz = zp - centz #sim.particles[0].z


	cx = cx[2:]
	cy = cy[2:]
	cz = cz[2:]
	
	# calculate point on surface
	e2 = (atarg ** 2 - btarg ** 2) / atarg ** 2
	p = np.sqrt(cx ** 2 + cy ** 2)
	if shapemodel == True:
		lat = np.arctan(cz / p)
	else:
		#print('pay attention', np.arctan(cz / ((1 - e2) * p)))
		lat = np.arctan(cz / ((1 - e2) * p))
                
	lon = np.arctan(cy / cx)  
	dist = np.sqrt(cx ** 2 + cy ** 2 + cz ** 2)
	print('dist=', dist)

	surf = []

	if shapemodel == True:
		x = []
		y = []
		z = []
		# read in data files of shape model vertices to get surface 
		with open(vertfile) as vertices:
			vertice = csv.reader(vertices, delimiter='\t')
			rownum = 0
			for row in vertice:
				#print(float(row[2]))
				#print(row)
				x.append(float(row[0]))
				y.append(float(row[1]))
				z.append(float(row[2]))
				rownum += 1
			vertx = 1e3 * np.asarray(x)
			verty = 1e3 * np.asarray(y)
			vertz = 1e3 * np.asarray(z)

			# Calculate longitudes and latitudes of each vertex
			vp = np.sqrt(vertx ** 2 + verty ** 2)
			vlat = np.arctan(vertz / vp)
			vlon = np.arctan(verty / vertx)

			for i in range(len(lat)):
				hav = (.5 * (1 - np.cos(lat[int(i)] - vlat)) + .5 * (1 - np.cos(lon[int(i)] - vlon)) * np.cos(lat[int(i)]) * np.cos(vlat))
				indmin = np.argmin(hav)
				surf.append(np.sqrt(vertx[int(indmin)]**2 + verty[int(indmin)]**2 + vertz[int(indmin)]**2))
					
				
			
			#shapesurf = np.sqrt(vertx ** 2 + verty ** 2 + vertz ** 2)
			#for d in dist:
			#	surf.append(np.min(np.abs(d - shapesurf)))
			surf = np.asarray(surf)

	else:
		N = atarg / np.sqrt(1 - e2 * (np.sin(lat))**2)

		xsurf = N * np.cos(lat) * np.cos(lon)
		ysurf = N * np.cos(lat) * np.sin(lon)
		zsurf = (1 - e2) * N * np.sin(lat)
		#A = ((np.cos(lat)) ** 2 * (np.cos(lon)) ** 2) / atarg ** 2
		#B = ((np.cos(lat)) ** 2 * (np.sin(lon) ** 2)) / btarg ** 2
		#C = ((np.sin(lat)) ** 2) / ctarg ** 2

		surf = np.sqrt(xsurf ** 2 + ysurf ** 2 + zsurf ** 2)  #(A + B + C) ** -1)
		bigdiff = dist - surf
	print('surf=', surf)	# np.where(dist<=surf))




	landx = []
	landy = []
	landz = []
	landbx, landby, landbz = ([], [], [])

	#inds = np.where(dist < np.mean(surf[2:]))[0]
	for i in reversed(np.arange(len(dist))):
		#print(i)
		if i > 2:
			#print('difference', dist[int(i)], surf[int(i)])
			if dist[int(i)] < surf[int(i)]:
                                #print('surf', surf)
                               

				#print('REMOVE ONE')
				landx.append(sim.particles[int(i+2)].x)
				landy.append(sim.particles[int(i+2)].y)
				landz.append(sim.particles[int(i+2)].z)

				sim.remove(int(i+2))
				landed += 1
				Nparts -= 1





	rmdatasave(landx, landy, landz, savefile + str(inttime) + condit + '.txt')


#	if binary == True:
#		 print('binaries are real')
#		#print indsb
#
#		bodies = np.linspace(0, sim.N - 1, sim.N)
#
#		xp = np.asarray([sim.particles[int(i)].x for i in bodies])
#		yp = np.asarray([sim.particles[int(i)].y for i in bodies])
#		zp = np.asarray([sim.particles[int(i)].z for i in bodies])
#
#		bsurf = abin
#		bx = xp - xp[1]	 # sim.particles[1].x
#		by = yp - yp[1]	 # sim.particles[1].y
#		bz = zp - zp[1]	 # sim.particles[1].z
#		bdist = np.sqrt(bx ** 2 + by ** 2 + bz ** 2)
#		print('binary=', bdist)
#		# np.concatenate((np.where(dist <= surf)[0], np.where(bdist <= bsurf)[0]))
#
 #		 for ib in reversed(np.arange(len(bdist))):
  #			 if ib > 2:
   #				 if bdist[int(ib)] < bsurf[int(ib)]:
    #					 
#					landbx.append(sim.particles[int(ib)].x)
#					landby.append(sim.particles[int(ib)].y)
#					landbz.append(sim.particles[int(ib)].z)
#
#					sim.remove(int(ib))
#
#					landed += 1
#					Nparts -= 1
#
#		rmdatasave(landbx, landby, landbz, 'rmblandpartdata' + str(inttime) + condit + '.txt')

	return sim, landed, Nparts





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



def veldist(sim, r, mtarg, rtarg, mex, mu, rhoT, Cvps=0., Ctg=0.8, Ybar=0.): #vi, a, mi, rhoT, rhoI, mu, K1, Cvps=0., Ctg=0.8, Ybar=0.):
	''' Defines velocity distribution based on Richardson (2007)'''

	G = sim.G

	g = (G * mtarg) / rtarg ** 2

	# transient crater radius

	# Rg = ((3. / np.pi) * (K1 * (mi / rhoT)) * (((g * a) / vi ** 2) * (rhoT / rhoI) ** (-1./3.)
	#										   + (Ybar / (rhoT * vi ** 2)) ** ((2 + mu) / 2)) ** (-3 * mu / (2 + mu))) ** (1. / 3.)

#	Rg = K1 * (rhoT / mi) ** (-1. / 3.) * (rhoT / rhoI) ** ((2 + mu - 6 * 0.4) / (6 + 3 * mu)) * (
#																								 (g * a) / vi ** 2) ** (
#																								 -mu / (2 + mu))

	Rg = ((3. / np.pi) * (mex / rhoT)) ** (1. / 3.)

	#print(Rg)

	# velocity distribution (y-axis)

	Cvpg = (np.sqrt(2) / Ctg) * (mu / (mu + 1))

	ve = Cvpg * np.sqrt(g * Rg) * (r / Rg) ** (-1 / mu)


	veff = np.sqrt(ve ** 2 - Cvpg ** 2 * (g * r) - Cvps ** 2 * (Ybar / rhoT))




	# Gaussian Distribution
	vmean = np.mean(veff)
	vstd  = np.std(veff)

	vnorm = np.random.normal(vmean, vstd, len(r))


	data = np.zeros((2, len(r)))

	data[0, :] = r
	data[1, :] = vnorm#veff
	np.savetxt('veldistinit.txt', data, fmt='%.18e', delimiter=' ', newline='\n')

	return veff
