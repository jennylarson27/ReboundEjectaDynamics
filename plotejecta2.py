import numpy as np
import matplotlib.pyplot as plt
from integratorparams import *
from bodyparams import *
import additionaleffects as ae

inttime = np.arange(1e5)

inner  = -1.3
innery = -1.3
outer  = 1.3
outery = 1.3

units = 1e3

scale = 1e0

partsize = .5*np.pi
sun = 2

rmx = np.zeros(2)
rmy = np.zeros(2)
rmz = np.zeros(2)
binx = []
biny = []
binz = []
#fig, ax = plt.subplots(1)


for t in inttime:

	seconds = t * dt

	fname = 'particledata' + str(int(t)) + condit + '.txt'
	rmfname = 'rmlandpartdata' + str(int(t)) + condit + '.txt'

	data   = np.loadtxt(fname)
	rmdata = np.loadtxt(rmfname)

#	planetorbs = np.loadtxt(planets, skiprows=3, delimiter='	')
	if binary == True:
		Nplanets = 3
	elif binary != True:
		Nplanets = 2


	rdist = np.sqrt(rtarg**2 * 9000.)

	datax0 = data[0, 0]
	datay0 = data[1, 0]
	dataz0 = data[2, 0]
	datavx0 = data[3, 0]
	datavy0 = data[4, 0]
	datavz0 = data[5, 0]

	data[0, :] -= data[0, 0]
	data[1, :] -= data[1, 0]
	data[2, :] -= data[2, 0]
	data[3, :] -= data[3, 0]
	data[4, :] -= data[4, 0]
	data[5, :] -= data[5, 0]



	ox = data[0, 0] / units
	oy = data[1, 0] / units
	oz = data[2, 0] / units



	if data[0, Nplanets-2:].shape != ():
		x = data[0, Nplanets:]
		y = data[1, Nplanets:]
		z = data[2, Nplanets:]


	if len(rmdata) != 0:
		rmdata[0] -= datax0
		rmdata[1] -= datay0
		rmdata[2] -= dataz0

		rmdatax = np.array(rmdata[0])
		rmdatay = np.array(rmdata[1])
		rmdataz = np.array(rmdata[2])

		if rmdatax.shape == ():
			rmdatax = np.array([rmdatax])
		if rmdatay.shape == ():
			rmdatay = np.array([rmdatay])
		if rmdataz.shape == ():
			rmdataz = np.array([rmdataz])

		rmx = np.concatenate((rmx, rmdatax), axis=0)
		rmy = np.concatenate((rmy, rmdatay), axis=0)
		rmz = np.concatenate((rmz, rmdataz), axis=0)
	print ('removed =', len(rmx))

	ux = x / np.sqrt(x ** 2 + y ** 2 + z ** 2)
	uy = y / np.sqrt(x ** 2 + y ** 2 + z ** 2)
	uz = z / np.sqrt(x ** 2 + y ** 2 + z ** 2)


	vx = data[3, Nplanets:]
	vy = data[4, Nplanets:]
	vz = data[5, Nplanets:]

#	if sizedist == True:
		# set up size distribution
#		r, counts, N = ae.sizedist(Nparts, Res, rmin, rmax, p)
#		radii = ae.partrad(r, counts, Nparts)

	# arrow towards sun
	dx = data[0, sun] / np.sqrt(data[0, sun]**2 + data[1, sun]**2 + data[2, sun]**2)
	dy = data[1, sun] / np.sqrt(data[0, sun]**2 + data[1, sun]**2 + data[2, sun]**2)
	dz = data[2, sun] / np.sqrt(data[0, sun]**2 + data[1, sun]**2 + data[2, sun]**2)


	print data[0, 0]
# Plot xy
	fig, ax = plt.subplots()

	fig.subplots_adjust(hspace=.4, bottom=.25, left=.25)

	if binary == True:
		bx = data[0, 1]/units
		by = data[1, 1]/units
		binx.append(bx)
		biny.append(by)
		circle3 = plt.Circle((bx, by), (rbin/units)*scale, alpha=.75, color='green')
		ax.add_artist(circle3)

		orbit1 = plt.Circle((ox, oy), abin / units, facecolor='None', edgecolor='green', linestyle='--')
		ax.add_artist(orbit1)
		#ax.arrow(bx, by, .1 * (by / np.sqrt(bx**2 + by**2)), .1 * (bx / np.sqrt(bx**2 + by**2)), head_width=.05, color='black', label='Binary')



	circle1 = plt.Circle((ox, oy), (rtarg/units)*scale, alpha=.65)
	ax.add_artist(circle1)  # scatter(0, 0, s=rtarg ** 2 * 9000, c='blue', alpha=0.75)

	ax.scatter(x/units, y/units, s=partsize, c='red', edgecolors='none')

	if len(rmx) != 0:
		ax.scatter(rmx[2:]/units, rmy[2:]/units, s=partsize, c='black', edgecolors='none')

	if radpress == True:
		ax.arrow(-100, 150, -50*dx, -50*dy, head_width=5, color='black', label='Sun')
		ax.annotate('Sun', xy=(-100, 175), fontsize=16)

	ax.set_xlim(inner, outer)
	ax.set_ylim(innery, outery)
	axissize = 18
	ax.tick_params(labelsize=axissize)
	font = 20
	ax.set_title('t=' + str(np.round(seconds, 3)) + 's', fontsize=font)
	ax.set_xlabel('x (km)', fontsize=font)
	ax.set_ylabel('y (km)', fontsize=font)

	plt.savefig('xy-zin' + str(int(t)) + condit + '.png')



# plot xz
	fig, ax = plt.subplots()

	fig.subplots_adjust(hspace=.4, bottom=.25, left=.25)

	if binary == True:
		bx = data[0, 1]/units
		bz = data[2, 1]/units
		binz.append(bz)
		circle4 = plt.Circle((bx, bz), (rbin/units)*scale, alpha=.65, color='green')
		ax.add_artist(circle4)

	circle2 = plt.Circle((ox, oz), (rtarg/units)*scale, alpha=.65)
	ax.add_artist(circle2)#scatter(0, 0, s=rtarg ** 2 * 9000, c='blue', alpha=0.75)

	ax.scatter(x/units, z/units, s=partsize, c='red', edgecolors='none')

	if len(rmx) != 0:
		ax.scatter(rmx[2:]/units, rmz[2:]/units, s=partsize, c='black', edgecolors='none')

	if radpress == True:
		ax.arrow(-100, 150, -50*dx, -50*dz, head_width=5, color='black', label='Sun')
		ax.annotate('Sun', xy=(-100, 175), fontsize=16)

	ax.set_xlim(inner, outer)
	ax.set_ylim(innery, outer)
	ax.tick_params(labelsize=axissize)
	ax.set_title('t=' + str(np.round(seconds, 3)) + 's', fontsize=font)
	ax.set_xlabel('x (km)', fontsize=font)
	ax.set_ylabel('z (km)', fontsize=font)

	plt.savefig('xz-zin' + str(int(t)) + condit + '.png')
	#print (np.sqrt(np.asarray(binx)**2 + np.asarray(biny)**2 + np.asarray(binz)**2))
