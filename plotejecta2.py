import numpy as np
import matplotlib.pyplot as plt
from integratorparams import *
from bodyparams import *
import additionaleffects as ae

inttime = np.asarray([3, 6, 10, 16, 21, 26])

fig, ax = plt.subplots(12)

for t in inttime:

	seconds = t * dt * days

	fname = 'particledata' + str(int(t)) + condit + '.txt'

	data = np.loadtxt(fname)

	Nplanets = len(planets) + 2


	rdist = np.sqrt(rtarg**2 * 9000.)

	x = data[0, Nplanets:]
	y = data[1, Nplanets:]
	z = data[2, Nplanets:]

	ux = x / np.sqrt(x ** 2 + y ** 2 + z ** 2)
	uy = y / np.sqrt(x ** 2 + y ** 2 + z ** 2)
	uz = z / np.sqrt(x ** 2 + y ** 2 + z ** 2)

	#x += ux * rdist
	#y += uy * rdist
	#z += uz * rdist

	vx = data[3, Nplanets:]
	vy = data[4, Nplanets:]
	vz = data[5, Nplanets:]

	if sizedist == True:
		# set up size distribution
		r, counts, N = ae.sizedist(Nparts, Res, rmin, rmax, p)
		radii = ae.partrad(r, counts, Nparts)

	# arrow towards sun
	dx = data[0, 1] / np.sqrt(data[0, 1]**2 + data[1, 1]**2 + data[2, 1]**2)
	dy = data[1, 1] / np.sqrt(data[0, 1]**2 + data[1, 1]**2 + data[2, 1]**2)
	dz = data[2, 1] / np.sqrt(data[0, 1]**2 + data[1, 1]**2 + data[2, 1]**2)





	fig.subplots_adjust(hspace=.4)

	if binary == True:
		bx = data[0, 2]
		by = data[1, 2]
		bz = data[2, 2]
		circle3 = plt.Circle((bx, by), rbin, alpha=.65)
		ax[0].add_artist(circle3)

	circle1 = plt.Circle((0, 0), rtarg, alpha=.65)
	ax[0].add_artist(circle1)  # scatter(0, 0, s=rtarg ** 2 * 9000, c='blue', alpha=0.75)
	ax[0].scatter(x, y, s=.75*np.pi * 1e0**2, c='red', edgecolors='none')
	ax[0].arrow(3.5, .25, .5*dx, .5*dy, head_width=.05, color='black', label='Sun')
	ax[0].annotate('Sun', xy=(3.5, 0), fontsize=16)
	ax[0].set_xlim(0, 4)
	ax[0].set_ylim(-2., .5)
	axissize = 18
	ax[0].tick_params(labelsize=axissize)
	font = 20
	ax[0].set_title('t=' + str(np.round(seconds, 3)) + 's', fontsize=font)
	ax[0].set_xlabel('x (km)', fontsize=font)
	ax[0].set_ylabel('y (km)', fontsize=font)

	plt.savefig('xydidymos-' + str(int(t)) + condit + '.png')


	fig, ax = plt.subplots(2)

	fig.subplots_adjust(hspace=.4)

	if binary == True:
		circle4 = plt.Circle((bx, bz), rbin, alpha=.65)
		ax[0].add_artist(circle4)

	circle2 = plt.Circle((0, 0), rtarg, alpha=.65)
	ax[0].add_artist(circle2)#scatter(0, 0, s=rtarg ** 2 * 9000, c='blue', alpha=0.75)
	ax[0].scatter(x, z, s=.75*np.pi * 1e0**2, c='red', edgecolors='none')
	ax[0].arrow(2.75, .25, .5*dx, .5*dz, head_width=.05, color='black', label='Sun')
	ax[0].annotate('Sun', xy=(2.5, 0.3), fontsize=16)
	ax[0].set_xlim(0, 3)
	ax[0].set_ylim(-.5, .5)
	ax[0].tick_params(labelsize=axissize)
	ax[0].set_title('t=' + str(np.round(seconds, 3)) + 's', fontsize=font)
	ax[0].set_xlabel('x (km)', fontsize=font)
	ax[0].set_ylabel('z (km)', fontsize=font)

	plt.savefig('xzdidymos-' + str(int(t)) + condit + '.png')