import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import Ellipse
from integratorparams import *
from bodyparams import *
import additionaleffects as ae

inttime = np.arange(0., 1e5)

inner  = -1.55
innery = -1.15
outer  = 1.55
outery = 1.15

bodycolor = 'steelblue'
partcolor = 'black'

units = 1e3

scale = 1

partsize = .5*np.pi
sun = 0	 # 1 for no binary; 2 for binary

rmx = np.zeros(2)
rmy = np.zeros(2)
rmz = np.zeros(2)
rmbx = np.zeros(2)
rmby = np.zeros(2)
rmbz = np.zeros(2)
binx = []
biny = []
binz = []
blastx, blasty, blastz = (0, 0, 0)
rmballx = []
rmbally = []
rmballz = []
rmt = []
#fig, ax = plt.subplots(1)


for t in inttime:
	omega = (2*np.pi) / perT
	seconds = t * dt
	print('timestep =', t)


	

	fname = 'particledata' + str(int(t)) + condit + '.txt'
	rmfname = 'rmlandpartdata' + str(int(t)) + condit + '.txt'
	if binary == True:
		rmbfname = 'rmblandpartdata' + str(int(t)) + condit + '.txt'
		rmbdata = np.loadtxt(rmbfname)

	data   = np.loadtxt(fname)
	rmdata = np.loadtxt(rmfname)

#	planetorbs = np.loadtxt(planets, skiprows=3, delimiter='	')
	if binary == True:
		Nplanets = 2
	elif binary != True:
		Nplanets = 1


	rdist = np.sqrt((aT**2 + bT**2 + cT**2) * 9000.)

	datax0	= data[0, 0]
	datay0	= data[1, 0]
	dataz0	= data[2, 0]
	datavx0 = data[3, 0]
	datavy0 = data[4, 0]
	datavz0 = data[5, 0]



	data[0, :] -= data[0, 0]
	data[1, :] -= data[1, 0]
	data[2, :] -= data[2, 0]
	data[3, :] -= data[3, 0]
	data[4, :] -= data[4, 0]
	data[5, :] -= data[5, 0]


#	print('dist =', np.sqrt(data[0,5]**2 + data[1,5]**2 + data[2,5]**2))



	ox = data[0, 0] / units
	oy = data[1, 0] / units
	oz = data[2, 0] / units



	if data[0, Nplanets-1:].shape != ():
		x = data[0, Nplanets:]
		y = data[1, Nplanets:]
		z = data[2, Nplanets:]

	if t == 0:
		startx = np.mean(x)
		starty = np.mean(y)


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

	if binary == True:
		if len(rmbdata) != 0:
			rmbdata[0] -= datax0
			rmbdata[1] -= datay0
			rmbdata[2] -= dataz0

			rmbdatax = np.array(rmbdata[0])
			rmbdatay = np.array(rmbdata[1])
			rmbdataz = np.array(rmbdata[2])

			if rmbdatax.shape == ():
				rmbx = np.zeros(2)
			elif rmbdatax.shape !=():
				rmbx = np.concatenate((np.zeros(2), rmbdatax), axis=0)
				#rmbdatax = np.array([rmbdatax])

			if rmbdatay.shape == ():
				rmby = np.zeros(2)
			elif rmbdatay.shape != ():
				rmby = np.concatenate((np.zeros(2), rmbdatay), axis=0)
				#rmbdatay = np.array([rmbdatay])

			if rmbdataz.shape == ():
				rmbz = np.zeros(2)
			elif rmbdataz.shape != ():
				rmbz = np.concatenate((np.zeros(2), rmbdataz), axis=0)
				#rmbdataz = np.array([rmbdataz])

			#rmbx = np.concatenate((rmbx, rmbdatax), axis=0)
			#rmby = np.concatenate((rmby, rmbdatay), axis=0)
			#rmbz = np.concatenate((rmbz, rmbdataz), axis=0)
			print('rmbx',rmbx)

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
	dx = -possun[0] / np.sqrt(possun[0]**2 + possun[1]**2 + possun[2]**2)  #-data[0, sun] / np.sqrt(data[0, sun]**2 + data[1, sun]**2 + data[2, sun]**2)
	dy = -possun[1] / np.sqrt(possun[0]**2 + possun[1]**2 + possun[2]**2)  #-data[1, sun] / np.sqrt(data[0, sun]**2 + data[1, sun]**2 + data[2, sun]**2)
	dz = -possun[2] / np.sqrt(possun[0]**2 + possun[1]**2 + possun[2]**2)  #-data[2, sun] / np.sqrt(data[0, sun]**2 + data[1, sun]**2 + data[2, sun]**2)

	bdx = data[0, 1] / np.sqrt(data[0, 1]**2 + data[1, 1]**2 + data[2, 1]**2)
	bdy = data[1, 1] / np.sqrt(data[0, 1]**2 + data[1, 1]**2 + data[2, 1]**2)
	bdz = data[2, 1] / np.sqrt(data[0, 1]**2 + data[1, 1]**2 + data[2, 1]**2)


	print(data[0, 0])
# Plot xy
	fig, ax = plt.subplots()

	fig.subplots_adjust(hspace=.4, bottom=.25, left=.25)

	if binary == True:
		bx = data[0, 1]/units
		by = data[1, 1]/units
		bz = data[2, 1]/units

		print('bin', np.sqrt(bx**2 + by**2 + bz**2))
		binx.append(bx)
		biny.append(by)
		circle3 = plt.Circle((bx, by), (rbin/units)*scale, alpha=.75, hatch="/")
		ax.add_artist(circle3)

		orbit1 = plt.Circle((ox, oy), abin / units, facecolor='None', edgecolor='black', linestyle='--')
		ax.add_artist(orbit1)
		#ax.arrow(bx, by, .1 * (by / np.sqrt(bx**2 + by**2)), .1 * (bx / np.sqrt(bx**2 + by**2)), head_width=.05, color='black', label='Binary')



	circle1 = plt.Circle((ox, oy), (rtarg/units)*scale, facecolor=bodycolor, alpha=.60)
	ax.add_artist(circle1)	# scatter(0, 0, s=rtarg ** 2 * 9000, c='blue', alpha=0.75)

	ax.scatter(x/units, y/units, s=partsize, c=partcolor, edgecolors=None)

	#if len(rmbx) != 0:
	if binary == True:
		if len(rmbdata) != 0:
			rmbinx = rmbx[2:]/units
			rmbiny = rmby[2:]/units
			rmbinz = rmbz[2:]/units

			rmballx = list(rmballx)
			rmbally = list(rmbally)
			rmballz = list(rmballz)

			for x in rmbinx:
				rmballx.append(x)
			for y in rmbiny:
				rmbally.append(y)
			for z in rmbinz:
				rmballz.append(z)
			rmt.append(t)

		if len(rmt) != 0:
			bdifx = bx - blastx
			bdify = by - blasty
			bdifz = bz - blastz

			rmballx += bdifx
			rmbally += bdify
			rmballz += bdifz

			ax.scatter(rmballx[2:], rmbally[2:], s=partsize, c='green', edgecolors=None)

			#
			# rmdist = np.sqrt((rmx[2:]/units)**2
			# 		 +(rmy[2:]/units)**2
			# 		 +(rmz[2:]/units)**2)
			#
			# bindist = np.sqrt((rmx[2:]/units - bx) ** 2
			# 				  + (rmy[2:]/units - by) ** 2
			# 				  + (rmz[2:]/units - bz) ** 2)
			#
			# for rmd in np.arange(len(rmdist)):
			# 	if bindist[rmd] <= aB:
			# 		bdifx = bx - blastx
			# 		bdify = by - blasty
			# 		bdifz = bz - blastz
			#
			# 		rmbx = rmx[rmd+2]/units + bdifx
			# 		rmby = rmy[rmd+2]/units + bdify
			# 		rmbz = rmz[rmd+2]/units + bdifz
			#
			# 		ax.scatter(rmbx, rmby, s=partsize, c='green', edgecolors=None)
			#
			# 	else:
			# 		ax.scatter(rmx[rmd + 2] / units, rmy[rmd + 2] / units, s=partsize, c='red', edgecolors=None)


	ax.scatter(rmx[2:] / units, rmy[2:] / units, s=partsize, c='red', edgecolors=None)
			# rmdist = np.sqrt((rmx[2:] / units) ** 2
			# 				 + (rmy[2:] / units) ** 2
			# 				 + (rmz[2:] / units) ** 2)
			# for rmd in np.arange(len(rmdist)):
			# 	if rmdist[rmd] <= aT:
			# 		ax.scatter(rmx[rmd + 2] / units, rmy[rmd + 2] / units, s=partsize, c='red', edgecolors=None)


	if radpress == True:
		ax.arrow(-1.2, -.75, .050*dx, .050*dy, head_width=.07, color='black', label='Sun')
		ax.annotate('Sun', xy=(-1.4, -1), fontsize=16)

	ax.annotate('t=' + str(np.round(seconds, 3)) + 's', xy=((inner+.05), (outery-.2)), fontsize=20)

	#ax.arrow(1, -.2, .1*bdx, .1*bdy, head_width=.05, color='black')

	ax.set_xlim(inner, outer)
	ax.set_ylim(innery, outery)
	axissize = 18
	ax.tick_params(labelsize=axissize)
	font = 20
	ax.set_title('Case 2', fontsize=font)
	ax.set_xlabel('x [km]', fontsize=font)
	ax.set_ylabel('y [km]', fontsize=font)

	plt.savefig('xy-zin' + str(int(t)) + condit + '.png')
	plt.close(fig)



# plot xz
	fig, ax = plt.subplots()

	fig.subplots_adjust(hspace=.4, bottom=.25, left=.25)

	if binary == True:
		bx = data[0, 1]/units
		bz = data[2, 1]/units
		binz.append(bz)
		circle4 = plt.Circle((bx, bz), (rbin/units)*scale, alpha=.65, hatch="/")
		ax.add_artist(circle4)

	circle2 = plt.Circle((ox, oz), (rtarg/units)*scale, facecolor=bodycolor, alpha=.60)
	ax.add_artist(circle2)#scatter(0, 0, s=rtarg ** 2 * 9000, c='blue', alpha=0.75)

	ax.scatter(x/units, z/units, s=partsize, c=partcolor, edgecolors=None)

	if len(rmx) != 0:
		if binary == True:
			if len(rmt) != 0:
				ax.scatter(rmballx, rmballz, s=partsize, c='green', edgecolors=None)

		# 	rmdist = np.sqrt((rmx[2:] / units) ** 2
		# 					 + (rmy[2:] / units) ** 2
		# 					 + (rmz[2:] / units) ** 2)
		#
		# 	bindist = np.sqrt((rmx[2:] / units - bx) ** 2
		# 					  + (rmy[2:] / units - by) ** 2
		# 					  + (rmz[2:] / units - bz) ** 2)
		#
		# 	for rmd in np.arange(len(rmdist)):
		# 		if rmdist[rmd] <= aT:
		# 			ax.scatter(rmx[rmd + 2] / units, rmy[rmd + 2] / units, s=partsize, c='red', edgecolors=None)
		# 		if bindist[rmd] <= aB:
		# 			bdifx = bx - blastx
		# 			bdify = by - blasty
		# 			bdifz = bz - blastz
		#
		# 			rmbx = rmx[rmd + 2] / units + bdifx
		# 			rmby = rmy[rmd + 2] / units + bdify
		# 			rmbz = rmz[rmd + 2] / units + bdifz
		#
		# 			ax.scatter(rmbx, rmby, s=partsize, c='red', edgecolors=None)
		#
		# else:
		# 	rmdist = np.sqrt((rmx[2:] / units) ** 2
		# 					 + (rmy[2:] / units) ** 2
		# 					 + (rmz[2:] / units) ** 2)
		# 	for rmd in np.arange(len(rmdist)):
		# 		if rmdist[rmd] <= aT:
		ax.scatter(rmx[2:] / units, rmz[2:] / units, s=partsize, c='red', edgecolors=None)
		# rmdist = np.sqrt((rmx[2:]/units)**2
		# 		 +(rmy[2:]/units)**2
		# 		 +(rmz[2:]/units)**2)
		# for rmd in np.arange(len(rmdist)):
		# 	if rmdist[rmd] <= rtarg:
		# 		ax.scatter(rmx[rmd+2]/units, rmz[rmd+2]/units, s=partsize, c='red', edgecolors=None)



	if radpress == True:
		ax.arrow(-1.2, -.75, .050*dx, .050*dz, head_width=.07, color='black', label='Sun')
		ax.annotate('Sun', xy=(-1.4, -1), fontsize=16)

	ax.annotate('t=' + str(np.round(seconds, 3)) + 's', xy=((inner+.05), (outery-.2)), fontsize=20)

	ax.set_xlim(inner, outer)
	ax.set_ylim(innery, outery)
	ax.tick_params(labelsize=axissize)
	#ax.set_title('t=' + str(np.round(seconds, 3)) + 's', fontsize=font)
	ax.set_xlabel('x [km]', fontsize=font)
	ax.set_ylabel('z [km]', fontsize=font)

	plt.savefig('xz-zin' + str(int(t)) + condit + '.png')
	#print (np.sqrt(np.asarray(binx)**2 + np.asarray(biny)**2 + np.asarray(binz)**2))
	plt.close(fig)



        
#plot yz
	fig, ax = plt.subplots()

	fig.subplots_adjust(hspace=.4, bottom=.25, left=.25)

	if binary == True:
		by = data[1, 1] / units
		bz = data[2, 1] / units
		binz.append(bz)
		circle4 = plt.Circle((by, bz), (rbin / units) * scale, alpha=.65, hatch="/")
		ax.add_artist(circle4)

	circle2 = plt.Circle((oy, oz), (rtarg / units) * scale, facecolor=bodycolor, alpha=.60)
	ax.add_artist(circle2)  # scatter(0, 0, s=rtarg ** 2 * 9000, c='blue', alpha=0.75)

	ax.scatter(y / units, z / units, s=partsize, c=partcolor, edgecolors=None)

	if len(rmx) != 0:
		if binary == True:
			if len(rmt) != 0:
				ax.scatter(rmbally, rmballz, s=partsize, c='green', edgecolors=None)

		# rmdist = np.sqrt((rmx[2:] / units) ** 2
		# 					 + (rmy[2:] / units) ** 2
		# 					 + (rmz[2:] / units) ** 2)
		#
		# 	bindist = np.sqrt((rmx[2:] / units - bx) ** 2
		# 					  + (rmy[2:] / units - by) ** 2
		# 					  + (rmz[2:] / units - bz) ** 2)
		#
		# 	for rmd in np.arange(len(rmdist)):
		# 		if rmdist[rmd] <= aT:
		# 			ax.scatter(rmx[rmd + 2] / units, rmy[rmd + 2] / units, s=partsize, c='red', edgecolors=None)
		# 		if bindist[rmd] <= aB:
		# 			bdifx = bx - blastx
		# 			bdify = by - blasty
		# 			bdifz = bz - blastz
		#
		# 			rmbx = rmx[rmd + 2] / units + bdifx
		# 			rmby = rmy[rmd + 2] / units + bdify
		# 			rmbz = rmz[rmd + 2] / units + bdifz
		#
		# 			ax.scatter(rmbx, rmby, s=partsize, c='red', edgecolors=None)
		#
		# else:
		# 	rmdist = np.sqrt((rmx[2:] / units) ** 2
		# 					 + (rmy[2:] / units) ** 2
		# 					 + (rmz[2:] / units) ** 2)
		# 	for rmd in np.arange(len(rmdist)):
		# 		if rmdist[rmd] <= aT:
		ax.scatter(rmy[2:] / units, rmz[2:] / units, s=partsize, c='red', edgecolors=None)
	# rmdist = np.sqrt((rmx[2:]/units)**2
	# 		 +(rmy[2:]/units)**2
	# 		 +(rmz[2:]/units)**2)
	# for rmd in np.arange(len(rmdist)):
	# 	if rmdist[rmd] <= rtarg:
	# 		ax.scatter(rmx[rmd+2]/units, rmz[rmd+2]/units, s=partsize, c='red', edgecolors=None)



	if radpress == True:
		ax.arrow(-1.2, -.75, .050 * dx, .050 * dz, head_width=.07, color='black', label='Sun')
		ax.annotate('Sun', xy=(-1.4, -1), fontsize=16)

	ax.annotate('t=' + str(np.round(seconds, 3)) + 's', xy=((inner + .05), (outery - .2)), fontsize=20)

	ax.set_xlim(inner, outer)
	ax.set_ylim(innery, outery)
	ax.tick_params(labelsize=axissize)
	# ax.set_title('t=' + str(np.round(seconds, 3)) + 's', fontsize=font)
	ax.set_xlabel('y [km]', fontsize=font)
	ax.set_ylabel('z [km]', fontsize=font)

	plt.savefig('yz-zin' + str(int(t)) + condit + '.png')
	# print (np.sqrt(np.asarray(binx)**2 + np.asarray(biny)**2 + np.asarray(binz)**2))
	plt.close(fig)



	

	# Plot 3D
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')

	fig.subplots_adjust(hspace=.4, bottom=.25, left=.25)

	u = np.linspace(0, 2*np.pi, 10)
	v = np.linspace(-(np.pi/2), (np.pi/2), 10)
	xc = (aT/units) * np.outer(np.cos(u), np.cos(v))
	yc = (bT/units) * np.outer(np.sin(u), np.cos(v))
	zc = (cT/units) * np.outer(np.ones(np.size(u)), np.sin(v))

	angles = np.linspace(0, 2*np.pi, 360)

	xeq = rtarg * np.cos(angles)
	yeq = rtarg * np.sin(angles)
	zeq = np.zeros(360)

	ax.plot(xeq/units,
		yeq/units,
		zeq/units,
		color='black',
		linewidth=.5)

	ax.plot_surface(xc, yc, zc,
			color=bodycolor,
			alpha=.2)

	ax.plot_wireframe(xc,yc,zc,
			  linewidths=.2,
			  edgecolor='black')

	axis = np.linspace(-(aT+50), (aT+50), 20)
	zeros = np.linspace(0, 0, 20)
	ax.plot(zeros,
		zeros,
		axis/units,
		color='black',
		linewidth=0.4)


		
	impact = (rtarg * np.cos((omega * seconds)+np.radians(lon)) * np.cos(np.radians(lat)) + 1,
		  rtarg * np.sin((omega * seconds)+np.radians(lon)) * np.cos(np.radians(lat)) + 1,
		  rtarg * np.sin(np.radians(lat)) + 1)
	ax.scatter(impact[0]/units,
		   impact[1]/units,
		   impact[2]/units,
		   color='black',
		   marker="x")
	

	if len(rmx) != 0:
		ax.scatter(rmx[2:]/units,
			   rmy[2:]/units,
			   rmz[2:]/units,
			   s=partsize,
			   c='red',
			   facecolor=None)
	

	ax.scatter(x/units,
		   y/units,
		   z/units,
		   s=partsize,
		   c=partcolor,
		   edgecolors=None)

	ax.set_title('t=' + str(np.round(seconds, 3)) + 's')

	ax.set_xlabel('x [km]')
	ax.set_ylabel('y [km]')
	ax.set_zlabel('z [km]')

	ax.set_xlim(-.2, -.2)
	ax.set_ylim(-.2, .2)
	ax.set_zlim(-.2, .2)

	plt.savefig('3d' + str(int(t)) + condit + '.png')
	plt.close(fig)

	if binary == True:
		blastx = bx
		blasty = by
		blastz = bz
