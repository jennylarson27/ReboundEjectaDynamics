### Package of additional effects such as size distribution, radiation pressure, binaries, etc ###

import rebound
import numpy as np
import math as ma
import basicsolarsystem1015 as bss



### Ellipsoidal Gravitational Acceleration ###

def spheregrav(sim, Mtarg, Nparts, binary=False):
    ''' ***Might not need this afterall*
    Calculates gravity of spherical body

    Currently in progress
    Oct 4, 2019 - testing
    '''

    G = sim.G

    inds = np.linspace(0, Nparts-1, Nparts)

    if binary == True:
        gravbod = 1
        bodies  = np.delete(inds, gravbod, 0)
    else:
        gravbod = 0
        bodies  = np.delete(inds, gravbod, 0)


    c = sim.particles[gravbod]

    # particle positions relative to gravitational body
    xp = np.asarray([sim.particles[int(i)].x for i in bodies]) - c.x
    yp = np.asarray([sim.particles[int(i)].y for i in bodies]) - c.y
    zp = np.asarray([sim.particles[int(i)].z for i in bodies]) - c.z

    r = np.sqrt(xp**2 + yp**2 + zp**2)

    amag = (G * Mtarg) / r**2

    ax = amag * (xp / r)
    ay = amag * (yp / r)
    az = amag * (zp / r)

    return ax, ay, az


def ellipgrav(sim, Mtarg, Nparts, a1, b1, c1, dt, omega, axdir, axtilt, binary=False):
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
    alpha = np.radians(180. - omega * dt)

    phix = phi * np.sin(theta)
    phiy = phi * np.cos(theta)

    xtilt, ytilt, ztilt = bss.rotmatrix3d(xp, yp, zp, alpha, phiy, phix)


    r = np.sqrt(xtilt ** 2 + ytilt ** 2 + ztilt ** 2)

    xhat = xtilt / r
    yhat = ytilt / r
    zhat = ztilt / r


    # Rotate hat vector back to global frame
    theta = np.radians(axdir)
    phi = np.radians(axtilt)
    alpha = np.radians(omega * dt)

    phix = phi * np.sin(theta)
    phiy = phi * np.cos(theta)

    xhatp, yhatp, zhatp = bss.rotmatrix3d(xhat, yhat, zhat, alpha, phiy, phix)


    # calculate gravitational acceleration from shell potential
    C1 = (a1 ** 2) / (b1 ** 2)
    C2 = (a1 ** 2) / (c1 ** 2)

    apos = np.sqrt(xtilt ** 2 + C1 * ytilt ** 2 + C2 * ztilt ** 2)

    h = np.sqrt(a1**2 - b1**2)
    k = np.sqrt(a1**2 - c1**2)

    aellip = (G * Mtarg) / np.sqrt((apos**2 - h**2) * (apos**2 - k**2))

    axellip = aellip * xhatp
    ayellip = aellip * yhatp
    azellip = aellip * zhatp

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



def addnetgrav(sim, Mtarg, a1, b1, c1, Nparts, dt, omega, axdir, axtilt, binary=False):
    ''' Add net gravitational acceleration to simulation
    Sept 17, 2019 - in progress
    Oct 4, 2019 - testing '''

    asphere = spheregrav(sim, Mtarg, Nparts, binary)

    aellip  = ellipgrav(sim, Mtarg, Nparts, a1, b1, c1, dt, omega, axdir, axtilt, binary)

    axnet, aynet, aznet = rmaddgrav(aellip, asphere)


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
    dist = np.sqrt(xp**2 + yp**2)
    arot = omega ** 2 * dist
    rotxhat = xp / dist
    rotyhat = yp / dist
### Don't forget to tilt the xhat and yhat
    axrot = arot * rotxhat
    ayrot = arot * rotyhat

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





# Do not use the below functions.



def tiltrot(sim, per, Nplanets, Nparts, axtilt, axdir, M, a1, b1, c1, binary=False):

    dt = sim.dt

    omega = 360. / per

    rot = omega * dt

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
        bodyrange = np.arange(Nplanets, Nplanets + Nparts)

        x = np.asarray([sim.particles[int(j)].x for j in bodyrange])
        y = np.asarray([sim.particles[int(j)].y for j in bodyrange])
        z = np.asarray([sim.particles[int(j)].z for j in bodyrange])

        cx = sim.particles[0].x
        cy = sim.particles[0].y
        cz = sim.particles[0].z

    x -= cx
    y -= cy
    z -= cz

    xtilt = np.radians(axtilt) * np.sin(np.radians(axdir + 180))
    ytilt = np.radians(axtilt) * np.cos(np.radians(axdir + 180))

    xp = x*np.cos(rot)*np.cos(xtilt) \
         + y*(np.sin(rot)*np.cos(ytilt)-np.cos(rot)*np.sin(xtilt)*np.sin(ytilt)) \
         + z*(np.sin(rot)*np.cos(ytilt)+np.cos(rot)*np.sin(xtilt)*np.cos(ytilt))

    yp = -x*np.sin(rot)*np.cos(xtilt) \
         + y*(np.cos(rot)*np.cos(ytilt)+np.sin(rot)*np.sin(xtilt)*np.sin(ytilt)) \
         + z*(np.sin(rot)*np.sin(ytilt)+np.cos(rot)*np.sin(xtilt)*np.cos(ytilt))

    zp = -x*np.sin(xtilt) \
         - y*np.cos(xtilt)*np.sin(ytilt) \
         + z*np.cos(xtilt)*np.cos(ytilt)

    # Locate lat/lon:
    lon = np.degrees(np.arctan(yp / xp))
    lat = np.degrees(np.sin(zp / np.sqrt(xp**2 + yp**2 + zp**2)))

    # Point on the surface
    # xsurf = aT * np.cos(np.radians(lon)) * np.cos(np.radians(lat))
    # ysurf = bT * np.sin(np.radians(lon)) * np.cos(np.radians(lat))
    # zsurf = cT * np.sin(np.radians(lat))

    C1 = (a1 ** 2) / (b1 ** 2)
    C2 = (a1 ** 2) / (c1 ** 2)

    apos = np.sqrt(xp ** 2 + C1 * yp ** 2 + C2 * zp ** 2)

    h = np.sqrt(a1**2 - b1**2)
    k = np.sqrt(a1**2 - c1**2)

    pot = (sim.G * M) / np.sqrt((apos**2 - h**2) * (apos**2 - k**2))

    gx = -1*pot * (xp / apos)#np.sqrt(xp**2 + (a1 / b1)**2 * yp**2 + (a1 / c1) * zp**2)
    gy = -1*pot * (yp / apos)#np.sqrt(xp**2 + (a1 / b1)**2 * yp**2 + (a1 / c1) * zp**2)
    gz = -1*pot * (zp / apos)#np.sqrt(xp**2 + (a1 / b1)**2 * yp**2 + (a1 / c1) * zp**2)

    #gx = sim.G * M * (xsurf / apos) * (1 / np.sqrt((apos ** 2 - h ** 2) * (apos ** 2 - k ** 2)))
    #gy = sim.G * M * (ysurf / apos) * (1 / np.sqrt((apos ** 2 - h ** 2) * (apos ** 2 - k ** 2)))
    #gz = sim.G * M * (zsurf / apos) * (1 / np.sqrt((apos ** 2 - h ** 2) * (apos ** 2 - k ** 2)))

    gx += omega ** 2 * xp
    gy += omega ** 2 * yp
    #print out rotation phase


    # Rotate grav vector back to global coordinates
    gxtilt = np.radians(axtilt) * np.sin(np.radians(axdir))
    gytilt = np.radians(axtilt) * np.cos(np.radians(axdir))

    ggx = gx * np.cos(rot + 180.) * np.cos(gxtilt) \
         + gy * (np.sin(rot + 180.) * np.cos(gytilt) - np.cos(rot + 180.) * np.sin(gxtilt) * np.sin(gytilt)) \
         + gz * (np.sin(rot + 180.) * np.cos(gytilt) + np.cos(rot + 180.) * np.sin(gxtilt) * np.cos(gytilt))

    ggy = -gx * np.sin(rot + 180.) * np.cos(gxtilt) \
         + gy * (np.cos(rot + 180.) * np.cos(gytilt) + np.sin(rot + 180.) * np.sin(gxtilt) * np.sin(gytilt)) \
         + gz * (np.sin(rot + 180.) * np.sin(gytilt) + np.cos(rot + 180.) * np.sin(gxtilt) * np.cos(gytilt))

    ggz = -gx * np.sin(gxtilt) \
         - gy * np.cos(gxtilt) * np.sin(gytilt) \
         + gz * np.cos(gxtilt) * np.cos(gytilt)

    #ax += ggx
    #ay += ggy
    #az += ggz
    print('shapeggx', ggx.shape)
    return ggx, ggy, ggz





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
	sim     : REBOUND simulation
	m       : float; binary mass
	r       : float; binary radius
	a       : float; distance between binaries
	e       : float; (default=0) binary eccentricity
	i       : float; (default=0) binary inclination
	periap  : float; (default=0) argument of periapsis for binary
	ascnode : float; (default=0) longitude of ascending node for binary
	f       : float; (default=0) true anomaly of binary

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
    r = np.linspace(rmin, rmax, Res)
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
	sim     : REBOUND simulation
	Nparts  : int; number of particles
	Nplanets: int; number of planets
	r       : array; particle diameters
	rho     : float; particle density
	L       : float; defaults to luminosity of the sun
	msun    : float; defaults to mass of sun

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
    Tsun = 5.778e3        # Kelvin
    Rsun = 6.95508e8      # meters Sun radius

    tsi = 1.36e3          # W/m**2
    eps = 1.5
    c = 3e8               # m/s
    Msun = 2e30           #kg
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


