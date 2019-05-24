### Package of additional effects such as size distribution, radiation pressure, binaries, etc ###

import rebound
import numpy as np
import math as ma



### Ellipsoidal Gravitational Acceleration ###

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
        x = np.asarray([sim.particles[int(j)].x for j in range(1, Nplanets + Nparts)])
        y = np.asarray([sim.particles[int(j)].y for j in range(1, Nplanets + Nparts)])
        z = np.asarray([sim.particles[int(j)].z for j in range(1, Nplanets + Nparts)])

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
    ggx = gmag * np.cos(np.radians(phi)) * np.sin(np.radians(theta))
    ggy = gmag * np.cos(np.radians(phi)) * np.cos(np.radians(theta))
    ggz = gmag * np.sin(np.radians(phi))

    return ggx, ggy, ggz


def gravity(sim, M, Nplanets, Nparts, atarg, btarg, ctarg, omega, tilt, axdir, time, binary=False):
    '''In progress
	Inputs gravitational acceleration vectors
	Applies acceleration vectors to system objects'''

    if binary == True:
        bodies = np.arange(0, Nplanets + Nparts)
        exclude = [1]
        bodyrange = np.delete(bodies, exclude)

        ax = np.asarray([sim.particles[int(j)].ax for j in bodyrange])
        ay = np.asarray([sim.particles[int(j)].ay for j in bodyrange])
        az = np.asarray([sim.particles[int(j)].az for j in bodyrange])

    else:
        ax = np.asarray([sim.particles[int(j)].ax for j in range(1, Nplanets + Nparts)])
        ay = np.asarray([sim.particles[int(j)].ay for j in range(1, Nplanets + Nparts)])
        az = np.asarray([sim.particles[int(j)].az for j in range(1, Nplanets + Nparts)])

    gx, gy, gz = localvector(sim, M, Nplanets, Nparts, atarg, btarg, ctarg, omega, tilt, axdir, time, binary)

    ggx, ggy, ggz = globalvector(sim, tilt, gx, gy, gz, Nplanets, Nparts, binary)

    ax += ggx
    ay += ggy
    az += ggz

    return sim




def localvector(sim, M, Nplanets, Nparts, atarg, btarg, ctarg, omega, tilt, axdir, time, binary=False):
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

    dist = np.sqrt((x - x0) ** 2 + (y - y0) ** 2 + (z - z0) ** 2)

    dx = dist * (np.sin(np.radians(tilt)) * np.cos(np.radians(axdir)) - np.cos(np.radians(omega * time)))
    dy = dist * (np.sin(np.radians(tilt)) * np.sin(np.radians(axdir)) - np.sin(np.radians(omega * time)))
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
    print (len(radii))
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
    ''' ***This is ok!***
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

    km = 1e3  # m in km
    days = 86400.  # seconds in day

    G = sim.G  # 6.67e-11 * km**3 * days**2

    Lum = L * days ** 3 * km ** -2

    c = 3e8 * (days / km)

    fpr = ((Lum * r ** 2) / (4 * c ** 2)) * np.sqrt((G * msun) / a ** 5)

    acc = fpr / ((4. / 3.) * np.pi * r ** 3 * rho)

    return acc

def radforce(sim, Nparts, Nplanets, r, rho, L=3.828e26, msun=2e30):
    ''' ***This is ok!***
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



def solarradpress(sim, Nplanets, type2=False, rho=2000):
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

    if type2 == True:
        sunind = 0
    else:
        sunind = Nplanets-1


    sun = sim.particles[int(sunind)]
    sunx = sun.x
    suny = sun.y
    sunz = sun.z

    for p in sim.particles[Nplanets:]:
        dist = np.sqrt((p.x-sunx)**2 + (p.y-suny)**2 + (p.z-sunz)**2)

        #Fs = SBconst * Tsun**4 * (Rsun**2 / dist**2)

        beta = (3 * eps * tsi * dist**2) / (4 * c * G * Msun * rho * p.r)
        amean = ((G * Msun) / (dist**2)) * (1 - beta)
        #amean = ((G * Msun) / (dist**2)) - (3 * eps * tsi) / (4 * c * rho * p.r)
        ax.append(amean * ((p.x-sunx) / dist))
        ay.append(amean * ((p.y-suny) / dist))
        az.append(amean * ((p.z-sunz) / dist))

    acc = (np.asarray(ax), np.asarray(ay), np.asarray(az))

    return acc


