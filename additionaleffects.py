### Package of additional effects such as size distribution, radiation pressure, binaries, etc ###

import rebound
import numpy as np
import math as ma


def sizedist(Nparts, Res, rmin, rmax, p):
    '''Calculates power law size distribution for particles. Based on O'Brian (2003).

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
    C = Nparts * ( 1 / (rmax**(1-pow) - rmin**(1-pow)) )
    const = 1 - Nparts * ((rmax**(1-pow)) / (rmax**(1-pow) - rmin**(1-pow)))
    r = np.linspace(rmin, rmax, Res)
    N = -(C * r**(1-pow)) + const

    counts = np.round(N, 0)

    if np.sum(counts) > Nparts:
        counts[0] = counts[0] - (np.sum(counts) - Nparts)

    return r, counts, N


def partrad(r, counts, Nparts):
    '''Creates array of particle radii

    Parameters
    ----------
    r      : array; possible radii
    counts : array; number of particles for each diameter
    Nparts :   int; total number of particles

    Returns
    -------
    radii : array; particle radii
    '''

    if sum(counts) != Nparts:
        counts[0] = counts[0] + (Nparts - sum(counts))

    #posrad = D / 2.
    rad = []

    for i in range(len(counts)):
        numparts = int(counts[i])
        for n in range(numparts):
            rad.append(r[i])

    radii = np.asarray(rad)

    return radii


def partmass(radii, rho):
    '''Calculates mass of particles

    Parameters
    ----------
    radii : array; particle radii
    rho   : float; particle density

    Returns
    -------
    mass : array; particle masses
    '''

    mass = (4./3.) * np.pi * radii**3 * rho

    return mass



def poyntingrobertson(r, a, rho, L=3.828e26, msun=2e30):
    '''Calculates radiation pressure on particle sizes

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

    G = 6.67e-11 * km**3 * days**2

    Lum = L * days**3 * km**-2

    c = 3e8 * (days/km)

    fpr = ((Lum * r**2) / (4 * c**2)) * np.sqrt((G * msun) / a**5)

    acc = fpr / ( (4./3.)* np.pi * r**3 * rho )

    return acc



def tidalforce(M, m, d, r):
    '''Calculation of tidal force

    Parameters
    ----------
    M : float; mass of perturbing body
    m : float; mass of perturbed body
    d : float; distance between bodies
    r : float; radius of perturbing body

    Returns
    -------
    F : float; tidal force
    '''

    G = 6.67e-11   # kg m^-2 s^-3

    F = (2 * G * M * m * r) / (d**3)

    return F



def binary(sim, m, r, a, e=0, i=0, periap=0, ascnode=0, f=0):
    '''Adds binary component to system (orbits central body)

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


def radforce(sim, Nparts, Nplanets, r, rho, L=3.828e26, msun=2e30):
    '''Calculates radiation pressure on particle sizes

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

    dist = np.sqrt(x**2 + y**2 + z**2)
    theta = np.arctan(y / x)
    phi   = np.arccos(z / dist)

    acc = poyntingrobertson(r, dist, rho, L, msun)

    return acc




def nonsymgrav(sim, semi, pos, rho):
    '''
    *** Do not use this function ***
    Calculates non-axisymmetric gravitational potential for ellipsoidal target body

    Parameters
    ----------
    sim  : REBOUND simulation
    semi : tuple; semimajor axes of ellipsoid
    pos  : tuple; position coordinates (Cartesian)
    rho  : float; density of target body

    Returns
    -------
    gravpot : float; gravitational potential at point given by pos
    '''
    G = sim.G

    a1 = semi[0]
    a2 = semi[1]
    a3 = semi[2]

    x1 = pos[0]
    x2 = pos[1]
    x3 = pos[2]

    a = x1**2 + x2**2 + x3**2 - a1**2 - a2**2 - a3**2
    b = -x1**2*a2**2 - x1**2*a3**2 - x2**2*a1**2 - x2**2*a3**2 - x3**2*a1**2 - x3**2*a2**2 + a1**2*a3**2 + a2**2*a3**2 + a1**2*a2**2
    c = x1**2*a2**2*a3**2 + x2**2*a1**2*a3**2 + x3**2*a1**2*a2**2 - a1**2*a2**2*a3**2

    root1 = (c**(1./3.)) / (-a * b - 1)**(1./3.)
    root2 = -((-c)**(1./3.)) / (-a * b - 1)**(1./3.)
    root3 = (-1)**(2./3.) * (c**(1./3.)) / (-a * b - 1)**(1./3.)

    roots = []
    if root1 > 0:
        roots.append(root1)
    if root2 > 0:
        roots.append(root2)
    if root3 > 0:
        roots.append(root3)

    maxroot = np.max(roots)

    delta = np.sqrt((a1**2 + maxroot) * (a2**2 + maxroot) * (a3**2 + maxroot))

    I = a1 * a2 * a3 * np.integrate(1. / delta, maxroot, np.inf)

    A1 = a1 * a2 * a3 * np.integrate(1. / (delta * (a1**2 + maxroot)), maxroot, np.inf)
    A2 = a1 * a2 * a3 * np.integrate(1. / (delta * (a2**2 + maxroot)), maxroot, np.inf)
    A3 = a1 * a2 * a3 * np.integrate(1. / (delta * (a3**2 + maxroot)), maxroot, np.inf)

    sumA = A1 * x1**2 + A2 * x2**2 + A3 * x3**2

    #gravpot = np.pi * G * rho * (I - sumA)




    gravpot =  (G * mtarg / (np.sqrt((x-x0)**2 + (y-y0)**2 + (z-z0)**2)))

    C   = []
    S   = []
    Pmn = []

    for n in range(0, 3):
        for m in range(0, n+1):

            if n == 0:
                if m == 0:
                    P = 1

            if n == 1:
                if m == 0:
                    P = np.sin(theta)

                if m == 1:
                    P = -np.sqrt(1 - (np.sin(theta))**2)

            if n == 2:
                if m == 0:
                    P = 0.5 * (3 * (np.sin(theta))**2 - 1)

                if m == 1:
                    P = -3 * np.sin(theta) * np.sqrt(1 - np.sin(theta)**2)

                if m == 2:
                    P = 3 * (1 - (np.sin(theta))**2)


            if m == n:
                delta = 1
            else:
                delta = 0

            stokesC = (2 - delta) * ( ma.factorial(n-m) / ma.factorial(n+m) ) * P

            C.append(stokesC * np.cos(m * theta))
            S.append(stokesC * np.sin(m * theta))
            Pmn.append(P)

    C   = np.asarray(C)
    S   = np.asarray(S)
    Pmn = np.asarray(Pmn)

    gravpot *= np.sum(Pmn * (C + S))

    return gravpot




def cosgamma(theta, phi, thetaprime, phiprime):
    cgamma = np.cos(theta) * np.cos(thetaprime) + np.sin(theta) * np.sin(thetaprime) * np.cos(phi - phiprime)
    return cgamma

def legendre(x, n):
    Pn =  ( x * n * (n+1) / (2**(n-1) * ma.factorial(n+1)) ) * (x**2 - 1)**(n-1)
    return Pn

def flattening(a, b, c, phiprime, thetaprime):
    fab = (a - b) / a

    fac = (a - c) / a

    fbc = (b - c) / b

    if -45 < phiprime < 45:
        if -45 < thetaprime < -135 or 45 < thetaprime < 135:
            return fab

        if -45 < thetaprime < 45 or -135 < thetaprime < 135:
            return fac

    if 45 < phiprime < 90 or -45 < phiprime < -90:
        return fbc

def ellipticity(a, b, c, phiprime, thetaprime):
    f = flattening(a, b, c, phiprime, thetaprime)

    e2 = 2 * f / (f - 3)

    return e2

def position(ar, a, b, c, theta, phi, thetaprime, phiprime, lmax):
    polysum = []

    cgamma = cosgamma(theta, phi, thetaprime, phiprime)

    for l in range(0, lmax):
        polysum.append(ellipticity(a, b, c, phiprime, thetaprime) * legendre(cgamma, l))

    r = ar * (1 + np.sum(polysum))

    return r

def rthetphipart(sim, Nplanets):
    x = []
    y = []
    z = []

    parts = sim.particles[:]
    for i in np.arange(len(parts)):
        x.append(parts[i].x)
        y.append(parts[i].y)
        z.append(parts[i].z)

    x = np.asarray(x)
    y = np.asarray(y)
    z = np.asarray(z)
    rprime = np.sqrt(x[Nplanets:]**2 + y[Nplanets:]**2 + z[Nplanets:]**2)

    thetaprime = np.arctan(y[Nplanets:] / x[Nplanets])

    phiprime = np.arctan(np.sqrt(x[Nplanets:]**2 + y[Nplanets:]**2) / z[Nplanets:])

    return rprime, thetaprime, phiprime

def gravpot(sim, theta, phi, thetaprime, phiprime, lmax, a, ap, rho):

    #G = sim.G
    #cgamma = cosgamma(theta, phi, thetaprime, phiprime)
    #r = position(a, theta, phi, lmax)
    #rprime = position(ap, thetaprime, phiprime, lmax)

    #el = ellipticity(a, b, c, phiprime, thetaprime)

    #summation = []
    #for l in range(0, lmax):
     #   Pl = legendre(cgamma, l)
     #   summation.append((3./(2 * l + 1)) * ((a - ap)**(l+3) / (r - rprime)**(l+1)) * el * Pl)

    #U = -(4. * np.pi / 3.) * G * rho * (((a - ap)**3 / (r - rprime)) + np.sum(summation))

    return U




def gravaccel(sim, a, b, c, theta, phi, thetaprime, phiprime, rprime, ar, rho, lmax, Nparts):
    '''Calculates gravitational acceleration of particles due to gravitational potential

    Parameters
    ----------

    Returns
    -------
    gx : float; x-component of gravitational acceleration
    gy : float; y-component of gravitational acceleration
    gz : float; z-component of gravitational acceleration
    '''

    gx = []
    gy = []
    gz = []

    for i in np.arange(Nparts):
        r = position(ar, a, b, c, theta, phi, thetaprime[i], phiprime[i], lmax)

        dcosdthet = (-np.cos(theta)*np.sin(thetaprime[i]) + np.sin(theta)*np.cos(thetaprime[i])*np.cos(phi-phiprime[i]))
        dcosdphi  = np.sin(theta)*np.sin(thetaprime[i])*np.sin(phi-phiprime[i])

        e = ellipticity(a, b, c, phiprime[i], thetaprime[i])
        cgamma = cosgamma(theta, phi, thetaprime[i], phiprime[i])
        summation = []
        for l in range(0, lmax):
            Pl  = legendre(cgamma, l)
            Pl1 = legendre(cgamma, l-1)
            summation.append(e * (l / (cgamma**2 - 1)) * (cgamma*Pl - Pl1))
        drdcos = ar * np.sum(summation)

        drdthet = dcosdthet * drdcos
        drdphi  = dcosdphi  * drdcos

        G = sim.G
        summation2 = []
        for l in range(0, lmax):
            Pl = legendre(cgamma, l)
            summation2.append((3./(2*l+1)) * (-l-3) * (a**(l+3) / (r-rprime[i])**(l+2)) * e * Pl)

        dUdr = -(4.*np.pi / 3.) * G * rho * (-(2*ar**3) / (r-rprime[i])**2 + np.sum(summation2))

        dUdrp   = -dUdr
        dUdthep = (1./rprime[i]) * dUdr * drdthet
        dUdphip = (1./(rprime[i]*np.sin(thetaprime[i]))) * dUdr * drdphi

    gx.append(dUdrp * np.sin(dUdthep) * np.cos(dUdphip))
    gy.append(dUdrp * np.sin(dUdthep) * np.sin(dUdphip))
    gz.append(dUdrp * np.cos(dUdthep))

    return np.asarray(gx), np.asarray(gy), np.asarray(gz)














def nonsymremove(sim, semi, outer, Nplanets, landed=0, gone=0):
    '''Create nonspherical body and remove particles on surface or outside system bounds

    Parameters
    ----------
    sim      : REBOUND simuluation
    semi     : tuple; semimajor axes
    outer    : float; outer limit that particles can go
    Nplanets :   int; number of planets
    landed   :   int; (default=0) number of particles that landed
    gone     :   int; (default=0) number of particles that left the simulation

    Returns
    -------
    sim    : REBOUND simulation
    rminds : list; indices of removed particles
    landed :  int; number of landed particles
    gone   :  int; number of particles that exited simulation
    '''

    partind = Nplanets
    rminds = []

    for p in sim.particles[Nplanets:]:

        x = p.x
        y = p.y
        z = p.z

        if (x**2 / semi[0]**2) + (y**2 / semi[1]**2) + (z**2 / semi[2]**2) <= 1:
            sim.remove(partind)
            rminds.append(partind)
            landed += 1

        if x**2 + y**2 + z**2 > outer:
            sim.remove(partind)
            rminds.append(partind)
            gone += 1

        partind += 1

    return sim, rminds, landed, gone
